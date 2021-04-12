#include <Eigen/Core>
#include <Eigen/LU>
#include <elasty/alembic-manager.hpp>
#include <elasty/fem.hpp>
#include <mathtoolbox/l-bfgs.hpp>
#include <timer.hpp>
#include <vector>

namespace
{
    constexpr std::size_t k_num_dims = 2;

    constexpr double k_youngs_modulus = 200.0;
    constexpr double k_poisson_ratio  = 0.45;

    constexpr double k_first_lame  = elasty::fem::calcFirstLame(k_youngs_modulus, k_poisson_ratio);
    constexpr double k_second_lame = elasty::fem::calcSecondLame(k_youngs_modulus, k_poisson_ratio);

    constexpr unsigned k_num_substeps = 5;
    constexpr double   k_delta_time   = 1.0 / 60.0;

    constexpr double k_damping_factor = 0.0;

    constexpr double k_spring_stiffness = 10000.0;

    enum class Model
    {
        CoRotational,
        StVenantKirchhoff
    };

    constexpr Model k_model = Model::CoRotational;
} // namespace

template <int N> struct Mesh
{
    using ElemList = Eigen::Matrix<std::int32_t, N + 1, Eigen::Dynamic>;

    ElemList elems;

    Eigen::VectorXd x_rest;
    Eigen::VectorXd x;
    Eigen::VectorXd v;
    Eigen::VectorXd f;

    double mass = 1.0;

    /// \details Should be precomputed
    Eigen::VectorXd lumped_mass;

    /// \details Should be precomputed
    std::vector<double> volume_array;

    /// \details Should be precomputed
    std::vector<Eigen::Matrix<double, N, N>> rest_shape_mat_inv_array;

    /// \details Should be precomputed
    std::vector<Eigen::Matrix<double, N * N, N*(N + 1)>> vec_PFPx_array;
};

using TriangleMesh = Mesh<2>;

struct Constraint
{
    std::size_t                              vert_index;
    std::function<Eigen::Vector2d(double t)> motion;
    double                                   stiffness;
};

template <typename Derived> typename Derived::Scalar calcEnergyDensity(const Eigen::MatrixBase<Derived>& deform_grad)
{
    switch (k_model)
    {
        case Model::StVenantKirchhoff:
            return elasty::fem::calcStVenantKirchhoffEnergyDensity(deform_grad, k_first_lame, k_second_lame);
        case Model::CoRotational:
            return elasty::fem::calcCoRotationalEnergyDensity(deform_grad, k_first_lame, k_second_lame);
    }
}

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, 2, 2> calcPiolaStress(const Eigen::MatrixBase<Derived>& deform_grad)
{
    switch (k_model)
    {
        case Model::StVenantKirchhoff:
            return elasty::fem::calcStVenantKirchhoffPiolaStress(deform_grad, k_first_lame, k_second_lame);
        case Model::CoRotational:
            return elasty::fem::calcCoRotationalPiolaStress(deform_grad, k_first_lame, k_second_lame);
    }
}

class VariationalImplicit2dEngine
{
public:
    VariationalImplicit2dEngine() {}

    void proceedFrame()
    {
        const std::size_t num_verts = m_mesh.x_rest.size() / k_num_dims;

        // Reset forces
        m_mesh.f = Eigen::VectorXd::Zero(2 * num_verts);

        // Apply gravity force
        for (std::size_t i = 0; i < num_verts; ++i)
        {
            m_mesh.f[i * 2 + 1] += m_mesh.lumped_mass(i * 2 + 1) * (-9.80665);
        }

        // Calculate the inverse lumped mass matrix
        const auto W = m_mesh.lumped_mass.cwiseInverse().asDiagonal();

        // Calculate the "inertia" position
        const double&         h = m_delta_physics_time;
        const Eigen::VectorXd y = m_mesh.x + h * m_mesh.v + h * h * W * m_mesh.f;

        const auto calcInternalPotential = [&](const Eigen::VectorXd& x) {
            double sum = 0.0;

            // Elastic potential
            for (std::size_t i = 0; i < m_mesh.elems.cols(); ++i)
            {
                const auto& indices = m_mesh.elems.col(i);

                // Retrieve precomputed values
                const auto& D_m_inv = m_mesh.rest_shape_mat_inv_array[i];
                const auto& area    = m_mesh.volume_array[i];

                // Calculate the deformation gradient $\mathbf{F}$
                const auto F = elasty::fem::calc2dTriangleDeformGrad(
                    x.segment<2>(2 * indices[0]), x.segment<2>(2 * indices[1]), x.segment<2>(2 * indices[2]), D_m_inv);

                sum += area * calcEnergyDensity(F);
            }

            // Attach-spring potential
            for (const auto& constraint : m_constraints)
            {
                const std::size_t vert_index   = constraint.vert_index;
                const double&     k            = constraint.stiffness;
                const auto        p            = x.segment<2>(vert_index * 2);
                const auto        q            = constraint.motion(m_physics_time + m_delta_physics_time);
                const double      squared_dist = (p - q).squaredNorm();

                sum += 0.5 * k * squared_dist;
            }

            return sum;
        };

        const auto calcInternalPotentialGrad = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::VectorXd sum = Eigen::VectorXd::Zero(x.size());
            for (std::size_t i = 0; i < m_mesh.elems.cols(); ++i)
            {
                const auto& indices = m_mesh.elems.col(i);

                // Retrieve precomputed values
                const auto& D_m_inv  = m_mesh.rest_shape_mat_inv_array[i];
                const auto& area     = m_mesh.volume_array[i];
                const auto& vec_PFPx = m_mesh.vec_PFPx_array[i];

                // Calculate the deformation gradient $\mathbf{F}$
                const auto F = elasty::fem::calc2dTriangleDeformGrad(
                    x.segment<2>(2 * indices[0]), x.segment<2>(2 * indices[1]), x.segment<2>(2 * indices[2]), D_m_inv);

                // Calculate $\frac{\partial \Phi}{\partial \mathbf{x}}$ and related values
                const auto P      = calcPiolaStress(F);
                const auto vec_P  = Eigen::Map<const Eigen::Vector4d>(P.data(), P.size());
                const auto PPsiPx = vec_PFPx.transpose() * vec_P;

                // Calculate $\frac{\partial E}{\partial \mathbf{x}}$
                const auto PEPx = area * PPsiPx;

                sum.segment<2>(2 * indices[0]) += PEPx.segment<2>(0 * 2);
                sum.segment<2>(2 * indices[1]) += PEPx.segment<2>(1 * 2);
                sum.segment<2>(2 * indices[2]) += PEPx.segment<2>(2 * 2);
            }

            for (const auto& constraint : m_constraints)
            {
                const std::size_t vert_index = constraint.vert_index;
                const double&     k          = constraint.stiffness;
                const auto        p          = x.segment<2>(vert_index * 2);
                const auto        q          = constraint.motion(m_physics_time + m_delta_physics_time);
                const auto        r          = p - q;

                sum.segment<2>(2 * vert_index) += k * r;
            }

            return sum;
        };

        const auto calcMomentumPotential = [&](const Eigen::VectorXd& x) {
            return (0.5 / (h * h)) * (x - y).transpose() * m_mesh.lumped_mass.asDiagonal() * (x - y);
        };

        const auto calcMomentumPotentialGrad = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return (1.0 / (h * h)) * m_mesh.lumped_mass.asDiagonal() * (x - y);
        };

        const auto calcObjective = [&](const Eigen::VectorXd& x) {
            const double momentum_potential = calcMomentumPotential(x);
            const double internal_potential = calcInternalPotential(x);

            assert(momentum_potential >= 0.0);
            assert(internal_potential >= 0.0);

            return momentum_potential + internal_potential;
        };

        const auto calcObjectiveGrad = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return calcMomentumPotentialGrad(x) + calcInternalPotentialGrad(x);
        };

        // Solve the minimization problem
        unsigned        num_iters;
        Eigen::VectorXd x_opt;
        mathtoolbox::optimization::RunLBfgs(y, calcObjective, calcObjectiveGrad, 1e-06, 100, x_opt, num_iters);

        // Update the internal state
        m_mesh.v = (1.0 / h) * (x_opt - m_mesh.x);
        m_mesh.x = x_opt;

        // Apply naive damping
        m_mesh.v *= std::exp(-k_damping_factor * m_delta_physics_time);

        // Update time counter
        m_physics_time += m_delta_physics_time;
    }

    void initializeScene()
    {
        constexpr std::size_t num_cols  = 20;
        constexpr std::size_t num_rows  = 5;
        constexpr std::size_t num_verts = (num_cols + 1) * (num_rows + 1);
        constexpr std::size_t num_elems = (num_cols * num_rows) * 2;
        constexpr double      size      = 1.0;

        m_mesh.elems.resize(3, num_elems);
        m_mesh.x_rest.resize(num_verts * k_num_dims);

        // Generate a triangle mesh
        for (std::size_t col = 0; col < num_cols; ++col)
        {
            for (std::size_t row = 0; row < num_rows; ++row)
            {
                const auto base = col * (num_rows + 1) + row;

                m_mesh.elems.col(2 * num_rows * col + 2 * row + 0) << 0 + base, 1 + base, (num_rows + 1) + 1 + base;
                m_mesh.elems.col(2 * num_rows * col + 2 * row + 1) << 0 + base, (num_rows + 1) + 1 + base,
                    (num_rows + 1) + base;

                m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * col + row), k_num_dims) =
                    Eigen::Vector2d{col * 1.0, -1.0 * row};
            }
            m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * col + num_rows), k_num_dims) =
                Eigen::Vector2d{col * 1.0, -1.0 * num_rows};
        }
        for (std::size_t row = 0; row < num_rows; ++row)
        {
            m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * num_cols + row), k_num_dims) =
                Eigen::Vector2d{num_cols * 1.0, -1.0 * row};
        }
        m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * num_cols + num_rows), k_num_dims) =
            Eigen::Vector2d{num_cols * 1.0, -1.0 * num_rows};

        // Set transform
        m_mesh.x_rest *= 1.0 / static_cast<double>(num_rows);
        for (std::size_t vert = 0; vert < num_verts; ++vert)
        {
            m_mesh.x_rest[2 * vert + 1] += 0.5;
        }
        m_mesh.x_rest *= size;

        // Initialize other values
        m_mesh.x = m_mesh.x_rest;
        m_mesh.v = Eigen::VectorXd::Zero(k_num_dims * num_verts);
        m_mesh.f = Eigen::VectorXd::Zero(k_num_dims * num_verts);

        // Set constraints
        for (std::size_t i = 0; i < num_rows + 1; ++i)
        {
            const auto motion = [&, i](double) -> Eigen::Vector2d { return m_mesh.x_rest.segment<2>(i * 2); };

            m_constraints.push_back(Constraint{i, motion, k_spring_stiffness});
        }
        for (std::size_t i = (num_rows + 1) * num_cols; i < (num_rows + 1) * (num_cols + 1); ++i)
        {
            const auto ease = [](double x) { return -(std::cos(3.14159265358979 * x) - 1.0) * 0.5; };

            const auto motion = [&, i](double t) -> Eigen::Vector2d {
                const auto   x_init = m_mesh.x_rest.segment<2>(i * 2);
                const auto   dir    = Eigen::Vector2d{3.0, 0.0};
                const double t_0    = 0.8;
                const double t_1    = t_0 + 1.0;
                const double a      = (t < t_0) ? 0.0 : ((t < t_1) ? t - t_0 : t_1 - t_0);
                const double b      = ease(a);

                return x_init + b * dir;
            };

            m_constraints.push_back(Constraint{i, motion, k_spring_stiffness});
        }

        // Perform precomputation
        m_mesh.volume_array.resize(m_mesh.elems.cols());
        m_mesh.rest_shape_mat_inv_array.resize(m_mesh.elems.cols());
        m_mesh.vec_PFPx_array.resize(m_mesh.elems.cols());
        for (std::size_t elem_index = 0; elem_index < m_mesh.elems.cols(); ++elem_index)
        {
            const auto& indices = m_mesh.elems.col(elem_index);

            m_mesh.volume_array[elem_index] = elasty::fem::calc2dTriangleArea(m_mesh.x_rest.segment<2>(2 * indices[0]),
                                                                              m_mesh.x_rest.segment<2>(2 * indices[1]),
                                                                              m_mesh.x_rest.segment<2>(2 * indices[2]));
            m_mesh.rest_shape_mat_inv_array[elem_index] =
                elasty::fem::calc2dShapeMatrix(m_mesh.x_rest.segment<2>(2 * indices[0]),
                                               m_mesh.x_rest.segment<2>(2 * indices[1]),
                                               m_mesh.x_rest.segment<2>(2 * indices[2]))
                    .inverse();
            m_mesh.vec_PFPx_array[elem_index] =
                elasty::fem::calcVecTrianglePartDeformGradPartPos(m_mesh.rest_shape_mat_inv_array[elem_index]);
        }
        m_mesh.lumped_mass = elasty::fem::calcLumpedMasses(m_mesh.x_rest, m_mesh.elems, m_mesh.mass);
    }

    /// \brief Getter of the delta physics time.
    ///
    /// \details The value equals to the delta frame time devided by the number of substeps.
    double getDeltaPhysicsTime() const { return m_delta_physics_time; }

    void setDeltaPhysicsTime(const double delta_physics_time) { m_delta_physics_time = delta_physics_time; }

    const TriangleMesh* getMesh() const { return &m_mesh; }

private:
    double m_delta_physics_time = 1.0 / 60.0;
    double m_physics_time       = 0.0;

    std::vector<Constraint> m_constraints;

    TriangleMesh m_mesh;
};

int main(int argc, char** argv)
{
    VariationalImplicit2dEngine engine;

    engine.initializeScene();
    engine.setDeltaPhysicsTime(k_delta_time / static_cast<double>(k_num_substeps));

    const auto        mesh      = engine.getMesh();
    const std::size_t num_verts = mesh->x_rest.size() / 2;
    const std::size_t num_elems = mesh->elems.cols();

    auto alembic_manager = elasty::createTriangleMesh2dAlembicManager(
        "./out.abc", k_delta_time, num_verts, num_elems, mesh->x.data(), mesh->elems.data());

    for (unsigned int frame = 0; frame < 240; ++frame)
    {
        timer::Timer t(std::to_string(frame));

        alembic_manager->submitCurrentStatus();

        for (int i = 0; i < k_num_substeps; ++i)
        {
            engine.proceedFrame();
        }
    }

    return 0;
}
