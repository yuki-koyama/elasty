#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <elasty/alembic-manager.hpp>
#include <elasty/fem.hpp>
#include <mathtoolbox/l-bfgs.hpp>
#include <timer.hpp>
#include <vector>

namespace
{
    constexpr double k_youngs_modulus = 600.0;
    constexpr double k_poisson_ratio  = 0.45;

    constexpr double k_first_lame  = elasty::fem::calcFirstLame(k_youngs_modulus, k_poisson_ratio);
    constexpr double k_second_lame = elasty::fem::calcSecondLame(k_youngs_modulus, k_poisson_ratio);

    constexpr unsigned k_num_substeps = 2;
    constexpr double   k_delta_time   = 1.0 / 60.0;

    constexpr double k_damping_factor = 0.0;

    constexpr double k_spring_stiffness    = 1000.0;
    constexpr double k_collision_stiffness = 100.0;

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

using TetraMesh = Mesh<3>;

struct Constraint
{
    std::size_t                              vert_index;
    std::function<Eigen::Vector3d(double t)> motion;
    double                                   stiffness;
};

struct Collision
{
    std::size_t     vert_index;
    Eigen::Vector3d surface_pos;
    Eigen::Vector3d surface_normal;
};

class HalfSpaceCollider
{
public:
    HalfSpaceCollider(const Eigen::Vector3d& representative_point, const Eigen::Vector3d& normal)
        : m_representative_point(representative_point), m_normal(normal.normalized())
    {
        assert(!normal.isZero());
    }

    bool testCollision(const Eigen::Vector3d& point, const double tol) const
    {
        return (point - m_representative_point).dot(m_normal) < tol;
    }

    void retrieveCollisionInfo(const Eigen::Vector3d& point,
                               Eigen::Vector3d&       surface_pos,
                               Eigen::Vector3d&       surface_normal) const
    {
        surface_pos    = point - (point - m_representative_point).dot(m_normal) * m_normal;
        surface_normal = m_normal;
    }

private:
    const Eigen::Vector3d m_representative_point;
    const Eigen::Vector3d m_normal;
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
Eigen::Matrix<typename Derived::Scalar, 3, 3> calcPiolaStress(const Eigen::MatrixBase<Derived>& deform_grad)
{
    switch (k_model)
    {
        case Model::StVenantKirchhoff:
            return elasty::fem::calcStVenantKirchhoffPiolaStress(deform_grad, k_first_lame, k_second_lame);
        case Model::CoRotational:
            return elasty::fem::calcCoRotationalPiolaStress(deform_grad, k_first_lame, k_second_lame);
    }
}

class VariationalImplicit3dEngine
{
public:
    VariationalImplicit3dEngine() {}

    void proceedFrame()
    {
        const std::size_t num_verts = m_mesh.x_rest.size() / 3;
        const std::size_t num_elems = m_mesh.elems.cols();

        // Reset forces
        m_mesh.f = Eigen::VectorXd::Zero(3 * num_verts);

        // Apply gravity force
        for (std::size_t i = 0; i < num_verts; ++i)
        {
            m_mesh.f[i * 3 + 1] += m_mesh.lumped_mass(i * 3 + 1) * (-9.80665);
        }

        // Calculate the inverse lumped mass matrix
        const auto W = m_mesh.lumped_mass.cwiseInverse().asDiagonal();

        // Calculate the "inertia" position
        const double&         h = m_delta_physics_time;
        const Eigen::VectorXd y = m_mesh.x + h * m_mesh.v + h * h * W * m_mesh.f;

        // Define colliders
        // TODO: Improve the scene setup
        const HalfSpaceCollider floor{Eigen::Vector3d{0.0, -1.0, 0.0}, Eigen::Vector3d{0.0, 1.0, 0.0}};

        // Detect collisions
        // TODO: Improve the collision detection algorithtm
        std::vector<Collision> collisions;
        for (std::size_t vert_index = 0; vert_index < num_verts; ++vert_index)
        {
            const Eigen::Vector3d predicted_pos = y.segment<3>(vert_index * 3);

            if (floor.testCollision(predicted_pos, 0.0))
            {
                Eigen::Vector3d surface_pos;
                Eigen::Vector3d surface_normal;
                floor.retrieveCollisionInfo(predicted_pos, surface_pos, surface_normal);

                const Collision collision{vert_index, surface_pos, surface_normal};
                collisions.push_back(collision);
            }
        }

        const auto calcInternalPotential = [&](const Eigen::VectorXd& x) {
            double sum = 0.0;

            // Elastic potential
            for (std::size_t elem_index = 0; elem_index < num_elems; ++elem_index)
            {
                const auto& indices = m_mesh.elems.col(elem_index);

                // Retrieve precomputed values
                const auto& D_m_inv = m_mesh.rest_shape_mat_inv_array[elem_index];
                const auto& vol     = m_mesh.volume_array[elem_index];

                // Calculate the deformation gradient $\mathbf{F}$
                const auto F = elasty::fem::calcTetrahedronDeformGrad(x.segment<3>(3 * indices[0]),
                                                                      x.segment<3>(3 * indices[1]),
                                                                      x.segment<3>(3 * indices[2]),
                                                                      x.segment<3>(3 * indices[3]),
                                                                      D_m_inv);

                sum += vol * calcEnergyDensity(F);
            }

            // Attach-spring potential
            for (const auto& constraint : m_constraints)
            {
                const std::size_t vert_index   = constraint.vert_index;
                const double&     k            = constraint.stiffness;
                const auto        p            = x.segment<3>(vert_index * 3);
                const auto        q            = constraint.motion(m_physics_time + m_delta_physics_time);
                const double      squared_dist = (p - q).squaredNorm();

                sum += 0.5 * k * squared_dist;
            }

            // Collision energy
            for (const auto& collision : collisions)
            {
                const auto&  vert_index = collision.vert_index;
                const auto   pos        = x.segment<3>(vert_index * 3);
                const auto&  normal     = collision.surface_normal;
                const double a          = (collision.surface_pos - pos).transpose() * normal;

                if (a > 0.0)
                {
                    sum += k_collision_stiffness * a * a;
                }
            }

            return sum;
        };

        const auto calcInternalPotentialGrad = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::VectorXd sum = Eigen::VectorXd::Zero(x.size());
            for (std::size_t elem_index = 0; elem_index < num_elems; ++elem_index)
            {
                const auto& indices = m_mesh.elems.col(elem_index);

                // Retrieve precomputed values
                const auto& D_m_inv  = m_mesh.rest_shape_mat_inv_array[elem_index];
                const auto& vol      = m_mesh.volume_array[elem_index];
                const auto& vec_PFPx = m_mesh.vec_PFPx_array[elem_index];

                // Calculate the deformation gradient $\mathbf{F}$
                const auto F = elasty::fem::calcTetrahedronDeformGrad(x.segment<3>(3 * indices[0]),
                                                                      x.segment<3>(3 * indices[1]),
                                                                      x.segment<3>(3 * indices[2]),
                                                                      x.segment<3>(3 * indices[3]),
                                                                      D_m_inv);

                // Calculate $\frac{\partial \Phi}{\partial \mathbf{x}}$ and related values
                const auto P      = calcPiolaStress(F);
                const auto vec_P  = Eigen::Map<const Eigen::Matrix<double, 9, 1>>(P.data(), P.size());
                const auto PPsiPx = vec_PFPx.transpose() * vec_P;

                // Calculate $\frac{\partial E}{\partial \mathbf{x}}$
                const auto PEPx = vol * PPsiPx;

                sum.segment<3>(3 * indices[0]) += PEPx.segment<3>(0 * 3);
                sum.segment<3>(3 * indices[1]) += PEPx.segment<3>(1 * 3);
                sum.segment<3>(3 * indices[2]) += PEPx.segment<3>(2 * 3);
                sum.segment<3>(3 * indices[3]) += PEPx.segment<3>(3 * 3);
            }

            for (const auto& constraint : m_constraints)
            {
                const std::size_t vert_index = constraint.vert_index;
                const double&     k          = constraint.stiffness;
                const auto        p          = x.segment<3>(vert_index * 3);
                const auto        q          = constraint.motion(m_physics_time + m_delta_physics_time);
                const auto        r          = p - q;

                sum.segment<3>(3 * vert_index) += k * r;
            }

            // Collision energy
            for (const auto& collision : collisions)
            {
                const auto&  vert_index = collision.vert_index;
                const auto   pos        = x.segment<3>(vert_index * 3);
                const auto&  normal     = collision.surface_normal;
                const double a          = (collision.surface_pos - pos).transpose() * normal;

                if (a > 0.0)
                {
                    sum.segment<3>(3 * vert_index) += 2.0 * k_collision_stiffness * a * (-normal);
                }
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
        mathtoolbox::optimization::RunLBfgs(y, calcObjective, calcObjectiveGrad, 1e-06, 50, x_opt, num_iters);

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
        constexpr std::size_t num_blocks_x = 16;
        constexpr std::size_t num_blocks_y = 4;
        constexpr std::size_t num_blocks_z = 4;
        constexpr std::size_t num_verts    = (num_blocks_x + 1) * (num_blocks_y + 1) * (num_blocks_z + 1);
        constexpr std::size_t num_elems    = num_blocks_x * num_blocks_y * num_blocks_z * 5;
        constexpr double      scale        = 0.5 / num_blocks_y;

        m_mesh.elems.resize(4, num_elems);
        m_mesh.x_rest.resize(num_verts * 3);

        for (std::size_t i_z = 0; i_z < num_blocks_z; ++i_z)
        {
            for (std::size_t i_y = 0; i_y < num_blocks_y; ++i_y)
            {
                for (std::size_t i_x = 0; i_x < num_blocks_x; ++i_x)
                {
                    const std::size_t block_index = i_x + num_blocks_x * i_y + num_blocks_x * num_blocks_y * i_z;
                    const std::size_t vert_base_index =
                        i_x + (num_blocks_x + 1) * i_y + (num_blocks_x + 1) * (num_blocks_y + 1) * i_z;
                    const std::size_t     elem_base_index = 5 * block_index;
                    const Eigen::Vector3d vert_base_pos{i_x, i_y, i_z};

                    const std::size_t indices[] = {vert_base_index,
                                                   vert_base_index + 1,
                                                   vert_base_index + (num_blocks_x + 1),
                                                   vert_base_index + (num_blocks_x + 1) + 1,
                                                   vert_base_index + (num_blocks_x + 1) * (num_blocks_y + 1),
                                                   vert_base_index + (num_blocks_x + 1) * (num_blocks_y + 1) + 1,
                                                   vert_base_index + (num_blocks_x + 1) * (num_blocks_y + 2),
                                                   vert_base_index + (num_blocks_x + 1) * (num_blocks_y + 2) + 1};

                    m_mesh.elems.col(elem_base_index + 0) << indices[0], indices[4], indices[1], indices[2];
                    m_mesh.elems.col(elem_base_index + 1) << indices[2], indices[7], indices[1], indices[3];
                    m_mesh.elems.col(elem_base_index + 2) << indices[6], indices[2], indices[7], indices[4];
                    m_mesh.elems.col(elem_base_index + 3) << indices[5], indices[1], indices[4], indices[7];
                    m_mesh.elems.col(elem_base_index + 4) << indices[1], indices[2], indices[4], indices[7];

                    m_mesh.x_rest.segment<3>(indices[0] * 3) << 0.0, 0.0, 0.0;
                    m_mesh.x_rest.segment<3>(indices[1] * 3) << 1.0, 0.0, 0.0;
                    m_mesh.x_rest.segment<3>(indices[2] * 3) << 0.0, 1.0, 0.0;
                    m_mesh.x_rest.segment<3>(indices[3] * 3) << 1.0, 1.0, 0.0;
                    m_mesh.x_rest.segment<3>(indices[4] * 3) << 0.0, 0.0, 1.0;
                    m_mesh.x_rest.segment<3>(indices[5] * 3) << 1.0, 0.0, 1.0;
                    m_mesh.x_rest.segment<3>(indices[6] * 3) << 0.0, 1.0, 1.0;
                    m_mesh.x_rest.segment<3>(indices[7] * 3) << 1.0, 1.0, 1.0;

                    for (std::size_t i = 0; i < 8; ++i)
                    {
                        m_mesh.x_rest.segment<3>(indices[i] * 3) += vert_base_pos;
                    }
                }
            }
        }

        // Set transform
        for (std::size_t vert = 0; vert < num_verts; ++vert)
        {
            m_mesh.x_rest[3 * vert + 1] -= 0.5 * num_blocks_y;
            m_mesh.x_rest[3 * vert + 2] -= 0.5 * num_blocks_z;
        }
        m_mesh.x_rest *= scale;

        // Initialize other values
        m_mesh.x = m_mesh.x_rest;
        m_mesh.v = Eigen::VectorXd::Zero(3 * num_verts);
        m_mesh.f = Eigen::VectorXd::Zero(3 * num_verts);

        // Set constraints
        for (std::size_t i_z = 0; i_z < num_blocks_z + 1; ++i_z)
        {
            for (std::size_t i_y = 0; i_y < num_blocks_y + 1; ++i_y)
            {
                const std::size_t i_x = 0;

                const std::size_t vert_index =
                    i_x + (num_blocks_x + 1) * i_y + (num_blocks_x + 1) * (num_blocks_y + 1) * i_z;

                const auto motion = [&, vert_index](double) -> Eigen::Vector3d {
                    return m_mesh.x_rest.segment<3>(vert_index * 3);
                };

                m_constraints.push_back(Constraint{vert_index, motion, k_spring_stiffness});
            }
        }
        for (std::size_t i_z = 0; i_z < num_blocks_z + 1; ++i_z)
        {
            for (std::size_t i_y = 0; i_y < num_blocks_y + 1; ++i_y)
            {
                const std::size_t i_x = num_blocks_x;

                const std::size_t vert_index =
                    i_x + (num_blocks_x + 1) * i_y + (num_blocks_x + 1) * (num_blocks_y + 1) * i_z;

                const auto motion = [&, vert_index](double t) -> Eigen::Vector3d {
                    constexpr double pi = 3.1415926535897932;

                    const auto ease = [](double x) {
                        return x < 0.5 ? 4.0 * x * x * x : 1.0 - 0.5 * std::pow(-2.0 * x + 2.0, 3.0);
                    };

                    const auto   x_init    = m_mesh.x_rest.segment<3>(vert_index * 3);
                    const auto   axis      = Eigen::Vector3d::UnitX();
                    const double t_0       = 0.5;
                    const double t_1       = t_0 + 2.0;
                    const double a         = (t < t_0) ? 0.0 : ((t < t_1) ? t - t_0 : t_1 - t_0);
                    const double b         = ease(a / 2.0);
                    const double total_rot = 1.5 * pi;

                    return Eigen::AngleAxisd(total_rot * b, axis) * x_init;
                };

                m_constraints.push_back(Constraint{vert_index, motion, k_spring_stiffness});
            }
        }

        // Perform precomputation
        m_mesh.volume_array.resize(num_elems);
        m_mesh.rest_shape_mat_inv_array.resize(num_elems);
        m_mesh.vec_PFPx_array.resize(num_elems);
        for (std::size_t elem_index = 0; elem_index < num_elems; ++elem_index)
        {
            const auto& indices = m_mesh.elems.col(elem_index);

            const auto& x_0 = m_mesh.x_rest.segment<3>(3 * indices[0]);
            const auto& x_1 = m_mesh.x_rest.segment<3>(3 * indices[1]);
            const auto& x_2 = m_mesh.x_rest.segment<3>(3 * indices[2]);
            const auto& x_3 = m_mesh.x_rest.segment<3>(3 * indices[3]);

            m_mesh.volume_array[elem_index]             = elasty::fem::calcTetrahedronVolume(x_0, x_1, x_2, x_3);
            m_mesh.rest_shape_mat_inv_array[elem_index] = elasty::fem::calc3dShapeMatrix(x_0, x_1, x_2, x_3).inverse();
            m_mesh.vec_PFPx_array[elem_index] =
                elasty::fem::calcVecTetrahedronPartDeformGradPartPos(m_mesh.rest_shape_mat_inv_array[elem_index]);

            assert(!m_mesh.rest_shape_mat_inv_array[elem_index].hasNaN());
        }
        m_mesh.lumped_mass = elasty::fem::calcTetraMeshLumpedMass(m_mesh.x_rest, m_mesh.elems, m_mesh.mass);
    }

    /// \brief Getter of the delta physics time.
    ///
    /// \details The value equals to the delta frame time devided by the number of substeps.
    double getDeltaPhysicsTime() const { return m_delta_physics_time; }

    void setDeltaPhysicsTime(const double delta_physics_time) { m_delta_physics_time = delta_physics_time; }

    const TetraMesh* getMesh() const { return &m_mesh; }

private:
    double m_delta_physics_time = 1.0 / 60.0;
    double m_physics_time       = 0.0;

    std::vector<Constraint> m_constraints;

    TetraMesh m_mesh;
};

int main(int argc, char** argv)
{
    VariationalImplicit3dEngine engine;

    engine.initializeScene();
    engine.setDeltaPhysicsTime(k_delta_time / static_cast<double>(k_num_substeps));

    const auto        mesh      = engine.getMesh();
    const std::size_t num_verts = mesh->x_rest.size() / 3;
    const std::size_t num_elems = mesh->elems.cols();

    auto alembic_manager = elasty::createTetraMeshAlembicManager(
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
