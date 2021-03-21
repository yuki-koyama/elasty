#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <elasty/fem.hpp>
#include <timer.hpp>
#include <vector>

namespace
{
    constexpr size_t k_num_dims = 2;

    constexpr double k_youngs_modulus = 200.0;
    constexpr double k_poisson_ratio  = 0.45;

    constexpr double k_first_lame  = elasty::fem::calcFirstLame(k_youngs_modulus, k_poisson_ratio);
    constexpr double k_second_lame = elasty::fem::calcSecondLame(k_youngs_modulus, k_poisson_ratio);

    constexpr unsigned k_num_substeps = 10;
    constexpr double   k_delta_time   = 1.0 / 60.0;
} // namespace

struct TriangleMesh
{
    using TriangleList = Eigen::Matrix<int32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

    TriangleList elems;

    Eigen::VectorXd x_rest;
    Eigen::VectorXd x;
    Eigen::VectorXd v;
    Eigen::VectorXd f;

    double mass = 1.0;
};

class AlembicManager
{
public:
    AlembicManager(const std::string& output_file_path, const TriangleMesh* mesh, const double delta_time = 1.0 / 60.0)
        : m_mesh(mesh), m_archive(Alembic::AbcCoreOgawa::WriteArchive(), output_file_path.c_str())
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const TimeSampling time_sampling(delta_time, 0);
        const uint32_t     time_sampling_index = m_archive.addTimeSampling(time_sampling);

        m_mesh_obj = OPolyMesh(OObject(m_archive, kTop), "mesh");
        m_mesh_obj.getSchema().setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const size_t num_verts = m_mesh->x_rest.size() / k_num_dims;

        const Eigen::VectorXf verts = [&]() {
            assert(k_num_dims == 2);

            Eigen::VectorXf temp = Eigen::VectorXf::Zero(3 * num_verts);

            for (size_t i = 0; i < num_verts; ++i)
            {
                temp.segment(3 * i, 2) = m_mesh->x.segment(2 * i, 2).cast<float>();
            }

            return temp;
        }();

        // If this is the first call, set a sample with full properties including vertex positions, indices, counts, and
        // UVs (if exists); otherwise, set a sample with only vertex positions.
        if (m_is_first)
        {
            const size_t num_indices = m_mesh->elems.size();

            const size_t               num_counts = num_indices / 3;
            const std::vector<int32_t> counts(num_counts, 3);

            // Ignore UV
            const OV2fGeomParam::Sample geom_param_sample = OV2fGeomParam::Sample();

            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), num_verts),
                                                 Int32ArraySample(m_mesh->elems.data(), num_indices),
                                                 Int32ArraySample(counts.data(), num_counts),
                                                 geom_param_sample);

            m_mesh_obj.getSchema().set(sample);

            m_is_first = false;
        }
        else
        {
            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), num_verts));
            m_mesh_obj.getSchema().set(sample);
        }
    }

private:
    bool                        m_is_first = true;
    const TriangleMesh*         m_mesh;
    Alembic::Abc::OArchive      m_archive;
    Alembic::AbcGeom::OPolyMesh m_mesh_obj;
};

template <typename Derived> typename Derived::Scalar calcEnergyDensity(const Eigen::MatrixBase<Derived>& deform_grad)
{
    return elasty::fem::calcStVenantKirchhoffEnergyDensity(deform_grad, k_first_lame, k_second_lame);
}

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, 2, 2> calcPiolaStress(const Eigen::MatrixBase<Derived>& deform_grad)
{
    return elasty::fem::calcStVenantKirchhoffPiolaStress(deform_grad, k_first_lame, k_second_lame);
}

// Reference:
// Dynamic Deformables: Implementation and Production Practicalities (Appendix D)
template <class Scalar>
Eigen::Matrix<Scalar, 4, 6> calcFlattenedPartDeformGradPartPos(const Eigen::Matrix<Scalar, 2, 2>& D_m_inv)
{
    // Analytics partial derivatives $\frac{\partial \mathbf{D}_\text{s}{\partial x_{i}}$
    Eigen::Matrix<Scalar, 2, 2> PDsPx[6];

    // x[0] (= x_0)
    PDsPx[0] << -1.0, -1.0, 0.0, 0.0;

    // x[1] (= y_0)
    PDsPx[1] << 0.0, 0.0, -1.0, -1.0;

    // x[2] (= x_1)
    PDsPx[2] << 1.0, 0.0, 0.0, 0.0;

    // x[3] (= y_1)
    PDsPx[3] << 0.0, 0.0, 1.0, 0.0;

    // x[4] (= x_2)
    PDsPx[4] << 0.0, 1.0, 0.0, 0.0;

    // x[5] (= y_2)
    PDsPx[5] << 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<Scalar, 4, 6> vec_PFPx;
    for (size_t i = 0; i < 6; ++i)
    {
        const Eigen::Matrix<Scalar, 2, 2> temp = PDsPx[i] * D_m_inv;

        vec_PFPx.col(i) = Eigen::Map<const Eigen::Matrix<Scalar, 4, 1>>(temp.data(), temp.size());
    }

    return vec_PFPx;
}

class Explicit2dEngine
{
public:
    Explicit2dEngine() {}

    void proceedFrame()
    {
        const size_t              num_verts         = m_mesh.x_rest.size() / k_num_dims;
        const std::vector<size_t> constrained_verts = {0, 1, 2, 3, 4};

        // TODO
        const auto M =
            (m_mesh.mass / static_cast<double>(num_verts)) * Eigen::VectorXd::Ones(k_num_dims * num_verts).asDiagonal();

        // Reset forces
        m_mesh.f = Eigen::VectorXd::Zero(2 * num_verts);

        // Calculate forces
        for (size_t i = 0; i < m_mesh.elems.rows(); ++i)
        {
            const auto& indices = m_mesh.elems.row(i);

            // TODO: Precompute this value
            const Eigen::Matrix2d D_m_inv = elasty::fem::calc2dShapeMatrix(m_mesh.x_rest.segment(2 * indices[0], 2),
                                                                           m_mesh.x_rest.segment(2 * indices[1], 2),
                                                                           m_mesh.x_rest.segment(2 * indices[2], 2))
                                                .inverse();

            // TODO: Precompute a part of this process
            const auto F = elasty::fem::calc2dTriangleDeformGrad(m_mesh.x.segment(2 * indices[0], 2),
                                                                 m_mesh.x.segment(2 * indices[1], 2),
                                                                 m_mesh.x.segment(2 * indices[2], 2),
                                                                 D_m_inv);

            // TODO: Precompute this value
            const auto area = elasty::fem::calc2dTriangleArea(m_mesh.x_rest.segment(2 * indices[0], 2),
                                                              m_mesh.x_rest.segment(2 * indices[1], 2),
                                                              m_mesh.x_rest.segment(2 * indices[2], 2));

            // TODO: Precompute this value
            // TODO: Assign a better name
            const auto flattened_PFPx = calcFlattenedPartDeformGradPartPos(D_m_inv);

            // TODO: Assign a better name
            const auto pk1 = calcPiolaStress(F);

            // TODO: Assign a better name
            const auto PPsiPx = flattened_PFPx.transpose() * Eigen::Map<const Eigen::Vector4d>(pk1.data(), pk1.size());

            const auto force = -area * PPsiPx;

            assert(indices.size() == 3);
            assert(F.rows() == 2 && F.cols() == 2);
            assert(flattened_PFPx.rows() == 4 && flattened_PFPx.cols() == 6);
            assert(PPsiPx.rows() == 6 && PPsiPx.cols() == 1);

            // Accumulate forces
            m_mesh.f.segment(2 * indices[0], 2) += force.segment(2 * 0, 2);
            m_mesh.f.segment(2 * indices[1], 2) += force.segment(2 * 1, 2);
            m_mesh.f.segment(2 * indices[2], 2) += force.segment(2 * 2, 2);

            // Validate the force calculation by comparing to numerical differentials
            const auto      x_orig = m_mesh.x;
            Eigen::VectorXd diff{6};
            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 2; ++j)
                {
                    constexpr double eps = 1e-06;

                    m_mesh.x[2 * indices[i] + j] += eps;
                    const auto F_p = elasty::fem::calc2dTriangleDeformGrad(m_mesh.x.segment(2 * indices[0], 2),
                                                                           m_mesh.x.segment(2 * indices[1], 2),
                                                                           m_mesh.x.segment(2 * indices[2], 2),
                                                                           D_m_inv);
                    const auto e_p = calcEnergyDensity(F_p);

                    m_mesh.x = x_orig;

                    m_mesh.x[2 * indices[i] + j] -= eps;
                    const auto F_m = elasty::fem::calc2dTriangleDeformGrad(m_mesh.x.segment(2 * indices[0], 2),
                                                                           m_mesh.x.segment(2 * indices[1], 2),
                                                                           m_mesh.x.segment(2 * indices[2], 2),
                                                                           D_m_inv);
                    const auto e_m = calcEnergyDensity(F_m);

                    m_mesh.x = x_orig;

                    diff[2 * i + j] = (e_p - e_m) / (2.0 * eps);
                }
            }
            assert(std::abs(PPsiPx.norm() - diff.norm()) < 1e-04);
        }

        // Apply gravity force
        for (size_t i = 0; i < num_verts; ++i)
        {
            m_mesh.f[i * 2 + 1] += M.diagonal()(i * 2 + 1) * (-9.80665);
        }

        // Note: Explicit Euler integration (for velocities)
        m_mesh.v = m_mesh.v + m_delta_physics_time * M.inverse() * m_mesh.f;

        // Naive constraints
        for (size_t i : constrained_verts)
        {
            m_mesh.v.segment(i * 2, 2) = Eigen::Vector2d::Zero();
        }

        // Note: Explicit Euler integration (for positions)
        m_mesh.x = m_mesh.x + m_delta_physics_time * m_mesh.v;

        // Apply naive damping
        constexpr double damping_factor = 0.5;

        m_mesh.v *= std::exp(-damping_factor * m_delta_physics_time);
    }

    void initializeScene()
    {
        if constexpr (false)
        {
            // A simple square

            constexpr size_t num_verts = 4;
            constexpr size_t num_elems = 2;

            m_mesh.elems.resize(num_elems, 3);
            m_mesh.elems.row(0) << 0, 1, 2;
            m_mesh.elems.row(1) << 1, 3, 2;

            m_mesh.x_rest.resize(num_verts * k_num_dims);
            m_mesh.x_rest.segment(k_num_dims * 0, k_num_dims) = Eigen::Vector2d{0.0, 0.0};
            m_mesh.x_rest.segment(k_num_dims * 1, k_num_dims) = Eigen::Vector2d{1.0, 0.0};
            m_mesh.x_rest.segment(k_num_dims * 2, k_num_dims) = Eigen::Vector2d{0.0, 1.0};
            m_mesh.x_rest.segment(k_num_dims * 3, k_num_dims) = Eigen::Vector2d{1.0, 1.0};

            m_mesh.x = m_mesh.x_rest;
            m_mesh.v = Eigen::VectorXd::Random(k_num_dims * num_verts);
            m_mesh.f = Eigen::VectorXd::Zero(k_num_dims * num_verts);
        }
        else
        {
            // A simple cantilever

            constexpr size_t num_cols  = 20;
            constexpr size_t num_rows  = 4;
            constexpr size_t num_verts = (num_cols + 1) * (num_rows + 1);
            constexpr size_t num_elems = (num_cols * num_rows) * 2;
            constexpr double size      = 1.0;

            m_mesh.elems.resize(num_elems, 3);
            m_mesh.x_rest.resize(num_verts * k_num_dims);

            for (size_t col = 0; col < num_cols; ++col)
            {
                for (size_t row = 0; row < num_rows; ++row)
                {
                    const auto base = col * (num_rows + 1) + row;

                    m_mesh.elems.row(2 * num_rows * col + 2 * row + 0) << 0 + base, 1 + base, (num_rows + 1) + 1 + base;
                    m_mesh.elems.row(2 * num_rows * col + 2 * row + 1) << 0 + base, (num_rows + 1) + 1 + base,
                        (num_rows + 1) + base;

                    m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * col + row), k_num_dims) =
                        Eigen::Vector2d{col * 1.0, -1.0 * row};
                }
                m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * col + num_rows), k_num_dims) =
                    Eigen::Vector2d{col * 1.0, -1.0 * num_rows};
            }
            size_t col = num_cols;
            for (size_t row = 0; row < num_rows; ++row)
            {
                m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * col + row), k_num_dims) =
                    Eigen::Vector2d{col * 1.0, -1.0 * row};
            }
            m_mesh.x_rest.segment(k_num_dims * ((num_rows + 1) * col + num_rows), k_num_dims) =
                Eigen::Vector2d{col * 1.0, -1.0 * num_rows};

            // Set transform
            m_mesh.x_rest *= 1.0 / static_cast<double>(num_rows);
            for (size_t vert = 0; vert < num_verts; ++vert)
            {
                m_mesh.x_rest[2 * vert + 1] += 0.5;
            }
            m_mesh.x_rest *= size;

            // Initialize other values
            m_mesh.x = m_mesh.x_rest;
            m_mesh.v = Eigen::VectorXd::Zero(k_num_dims * num_verts);
            m_mesh.f = Eigen::VectorXd::Zero(k_num_dims * num_verts);
        }
    }

    /// \brief Getter of the delta physics time.
    ///
    /// \details The value equals to the delta frame time devided by the number of substeps.
    double getDeltaPhysicsTime() const { return m_delta_physics_time; }

    void setDeltaPhysicsTime(const double delta_physics_time) { m_delta_physics_time = delta_physics_time; }

    const TriangleMesh* getMesh() const { return &m_mesh; }

private:
    double m_delta_physics_time = 1.0 / 60.0;

    TriangleMesh m_mesh;
};

int main(int argc, char** argv)
{
    Explicit2dEngine engine;

    engine.initializeScene();
    engine.setDeltaPhysicsTime(k_delta_time / static_cast<double>(k_num_substeps));

    auto alembic_manager = AlembicManager("./out.abc", engine.getMesh(), k_delta_time);

    for (unsigned int frame = 0; frame < 240; ++frame)
    {
        timer::Timer t(std::to_string(frame));

        alembic_manager.submitCurrentStatus();

        for (int i = 0; i < k_num_substeps; ++i)
        {
            engine.proceedFrame();
        }
    }

    return 0;
}
