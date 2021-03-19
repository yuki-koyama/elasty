#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <timer.hpp>
#include <vector>

namespace
{
    constexpr size_t num_dims = 2;
}

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

        const size_t num_verts = m_mesh->x_rest.size() / num_dims;

        const Eigen::VectorXf verts = [&]() {
            assert(num_dims == 2);

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

// Reference:
// Dynamic Deformables: Implementation and Production Practicalities (Appendix D)
template <class Scalar>
Eigen::Matrix<Scalar, 2, 2> calc2dDeformGrad(const Scalar* x_0_data,
                                             const Scalar* x_1_data,
                                             const Scalar* x_2_data,
                                             const Scalar* x_rest_0_data,
                                             const Scalar* x_rest_1_data,
                                             const Scalar* x_rest_2_data)
{
    using Vec = Eigen::Matrix<Scalar, 2, 1>;
    using Mat = Eigen::Matrix<Scalar, 2, 2>;

    const auto x_0      = Eigen::Map<const Vec>(x_0_data);
    const auto x_1      = Eigen::Map<const Vec>(x_1_data);
    const auto x_2      = Eigen::Map<const Vec>(x_2_data);
    const auto x_rest_0 = Eigen::Map<const Vec>(x_rest_0_data);
    const auto x_rest_1 = Eigen::Map<const Vec>(x_rest_1_data);
    const auto x_rest_2 = Eigen::Map<const Vec>(x_rest_2_data);

    const auto o_1 = x_1 - x_0;
    const auto o_2 = x_2 - x_0;

    const auto o_rest_1 = x_rest_1 - x_rest_0;
    const auto o_rest_2 = x_rest_2 - x_rest_0;

    // Note: m denotes material (rest)
    Mat D_m;
    D_m.col(0) = o_rest_1;
    D_m.col(1) = o_rest_2;

    // Note: s denotes spatial (deformed)
    Mat D_s;
    D_s.col(0) = o_1;
    D_s.col(1) = o_2;

    const auto F = D_s * D_m.inverse();

    return F;
}

class Explicit2dEngine
{
public:
    Explicit2dEngine() {}

    void proceedFrame()
    {
        const size_t num_verts = m_mesh.x_rest.size() / num_dims;

        // TODO
        const auto M =
            (m_mesh.mass / static_cast<double>(num_verts)) * Eigen::VectorXd::Ones(num_dims * num_verts).asDiagonal();

        // TODO
        const auto K = 10.0 * Eigen::MatrixXd::Identity(num_dims * num_verts, num_dims * num_verts);

        // TODO
        m_mesh.f = -K * (m_mesh.x - m_mesh.x_rest);

        // Note: Explicit Euler integration
        m_mesh.v = m_mesh.v + m_delta_physics_time * M.inverse() * m_mesh.f;
        m_mesh.x = m_mesh.x + m_delta_physics_time * m_mesh.v;
    }

    void initializeScene()
    {
        constexpr size_t num_verts = 3;
        constexpr size_t num_elems = 1;

        m_mesh.elems.resize(num_elems, 3);
        m_mesh.elems(0, 0) = 0;
        m_mesh.elems(0, 1) = 1;
        m_mesh.elems(0, 2) = 2;

        m_mesh.x_rest.resize(num_verts * num_dims);
        m_mesh.x_rest.segment(num_dims * 0, num_dims) = Eigen::Vector2d{0.0, 0.0};
        m_mesh.x_rest.segment(num_dims * 1, num_dims) = Eigen::Vector2d{1.0, 0.0};
        m_mesh.x_rest.segment(num_dims * 2, num_dims) = Eigen::Vector2d{0.0, 1.0};

        m_mesh.x = m_mesh.x_rest;
        m_mesh.v = Eigen::VectorXd::Random(num_dims * num_verts);
        m_mesh.f = Eigen::VectorXd::Zero(num_dims * num_verts);
    }

    /// \brief Getter of the delta physics time.
    ///
    /// \details The value equals to the delta frame time devided by the number of substeps.
    double getDeltaPhysicsTime() const { return m_delta_physics_time; }

    const TriangleMesh* getMesh() const { return &m_mesh; }

private:
    const double m_delta_physics_time = 1.0 / 60.0;

    double m_current_physics_time = 0.0;

    TriangleMesh m_mesh;
};

int main(int argc, char** argv)
{
    Explicit2dEngine engine;
    engine.initializeScene();

    auto alembic_manager = AlembicManager("./out.abc", engine.getMesh(), engine.getDeltaPhysicsTime());

    for (unsigned int frame = 0; frame < 90; ++frame)
    {
        timer::Timer t(std::to_string(frame));

        alembic_manager.submitCurrentStatus();

        engine.proceedFrame();
    }

    return 0;
}
