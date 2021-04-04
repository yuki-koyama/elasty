#include <Alembic/AbcCoreOgawa/All.h>
#include <Alembic/AbcGeom/All.h>
#include <elasty/alembic-manager.hpp>
#include <elasty/cloth-sim-object.hpp>
#include <elasty/utils.hpp>

/// \brief A class for implementing common operations for alembic manager classes
class AlembicManagerBase : public elasty::AbstractAlembicManager
{
public:
    AlembicManagerBase(const std::string& output_file_path, const double delta_time, const std::string& object_name)
        : m_archive(Alembic::AbcCoreOgawa::WriteArchive(), output_file_path.c_str())
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const TimeSampling time_sampling(delta_time, 0);
        const uint32_t     time_sampling_index = m_archive.addTimeSampling(time_sampling);

        m_mesh_obj = OPolyMesh(OObject(m_archive, kTop), object_name.c_str());
        m_mesh_obj.getSchema().setTimeSampling(time_sampling_index);
    }

protected:
    Alembic::Abc::OArchive      m_archive;
    Alembic::AbcGeom::OPolyMesh m_mesh_obj;
};

class TriangleMesh2dAlembicManager : public AlembicManagerBase
{
public:
    TriangleMesh2dAlembicManager(const std::string&  output_file_path,
                                 const double        delta_time,
                                 const std::size_t   num_verts,
                                 const std::size_t   num_triangles,
                                 const double*       positions,
                                 const std::int32_t* indices)
        : AlembicManagerBase(output_file_path, delta_time, "Mesh"),
          m_num_verts(num_verts),
          m_num_triangles(num_triangles),
          m_positions(positions),
          m_indices(indices)
    {
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        // Extend the data into 3D by inserting zeros
        Eigen::MatrixXf verts = Eigen::Map<const Eigen::MatrixXd>(m_positions, 2, m_num_verts).cast<float>();
        verts.conservativeResize(3, Eigen::NoChange);
        verts.block(2, 0, 1, m_num_verts) = Eigen::MatrixXf::Zero(1, m_num_verts);

        // If this is the first call, set a sample with full properties including vertex positions, indices, counts;
        // otherwise, set a sample with only vertex positions.
        if (m_is_first)
        {
            const std::size_t               num_counts = m_num_triangles;
            const std::vector<std::int32_t> counts(num_counts, 3);

            // Ignore UV
            const OV2fGeomParam::Sample geom_param_sample = OV2fGeomParam::Sample();

            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), m_num_verts),
                                                 Int32ArraySample(m_indices, m_num_triangles * 3),
                                                 Int32ArraySample(counts.data(), num_counts),
                                                 geom_param_sample);

            m_mesh_obj.getSchema().set(sample);

            m_is_first = false;
        }
        else
        {
            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), m_num_verts));
            m_mesh_obj.getSchema().set(sample);
        }
    }

private:
    bool m_is_first = true;

    const std::size_t m_num_verts;
    const std::size_t m_num_triangles;

    const double*       m_positions;
    const std::int32_t* m_indices;
};

class ClothAlembicManager : public AlembicManagerBase
{
public:
    ClothAlembicManager(const std::string&                            output_file_path,
                        const double                                  delta_time,
                        const std::shared_ptr<elasty::ClothSimObject> cloth_sim_object)
        : AlembicManagerBase(output_file_path, delta_time, "Cloth"), m_cloth_sim_object(cloth_sim_object)
    {
    }

    void submitCurrentStatus() override
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const std::size_t        num_verts = m_cloth_sim_object->m_particles.size();
        const std::vector<float> verts     = packParticlePositions(m_cloth_sim_object->m_particles);

        // If this is the first call, set a sample with full properties including vertex positions, indices, counts, and
        // UVs (if exists); otherwise, set a sample with only vertex positions.
        if (m_is_first)
        {
            const std::size_t               num_indices = m_cloth_sim_object->m_triangle_list.size();
            const std::vector<std::int32_t> indices     = [&]() {
                std::vector<std::int32_t> indices;
                indices.reserve(num_indices);
                for (unsigned int i = 0; i < m_cloth_sim_object->m_triangle_list.rows(); ++i)
                {
                    indices.push_back(m_cloth_sim_object->m_triangle_list(i, 0));
                    indices.push_back(m_cloth_sim_object->m_triangle_list(i, 1));
                    indices.push_back(m_cloth_sim_object->m_triangle_list(i, 2));
                }
                return indices;
            }();
            const std::size_t               num_counts = m_cloth_sim_object->m_triangle_list.rows();
            const std::vector<std::int32_t> counts(num_counts, 3);

            // If the model has UV info, include it into the geometry sample
            const OV2fGeomParam::Sample geom_param_sample =
                m_cloth_sim_object->hasUv()
                    ? OV2fGeomParam::Sample(
                          V2fArraySample((const V2f*) m_cloth_sim_object->m_uv_list.data(), num_indices),
                          kFacevaryingScope)
                    : OV2fGeomParam::Sample();

            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), num_verts),
                                                 Int32ArraySample(indices.data(), num_indices),
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
    bool m_is_first = true;

    const std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;
};

std::shared_ptr<elasty::AbstractAlembicManager>
elasty::createClothAlembicManager(const std::string&                    file_path,
                                  const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                  const double                          delta_time)
{
    return std::make_shared<ClothAlembicManager>(file_path, delta_time, cloth_sim_object);
}

std::shared_ptr<elasty::AbstractAlembicManager>
elasty::createTriangleMesh2dAlembicManager(const std::string&  file_path,
                                           const double        delta_time,
                                           const std::size_t   num_verts,
                                           const std::size_t   num_triangles,
                                           const double*       positions,
                                           const std::int32_t* indices)
{
    return std::make_shared<TriangleMesh2dAlembicManager>(
        file_path, delta_time, num_verts, num_triangles, positions, indices);
}
