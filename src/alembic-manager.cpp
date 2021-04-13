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

        const TimeSampling  time_sampling(delta_time, 0);
        const std::uint32_t time_sampling_index = m_archive.addTimeSampling(time_sampling);

        m_mesh_obj = OPolyMesh(OObject(m_archive, kTop), object_name.c_str());
        m_mesh_obj.getSchema().setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        if (m_is_first)
        {
            submitCurrentStatusFirstTime();

            m_is_first = false;

            return;
        }

        const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) getAlembicPositionArray(), getNumVerts()));

        m_mesh_obj.getSchema().set(sample);
    }

    /// \brief A function that is called at the first frame of the alembic export.
    ///
    /// \details At the first frame, the sample needs to include indices, UVs, etc., so special care is necessary for
    /// each individual inherited class.
    virtual void submitCurrentStatusFirstTime() = 0;

protected:
    virtual const float* getAlembicPositionArray() = 0;
    virtual std::int32_t getNumVerts() const       = 0;

    bool m_is_first = true;

    Alembic::Abc::OArchive      m_archive;
    Alembic::AbcGeom::OPolyMesh m_mesh_obj;
};

class TriangleMesh2dAlembicManager : public AlembicManagerBase
{
public:
    TriangleMesh2dAlembicManager(const std::string&  output_file_path,
                                 const double        delta_time,
                                 const std::size_t   num_verts,
                                 const std::size_t   num_elems,
                                 const double*       positions,
                                 const std::int32_t* indices)
        : AlembicManagerBase(output_file_path, delta_time, "Mesh"),
          m_num_verts(num_verts),
          m_num_elems(num_elems),
          m_positions(positions),
          m_indices(indices)
    {
        m_alembic_position_buffer = Eigen::Matrix<float, 3, Eigen::Dynamic>::Zero(3, m_num_verts);
    }

    void submitCurrentStatusFirstTime() override
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const std::size_t               num_counts = m_num_elems;
        const std::vector<std::int32_t> counts(num_counts, 3);

        // Ignore UV
        const OV2fGeomParam::Sample geom_param_sample = OV2fGeomParam::Sample();

        const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) getAlembicPositionArray(), getNumVerts()),
                                             Int32ArraySample(m_indices, m_num_elems * 3),
                                             Int32ArraySample(counts.data(), num_counts),
                                             geom_param_sample);

        m_mesh_obj.getSchema().set(sample);
    }

private:
    const float* getAlembicPositionArray() override
    {
        m_alembic_position_buffer.block(0, 0, 2, m_num_verts) =
            Eigen::Map<const Eigen::MatrixXd>(m_positions, 2, m_num_verts).cast<float>();

        return m_alembic_position_buffer.data();
    }

    std::int32_t getNumVerts() const override { return m_num_verts; }

    const std::size_t m_num_verts;
    const std::size_t m_num_elems;

    const double*       m_positions;
    const std::int32_t* m_indices;

    /// \brief A buffer object to store the vertex position array in the alembic format.
    Eigen::Matrix<float, 3, Eigen::Dynamic> m_alembic_position_buffer;
};

/// \details This manager assumes that the number of particles does not change during simulation
class ClothAlembicManager : public AlembicManagerBase
{
public:
    ClothAlembicManager(const std::string&                            output_file_path,
                        const double                                  delta_time,
                        const std::shared_ptr<elasty::ClothSimObject> cloth_sim_object)
        : AlembicManagerBase(output_file_path, delta_time, "Cloth"), m_cloth_sim_object(cloth_sim_object)
    {
    }

    void submitCurrentStatusFirstTime() override
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const std::size_t               num_indices = m_cloth_sim_object->m_triangle_list.size();
        const std::size_t               num_counts  = m_cloth_sim_object->m_triangle_list.rows();
        const std::vector<std::int32_t> counts(num_counts, 3);

        // If the model has UV info, include it into the geometry sample
        const OV2fGeomParam::Sample geom_param_sample =
            m_cloth_sim_object->hasUv()
                ? OV2fGeomParam::Sample(V2fArraySample((const V2f*) m_cloth_sim_object->m_uv_list.data(), num_indices),
                                        kFacevaryingScope)
                : OV2fGeomParam::Sample();

        const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) getAlembicPositionArray(), getNumVerts()),
                                             Int32ArraySample(m_cloth_sim_object->m_triangle_list.data(), num_indices),
                                             Int32ArraySample(counts.data(), num_counts),
                                             geom_param_sample);

        m_mesh_obj.getSchema().set(sample);
    }

private:
    const float* getAlembicPositionArray() override
    {
        m_alembic_position_buffer = packParticlePositions(m_cloth_sim_object->m_particles);

        return m_alembic_position_buffer.data();
    }

    std::int32_t getNumVerts() const override { return m_cloth_sim_object->m_particles.size(); }

    bool m_is_first = true;

    const std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;

    /// \brief A buffer object to store the vertex position array in the alembic format.
    std::vector<float> m_alembic_position_buffer;
};

std::shared_ptr<elasty::AbstractAlembicManager>
elasty::createClothAlembicManager(const std::string&                    file_path,
                                  const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                  const double                          delta_time)
{
    return std::make_shared<ClothAlembicManager>(file_path, delta_time, cloth_sim_object);
}

std::shared_ptr<elasty::AbstractAlembicManager> elasty::createTriangleMesh2dAlembicManager(const std::string& file_path,
                                                                                           const double      delta_time,
                                                                                           const std::size_t num_verts,
                                                                                           const std::size_t num_elems,
                                                                                           const double*     positions,
                                                                                           const std::int32_t* indices)
{
    return std::make_shared<TriangleMesh2dAlembicManager>(
        file_path, delta_time, num_verts, num_elems, positions, indices);
}
