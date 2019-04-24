#include <elasty/utils.hpp>
#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/particle.hpp>
#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>

void elasty::generateFixedPointConstraints(const Eigen::Vector3d& search_position,
                                           const Eigen::Vector3d& fixed_position,
                                           const std::vector<std::shared_ptr<Particle>>& particles,
                                           std::vector<std::shared_ptr<Constraint>>& constraints)
{
    for (const auto& particle : particles)
    {
        if (particle->x.isApprox(search_position))
        {
            constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(particle, 1.0, fixed_position));
        }
    }
}

class elasty::AlembicManager
{
public:
    AlembicManager(const std::string& file_path,
                   const std::shared_ptr<ClothSimObject> cloth_sim_object,
                   const double dt = 1.0 / 60.0) :
    m_cloth_sim_object(cloth_sim_object),
    m_archive(Alembic::AbcCoreOgawa::WriteArchive(), file_path.c_str())
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const TimeSampling time_sampling(dt, 0);
        const uint32_t time_sampling_index = m_archive.addTimeSampling(time_sampling);

        m_mesh_obj = OPolyMesh(OObject(m_archive, kTop), "cloth");
        m_mesh_obj.getSchema().setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        const size_t num_verts = m_cloth_sim_object->m_particles.size();
        const std::vector<float> verts = [&]()
        {
            std::vector<float> verts;
            verts.reserve(3 * num_verts);
            for (const auto& particle : m_cloth_sim_object->m_particles)
            {
                verts.push_back(particle->x(0));
                verts.push_back(particle->x(1));
                verts.push_back(particle->x(2));
            }
            return verts;
        }();
        const size_t num_indices = m_cloth_sim_object->m_triangle_list.size();
        const std::vector<int32_t> indices = [&]()
        {
            std::vector<int32_t> indices;
            indices.reserve(num_indices);
            for (unsigned int i = 0; i < m_cloth_sim_object->m_triangle_list.rows(); ++ i)
            {
                indices.push_back(m_cloth_sim_object->m_triangle_list(i, 0));
                indices.push_back(m_cloth_sim_object->m_triangle_list(i, 1));
                indices.push_back(m_cloth_sim_object->m_triangle_list(i, 2));
            }
            return indices;
        }();
        const size_t num_counts = m_cloth_sim_object->m_triangle_list.rows();
        const std::vector<int32_t> counts(num_counts, 3);

        if (m_is_first)
        {
            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts.data(), num_verts),
                                                 Int32ArraySample(indices.data(), num_indices),
                                                 Int32ArraySample(counts.data(), num_counts));
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
    const std::shared_ptr<ClothSimObject> m_cloth_sim_object;
    Alembic::Abc::OArchive m_archive;
    Alembic::AbcGeom::OPolyMesh m_mesh_obj;
};

std::shared_ptr<elasty::AlembicManager> elasty::createAlembicManager(const std::string& file_path,
                                                                     const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                                                     const double dt)
{
    return std::make_shared<AlembicManager>(file_path, cloth_sim_object, dt);
}

void elasty::submitCurrentStatus(const std::shared_ptr<AlembicManager> alembic_manager)
{
    alembic_manager->submitCurrentStatus();
}
