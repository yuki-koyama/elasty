#include <elasty/utils.hpp>
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
        const Alembic::Abc::TimeSampling time_sampling(dt, 0);
        const uint32_t time_sampling_index = m_archive.addTimeSampling(time_sampling);

        Alembic::AbcGeom::OPolyMesh mesh_obj(Alembic::Abc::OObject(m_archive, Alembic::Abc::kTop), "mesh");
        m_mesh_ptr = &(mesh_obj.getSchema());

        m_mesh_ptr->setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        using namespace Alembic::Abc;
        using namespace Alembic::AbcGeom;

        constexpr size_t num_verts = 8;
        constexpr float verts[] =
        {
            - 1.0f, - 1.0f, - 1.0f,
            + 1.0f, - 1.0f, - 1.0f,
            - 1.0f, + 1.0f, - 1.0f,
            + 1.0f, + 1.0f, - 1.0f,
            - 1.0f, - 1.0f, + 1.0f,
            + 1.0f, - 1.0f, + 1.0f,
            - 1.0f, + 1.0f, + 1.0f,
            + 1.0f, + 1.0f, + 1.0f,
        };
        constexpr size_t num_indices = 24;
        constexpr int32_t indices[] =
        {
            0, 4, 6, 2,
            5, 1, 3, 7,
            0, 1, 5, 4,
            6, 7, 3, 2,
            1, 0, 2, 3,
            4, 5, 7, 6,
        };
        constexpr size_t num_counts = 6;
        constexpr int32_t counts[] = { 4, 4, 4, 4, 4, 4 };

        if (m_is_first)
        {
            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts, num_verts),
                                                 Int32ArraySample(indices, num_indices),
                                                 Int32ArraySample(counts, num_counts));
            m_mesh_ptr->set(sample);

            m_is_first = false;
        }
        else
        {
            const OPolyMeshSchema::Sample sample(V3fArraySample((const V3f*) verts, num_verts));
            m_mesh_ptr->set(sample);
        }
    }

private:

    bool m_is_first = true;
    const std::shared_ptr<ClothSimObject> m_cloth_sim_object;
    Alembic::Abc::OArchive m_archive;
    Alembic::AbcGeom::OPolyMeshSchema* m_mesh_ptr;
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
