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
    m_cloth_sim_object(cloth_sim_object)
    {
        m_archive = std::make_unique<Alembic::Abc::OArchive>(Alembic::AbcCoreOgawa::WriteArchive(), file_path.c_str());

        const Alembic::Abc::TimeSampling time_sampling(dt, 0);
        const uint32_t time_sampling_index = m_archive->addTimeSampling(time_sampling);

        Alembic::AbcGeom::OPolyMesh mesh_obj(Alembic::Abc::OObject(*(m_archive.get()), Alembic::Abc::kTop), "mesh");
        Alembic::AbcGeom::OPolyMeshSchema& mesh = mesh_obj.getSchema();

        mesh.setTimeSampling(time_sampling_index);
    }

    void submitCurrentStatus()
    {
        // TODO
        m_is_first = false;
    }

    bool m_is_first = true;
    std::unique_ptr<Alembic::Abc::OArchive> m_archive;
    std::shared_ptr<ClothSimObject> m_cloth_sim_object;
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
