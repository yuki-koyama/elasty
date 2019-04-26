#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <elasty/particle.hpp>
#include <elasty/utils.hpp>
#include <timer.hpp>

class SimpleEngine final : public elasty::Engine
{
public:

    void initializeScene() override
    {
        m_num_iterations = 30;

        // Instantiate a cloth object
        const double cloth_distance_stiffness = 0.95;
        const double cloth_bending_stiffness = 0.03;
        const std::string cloth_obj_path = "./models/cloths/0.10.obj";
#if 0 // Drape
        const Eigen::Affine3d cloth_import_transform = Eigen::Translation3d(0.0, 1.0, 0.0) * Eigen::AngleAxisd(0.5 * elasty::pi(), Eigen::Vector3d::UnitX());
#else // Fall
        const Eigen::Affine3d cloth_import_transform = Eigen::Affine3d(Eigen::Translation3d(0.0, 2.0, 1.0));
#endif
        m_cloth_sim_object = std::make_shared<elasty::ClothSimObject>(cloth_obj_path,
                                                                      cloth_distance_stiffness,
                                                                      cloth_bending_stiffness,
                                                                      cloth_import_transform,
                                                                      elasty::ClothSimObject::Strategy::IsometricBending);

        // Register the cloth object
        std::copy(m_cloth_sim_object->m_particles.begin(),
                  m_cloth_sim_object->m_particles.end(),
                  std::back_inserter(m_particles));
        std::copy(m_cloth_sim_object->m_constraints.begin(),
                  m_cloth_sim_object->m_constraints.end(),
                  std::back_inserter(m_constraints));

        // Pin two of the corners of the cloth
        constexpr double range_radius = 0.1;
        for (const auto& particle : m_particles)
        {
            if ((particle->x - Eigen::Vector3d(+ 1.0, 2.0, 0.0)).norm() < range_radius)
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(particle, 1.0, particle->x));
            }
            if ((particle->x - Eigen::Vector3d(- 1.0, 2.0, 0.0)).norm() < range_radius)
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(particle, 1.0, particle->x));
            }
        }
    }

    void setExternalForces() override
    {
        const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, - 9.8, 0.0);

        for (auto particle : m_particles)
        {
            particle->f = particle->m * gravity;
        }
    }

    void generateCollisionConstraints() override {}
    void updateVelocities() override {}

    std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;
};

int main(int argc, char** argv)
{
    SimpleEngine engine;
    engine.initializeScene();

    auto alembic_manager = elasty::createAlembicManager("./cloth.abc", engine.m_cloth_sim_object, engine.m_dt);

    for (unsigned int frame = 0; frame < 180; ++ frame)
    {
        timer::Timer t(std::to_string(frame));
        elasty::submitCurrentStatus(alembic_manager);
        engine.stepTime();
    }

    return 0;
}
