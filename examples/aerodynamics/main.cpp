#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <elasty/particle.hpp>
#include <elasty/utils.hpp>
#include <timer.hpp>

class SimpleEngine final : public elasty::AbstractEngine
{
public:
    SimpleEngine(const unsigned num_constraint_iters, const elasty::AlgorithmType type)
        : elasty::AbstractEngine(1.0 / 60.0, num_constraint_iters, 1, type)
    {
    }

    void initializeScene() override
    {
        // Instantiate a cloth object
        constexpr double   cloth_in_plane_stiffness      = 1.000; ///< PBD
        constexpr double   cloth_in_plane_compliance     = 5e-02; ///< XPBD
        constexpr double   cloth_out_of_plane_stiffness  = 0.100; ///< PBD
        constexpr double   cloth_out_of_plane_compliance = 5e+04; ///< XPBD
        constexpr unsigned cloth_resolution              = 40;

        const Eigen::Affine3d cloth_import_transform =
            Eigen::Translation3d(0.0, 1.0, 0.0) * Eigen::AngleAxisd(0.5 * elasty::pi(), Eigen::Vector3d::UnitX());

        m_cloth_sim_object =
            std::make_shared<elasty::ClothSimObject>(cloth_resolution,
                                                     cloth_in_plane_stiffness,
                                                     cloth_in_plane_compliance,
                                                     cloth_out_of_plane_stiffness,
                                                     cloth_out_of_plane_compliance,
                                                     getDeltaPhysicsTime(),
                                                     cloth_import_transform,
                                                     elasty::ClothSimObject::InPlaneStrategy::EdgeDistance,
                                                     elasty::ClothSimObject::OutOfPlaneStrategy::IsometricBending);

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
            if ((particle->x - Eigen::Vector3d(+1.0, 2.0, 0.0)).norm() < range_radius)
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(
                    particle, 1.0, 0.0, getDeltaPhysicsTime(), particle->x));
            }
            if ((particle->x - Eigen::Vector3d(-1.0, 2.0, 0.0)).norm() < range_radius)
            {
                m_constraints.push_back(std::make_shared<elasty::FixedPointConstraint>(
                    particle, 1.0, 0.0, getDeltaPhysicsTime(), particle->x));
            }
        }
    }

    void setExternalForces() override
    {
        const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, -9.8, 0.0);

        for (auto particle : m_particles)
        {
            particle->f = particle->m * gravity;
        }
    }

    void generateCollisionConstraints() override {}

    void updateVelocities() override
    {
        for (auto particle : m_particles)
        {
            particle->v *= 0.95;
        }
    }

    std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;

    int m_count = 0;
};

int main(int argc, char** argv)
{
    constexpr unsigned              num_iters_conditions[] = {2, 4, 8, 16, 32, 64, 128, 256, 512};
    constexpr elasty::AlgorithmType type_conditions[]      = {
        elasty::AlgorithmType::Pbd,
        elasty::AlgorithmType::Xpbd,
    };

    for (const auto type : type_conditions)
    {
        const std::string type_name = (type == elasty::AlgorithmType::Pbd) ? "pbd" : "xpbd";

        for (const auto num_iters : num_iters_conditions)
        {
            const std::string name = type_name + "-" + std::to_string(num_iters);

            const timer::Timer t(name);

            SimpleEngine engine(num_iters, type);
            engine.initializeScene();

            auto alembic_manager =
                elasty::createAlembicManager(name + ".abc", engine.m_cloth_sim_object, engine.getDeltaPhysicsTime());

            for (unsigned int frame = 0; frame < 300; ++frame)
            {
                elasty::submitCurrentStatus(alembic_manager);

                engine.m_count = frame;
                engine.proceedFrame();
            }

            elasty::exportCurrentClothStateAsObj(name + ".obj", engine.m_cloth_sim_object);
        }
    }

    return 0;
}
