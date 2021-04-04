#include <elasty/cloth-sim-object.hpp>
#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <elasty/particle.hpp>
#include <elasty/utils.hpp>
#include <timer.hpp>

class SimpleEngine final : public elasty::AbstractEngine
{
public:
    SimpleEngine(const double draf_coeff, const double lift_coeff, const double wind_velocity)
        : elasty::AbstractEngine(1.0 / 60.0, 16, 4, elasty::AlgorithmType::Xpbd),
          m_wind_velocity(0.0, 0.0, wind_velocity),
          m_drag_coeff(draf_coeff),
          m_lift_coeff(lift_coeff)
    {
    }

    void initializeScene() override
    {
        // Instantiate a cloth object
        constexpr double   cloth_in_plane_stiffness      = 1.000; ///< PBD
        constexpr double   cloth_in_plane_compliance     = 5e-02; ///< XPBD
        constexpr double   cloth_out_of_plane_stiffness  = 0.100; ///< PBD
        constexpr double   cloth_out_of_plane_compliance = 5e+04; ///< XPBD
        constexpr unsigned cloth_resolution              = 60;

        const Eigen::Affine3d cloth_import_transform = Eigen::Affine3d(Eigen::Translation3d(0.0, 2.0, 1.0));

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

        m_cloth_sim_object->applyAerodynamicForces(m_wind_velocity, m_drag_coeff, m_lift_coeff);
    }

    void generateCollisionConstraints() override
    {
        // Collision with a sphere
        const Eigen::Vector3d center(0.0, 1.0, 0.0);
        constexpr double      tolerance  = 0.05;
        constexpr double      radius     = 0.50 + 0.02;
        constexpr double      stiffness  = 1.00;
        constexpr double      compliance = 0.00;
        for (auto particle : m_particles)
        {
            const Eigen::Vector3d direction = particle->x - center;
            if (direction.norm() < radius + tolerance)
            {
                const Eigen::Vector3d normal   = direction.normalized();
                const double          distance = center.transpose() * normal + radius;
                m_instant_constraints.push_back(std::make_shared<elasty::EnvironmentalCollisionConstraint>(
                    particle, stiffness, compliance, getDeltaPhysicsTime(), normal, distance));
            }
        }
    }

    void updateVelocities() override
    {
        const double decay_rate = std::exp(std::log(0.95) * getDeltaPhysicsTime());

        for (auto particle : m_particles)
        {
            particle->v *= decay_rate;
        }
    }

    std::shared_ptr<elasty::ClothSimObject> m_cloth_sim_object;

private:
    const Eigen::Vector3d m_wind_velocity;
    const double          m_drag_coeff;
    const double          m_lift_coeff;
};

int main(int argc, char** argv)
{
    constexpr std::tuple<double, double, double, const char*> conditions[] = {
        {0.000, 0.000, 0.0, "without-aerodynamics"},
        {0.060, 0.030, 0.0, "with-aerodynamics"},
        {0.080, 0.030, 8.0, "wind"},
        {0.080, 0.000, 8.0, "wind-high-drag"},
        {0.080, 0.080, 8.0, "wind-high-lift"},
    };

    for (const auto& condition : conditions)
    {
        const auto& drag_coeff     = std::get<0>(condition);
        const auto& lift_coeff     = std::get<1>(condition);
        const auto& wind_velocity  = std::get<2>(condition);
        const auto& condition_name = std::get<3>(condition);

        const std::string name{condition_name};

        const timer::Timer t(name);

        SimpleEngine engine(drag_coeff, lift_coeff, wind_velocity);
        engine.initializeScene();

        auto alembic_manager =
            elasty::createClothAlembicManager(name + ".abc", engine.m_cloth_sim_object, engine.getDeltaFrameTime());

        for (unsigned int frame = 0; frame < 300; ++frame)
        {
            alembic_manager->submitCurrentStatus();
            engine.proceedFrame();
        }
    }

    return 0;
}
