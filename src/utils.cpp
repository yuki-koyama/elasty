#include <elasty/utils.hpp>
#include <elasty/constraint.hpp>
#include <elasty/particle.hpp>

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
