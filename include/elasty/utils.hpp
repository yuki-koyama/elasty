#ifndef elasty_utils_hpp
#define elasty_utils_hpp

#include <vector>
#include <Eigen/Core>

namespace elasty
{
    struct Particle;
    class Constraint;

    void generateFixedPointConstraints(const Eigen::Vector3d& search_position,
                                       const Eigen::Vector3d& fixed_position,
                                       const std::vector<std::shared_ptr<Particle>>& particles,
                                       std::vector<std::shared_ptr<Constraint>>& constraints);
}

#endif /* elasty_utils_hpp */
