#ifndef elasty_utils_hpp
#define elasty_utils_hpp

#include <string>
#include <vector>
#include <Eigen/Core>

namespace elasty
{
    struct Particle;
    class Constraint;
    class ClothSimObject;
    class AlembicManager;

    void generateFixedPointConstraints(const Eigen::Vector3d& search_position,
                                       const Eigen::Vector3d& fixed_position,
                                       const std::vector<std::shared_ptr<Particle>>& particles,
                                       std::vector<std::shared_ptr<Constraint>>& constraints);

    std::shared_ptr<AlembicManager> createAlembicManager(const std::string& file_path,
                                                         const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                                         const double dt);

    void submitCurrentStatus(const std::shared_ptr<AlembicManager> alembic_manager);
}

#endif /* elasty_utils_hpp */
