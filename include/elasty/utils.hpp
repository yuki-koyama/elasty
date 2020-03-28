#ifndef ELASTY_UTILS_HPP
#define ELASTY_UTILS_HPP

#include <Eigen/Core>
#include <memory>
#include <string>
#include <vector>

namespace elasty
{
    struct Particle;
    class AbstractConstraint;
    class ClothSimObject;
    class AlembicManager;

    constexpr double pi() { return 3.14159265358979323846264338327950288; }

    /// \param particles an array of N particles
    ///
    /// \return an array of 3 * N values:
    /// std::vector<float>
    /// {
    ///     p[0].x, p[0].y, p[0].z,
    ///     p[1].x, p[1].y, p[1].z,
    ///     ...
    /// }
    std::vector<float> packParticlePositions(const std::vector<std::shared_ptr<Particle>>& particles);

    void setRandomVelocities(const std::vector<std::shared_ptr<Particle>>& particles, const double scale = 1.0);

    void generateFixedPointConstraints(const Eigen::Vector3d&                            search_position,
                                       const Eigen::Vector3d&                            fixed_position,
                                       const std::vector<std::shared_ptr<Particle>>&     particles,
                                       std::vector<std::shared_ptr<AbstractConstraint>>& constraints);

    std::shared_ptr<AlembicManager> createAlembicManager(const std::string&                    file_path,
                                                         const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                                         const double                          delta_time);

    void submitCurrentStatus(const std::shared_ptr<AlembicManager> alembic_manager);
} // namespace elasty

#endif // ELASTY_UTILS_HPP
