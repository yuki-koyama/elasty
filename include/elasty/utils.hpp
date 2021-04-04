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

    /// \brief An abstract class for making alembic export easier.
    class AlembicManager
    {
    public:
        /// \brief Submit the current frame to the managed alembic file.
        virtual void submitCurrentStatus() = 0;
    };

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

    /// \brief A utility function to instantiating an alembic manager tailered for cloth objects
    std::shared_ptr<AlembicManager> createClothAlembicManager(const std::string&                    output_file_path,
                                                              const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                                              const double                          delta_time);

    /// \brief Export the current cloth state as an obj file.
    void exportCurrentClothStateAsObj(const std::string&                    output_file_path,
                                      const std::shared_ptr<ClothSimObject> cloth_sim_object);

    /// \brief Generate rectangular plane mesh for cloth simulation in the obj string format.
    ///
    /// \param vertical_resolution Need to be an even number.
    std::string generateClothMeshObjData(const double   width                 = 2.0,
                                         const double   height                = 2.0,
                                         const unsigned horizontal_resolution = 50,
                                         const unsigned vertical_resolution   = 50);
} // namespace elasty

#endif // ELASTY_UTILS_HPP
