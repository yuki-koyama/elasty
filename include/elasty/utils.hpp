#ifndef ELASTY_UTILS_HPP
#define ELASTY_UTILS_HPP

#include <Eigen/Core>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
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

    inline std::string generateClothMeshObjData(const double   width                 = 1.0,
                                                const double   height                = 1.0,
                                                const unsigned horizontal_resolution = 10,
                                                const unsigned vertical_resolution   = 10)
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3i> triangles;

        const double h_step = width / static_cast<double>(horizontal_resolution + 1);
        const double v_step = height / static_cast<double>(vertical_resolution + 1);

        // Vertices
        for (unsigned h_index = 0; h_index <= horizontal_resolution; ++h_index)
        {
            for (unsigned v_index = 0; v_index <= vertical_resolution; ++v_index)
            {
                const double x = [&]() {
                    if (v_index % 2 == 0 || h_index == 0)
                    {
                        return static_cast<double>(h_index) * h_step - 0.5 * width;
                    }
                    else
                    {
                        return (static_cast<double>(h_index) - 0.5) * h_step - 0.5 * width;
                    }
                }();
                const double y = static_cast<double>(v_index) * v_step - 0.5 * height;

                vertices.push_back(Eigen::Vector3d(x, y, 0.0));

                // Additional vetex at the even-indexed row
                if (v_index % 2 == 0 || h_index == horizontal_resolution)
                {
                    vertices.push_back(Eigen::Vector3d(0.5 * width, y, 0.0));
                }
            }
        }

        // Triangles
        // TODO

        // Output
        std::stringstream sstream;
        for (const auto& vertex : vertices)
        {
            sstream << "v " << vertex(0) << " " << vertex(1) << " " << vertex(2) << std::endl;
        }
        for (const auto& triangle : triangles)
        {
            sstream << "f " << triangle(0) << " " << triangle(1) << " " << triangle(2) << std::endl;
        }

        return sstream.str();
    }

} // namespace elasty

#endif // ELASTY_UTILS_HPP
