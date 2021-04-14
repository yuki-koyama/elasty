#ifndef ELASTY_ALEMBIC_MANAGER_HPP
#define ELASTY_ALEMBIC_MANAGER_HPP

#include <Eigen/Core>
#include <cstddef>
#include <memory>
#include <string>

namespace elasty
{
    class ClothSimObject;

    /// \brief An abstract class for making alembic export easier.
    class AbstractAlembicManager
    {
    public:
        /// \brief Submit the current frame to the managed alembic file.
        virtual void submitCurrentStatus() = 0;
    };

    /// \brief A utility function to instantiating an alembic manager tailered for cloth objects
    std::shared_ptr<AbstractAlembicManager>
    createClothAlembicManager(const std::string&                    output_file_path,
                              const std::shared_ptr<ClothSimObject> cloth_sim_object,
                              const double                          delta_time);

    /// \param positions Positions of the vertices in the format: [x_{0}, y_{0}, ..., x_{n - 1}, y_{n - 1}].
    std::shared_ptr<elasty::AbstractAlembicManager> createTriangleMesh2dAlembicManager(const std::string&  file_path,
                                                                                       const double        delta_time,
                                                                                       const std::size_t   num_verts,
                                                                                       const std::size_t   num_elems,
                                                                                       const double*       positions,
                                                                                       const std::int32_t* indices);

    std::shared_ptr<elasty::AbstractAlembicManager> createTetraMeshAlembicManager(const std::string&  file_path,
                                                                                  const double        delta_time,
                                                                                  const std::size_t   num_verts,
                                                                                  const std::size_t   num_elems,
                                                                                  const double*       positions,
                                                                                  const std::int32_t* indices);
} // namespace elasty

#endif // ELASTY_ALEMBIC_MANAGER_HPP
