#ifndef ELASTY_ALEMBIC_MANAGER_HPP
#define ELASTY_ALEMBIC_MANAGER_HPP

#include <Eigen/Core>
#include <memory>
#include <string>

namespace elasty
{
    class ClothSimObject;

    /// \brief An abstract class for making alembic export easier.
    class AlembicManager
    {
    public:
        /// \brief Submit the current frame to the managed alembic file.
        virtual void submitCurrentStatus() = 0;
    };

    /// \brief A utility function to instantiating an alembic manager tailered for cloth objects
    std::shared_ptr<AlembicManager> createClothAlembicManager(const std::string&                    output_file_path,
                                                              const std::shared_ptr<ClothSimObject> cloth_sim_object,
                                                              const double                          delta_time);
} // namespace elasty

#endif // ELASTY_ALEMBIC_MANAGER_HPP
