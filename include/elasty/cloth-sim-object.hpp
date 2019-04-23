#ifndef cloth_sim_object_hpp
#define cloth_sim_object_hpp

#include <elasty/sim-object.hpp>
#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace elasty
{
    class ClothSimObject : public SimObject
    {
    public:

        enum class Strategy
        {
            Bending,
            IsometricBending,
            Cross,
        };

        using TriangleList = Eigen::Matrix<int32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

        ClothSimObject(const std::string& obj_path,
                       const double distance_stiffness = 0.90,
                       const double bending_stiffness = 0.50,
                       const Eigen::Affine3d& transform = Eigen::Affine3d::Identity(),
                       const Strategy strategy = Strategy::IsometricBending);

        TriangleList m_triangle_list;
    };
}

#endif /* cloth_sim_object_hpp */
