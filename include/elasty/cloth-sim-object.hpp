#ifndef cloth_sim_object_hpp
#define cloth_sim_object_hpp

#include <elasty/sim-object.hpp>
#include <string>
#include <Eigen/Core>

namespace elasty
{
    class ClothSimObject : public SimObject
    {
    public:

        using TriangleList = Eigen::Matrix<int32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;

        ClothSimObject(const std::string& obj_path);

        TriangleList m_triangle_indices;
    };
}

#endif /* cloth_sim_object_hpp */
