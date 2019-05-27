#ifndef cloth_sim_object_hpp
#define cloth_sim_object_hpp

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <elasty/sim-object.hpp>
#include <string>

namespace elasty
{
    class ClothSimObject : public SimObject
    {
    public:
        enum class InPlaneStrategy
        {
            EdgeDistance,
            ContinuumTriangle,
            Both,
        };

        enum class OutOfPlaneStrategy
        {
            Bending,
            IsometricBending,
            Cross,
        };

        using TriangleList =
            Eigen::Matrix<int32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
        using UvList =
            Eigen::Matrix<float, Eigen::Dynamic, 2 * 3, Eigen::RowMajor>;

        ClothSimObject(
            const std::string&       obj_path,
            const double             in_plane_stiffness      = 1.000, ///< PBD
            const double             in_plane_compliance     = 0.001, ///< XPBD
            const double             out_of_plane_stiffness  = 0.100, ///< PBD
            const double             out_of_plane_compliance = 0.010, ///< XPBD
            const double             dt        = 1.0 / 60.0,          ///< XPBD
            const Eigen::Affine3d&   transform = Eigen::Affine3d::Identity(),
            const InPlaneStrategy    in_plane_strategy = InPlaneStrategy::Both,
            const OutOfPlaneStrategy out_of_plane_strategy =
                OutOfPlaneStrategy::IsometricBending);

        TriangleList m_triangle_list;
        UvList       m_uv_list;

        bool hasUv() const { return m_uv_list.size() != 0; }
    };
} // namespace elasty

#endif /* cloth_sim_object_hpp */
