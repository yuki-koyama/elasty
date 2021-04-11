#ifndef ELASTY_CLOTH_SIM_OBJECT_HPP
#define ELASTY_CLOTH_SIM_OBJECT_HPP

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

        using TriangleList = Eigen::Matrix<std::int32_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
        using UvList       = Eigen::Matrix<float, Eigen::Dynamic, 2 * 3, Eigen::RowMajor>;

        ClothSimObject(const unsigned           resolution,
                       const double             in_plane_stiffness      = 1.000,      ///< PBD
                       const double             in_plane_compliance     = 0.001,      ///< XPBD
                       const double             out_of_plane_stiffness  = 0.100,      ///< PBD
                       const double             out_of_plane_compliance = 0.010,      ///< XPBD
                       const double             dt                      = 1.0 / 60.0, ///< XPBD
                       const Eigen::Affine3d&   transform               = Eigen::Affine3d::Identity(),
                       const InPlaneStrategy    in_plane_strategy       = InPlaneStrategy::Both,
                       const OutOfPlaneStrategy out_of_plane_strategy   = OutOfPlaneStrategy::IsometricBending);

        TriangleList m_triangle_list;
        UvList       m_uv_list;

        bool hasUv() const { return m_uv_list.size() != 0; }

        /// \brief Getter of the areas of the triangles.
        const Eigen::VectorXd& getAreaList() const { return m_area_list; }

        /// \brief Apply aerodynamic forces to the relevant particles.
        ///
        /// \details The aerodynamics model is based on [Wilson+14]. This member function is intended to be called in
        /// the setExternalForces() routine in engines.
        ///
        /// \param drag_coeff The coefficient for drag forces. This value should be larger than lift_coeff.
        ///
        /// \param lift_coeff The coefficient for lift forces. This value should be smaller than drag_coeff.
        void applyAerodynamicForces(const Eigen::Vector3d& global_velocity = Eigen::Vector3d::Zero(),
                                    const double           drag_coeff      = 0.100,
                                    const double           lift_coeff      = 0.060);

    private:
        Eigen::VectorXd m_area_list;

        /// \brief Calculate areas of the triangles.
        ///
        /// \details The result will be stored in m_area_list. This function needs to be called only once.
        void calculateAreas();
    };
} // namespace elasty

#endif // ELASTY_CLOTH_SIM_OBJECT_HPP
