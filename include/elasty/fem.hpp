#ifndef ELASTY_FEM_HPP
#define ELASTY_FEM_HPP

#include <Eigen/Core>

// References:
//
// [1] Eftychios Sifakis and Jernej Barbic. 2012. FEM simulation of 3D deformable solids: a practitioner's guide to
// theory, discretization and model reduction. In ACM SIGGRAPH 2012 Courses (SIGGRAPH '12). Association for Computing
// Machinery, New York, NY, USA, Article 20, 1–50. DOI:https://doi.org/10.1145/2343483.2343501
//
// [2] Theodore Kim and David Eberle. 2020. Dynamic deformables: implementation and production practicalities. In ACM
// SIGGRAPH 2020 Courses (SIGGRAPH '20). Association for Computing Machinery, New York, NY, USA, Article 23, 1–182.
// DOI:https://doi.org/10.1145/3388769.3407490

namespace elasty::fem
{
    /// \details Reference: [1]
    template <typename Scalar> constexpr Scalar calcFirstLame(const Scalar youngs_modulus, const Scalar poisson_ratio)
    {
        return youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    }

    /// \details Reference: [1]
    template <typename Scalar> constexpr Scalar calcSecondLame(const Scalar youngs_modulus, const Scalar poisson_ratio)
    {
        return youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    }

    /// \brief Calculate either the "deformed" shape matrix (D_s) or "reference" shape matrix (D_m) of a triangle.
    ///
    /// \param x_0 A 2D vector.
    ///
    /// \param x_1 A 2D vector.
    ///
    /// \param x_2 A 2D vector.
    ///
    /// \details Reference: [1]
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 2, 2> calc2dShapeMatrix(const Eigen::MatrixBase<Derived>& x_0,
                                                                    const Eigen::MatrixBase<Derived>& x_1,
                                                                    const Eigen::MatrixBase<Derived>& x_2)
    {
        using Mat = Eigen::Matrix<typename Derived::Scalar, 2, 2>;

        Mat shape_matrix;
        shape_matrix.col(0) = x_1 - x_0;
        shape_matrix.col(1) = x_2 - x_0;

        return shape_matrix;
    }

    /// \brief Calculate the area of a 2D triangle
    template <typename Derived>
    typename Derived::Scalar calc2dTriangleArea(const Eigen::MatrixBase<Derived>& x_0,
                                                const Eigen::MatrixBase<Derived>& x_1,
                                                const Eigen::MatrixBase<Derived>& x_2)
    {
        const auto r_1 = x_1 - x_0;
        const auto r_2 = x_2 - x_0;

        return 0.5 * std::abs(r_1(0) * r_2(1) - r_2(0) * r_1(1));
    }

    /// \brief Calculate the Green strain tensor for a 2D element.
    ///
    /// \param deform_grad The deformation gradient, which is a 2x2 matrix.
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 2, 2> calc2dGreenStrain(const Eigen::MatrixBase<Derived>& deform_grad)
    {
        using Mat = Eigen::Matrix<typename Derived::Scalar, 2, 2>;

        return 0.5 * (deform_grad.transpose() * deform_grad - Mat::Identity());
    }

    template <typename Derived>
    typename Derived::Scalar calcStVenantKirchhoffEnergyDensity(const Eigen::MatrixBase<Derived>& deform_grad,
                                                                const typename Derived::Scalar    first_lame,
                                                                const typename Derived::Scalar    second_lame)
    {
        const auto E     = calc2dGreenStrain(deform_grad);
        const auto trace = E.trace();

        return first_lame * E.squaredNorm() + 0.5 * second_lame * trace * trace;
    }

    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 2, 2>
    calcStVenantKirchhoffPiolaStress(const Eigen::MatrixBase<Derived>& deform_grad,
                                     const typename Derived::Scalar    first_lame,
                                     const typename Derived::Scalar    second_lame)
    {
        const auto E = calc2dGreenStrain(deform_grad);

        return 2.0 * first_lame * deform_grad * E + second_lame * E.trace() * deform_grad;
    }

    template <typename DerivedVec, typename DerivedMat>
    Eigen::Matrix<typename DerivedVec::Scalar, 2, 2>
    calc2dTriangleDeformGrad(const Eigen::MatrixBase<DerivedVec>& x_0,
                             const Eigen::MatrixBase<DerivedVec>& x_1,
                             const Eigen::MatrixBase<DerivedVec>& x_2,
                             const Eigen::MatrixBase<DerivedMat>& rest_shape_mat_inv)
    {
        const auto D_s = elasty::fem::calc2dShapeMatrix(x_0, x_1, x_2);
        const auto F   = D_s * rest_shape_mat_inv;

        return F;
    }
} // namespace elasty::fem

#endif // ELASTY_FEM_HPP
