#ifndef ELASTY_FEM_HPP
#define ELASTY_FEM_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>

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
    /// \brief Extract the rotational part of the given square matrix by performing polar decomposion.
    ///
    /// \details This implementation is based on SVD. It checks the determinant to avoid any reflection.
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>
    extractRotation(const Eigen::MatrixBase<Derived>& F)
    {
        const auto svd   = F.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
        const auto Sigma = svd.singularValues();
        const auto U     = svd.matrixU();
        const auto V     = svd.matrixV();
        const auto R     = (U * V.transpose()).eval();

        assert(std::abs(std::abs(R.determinant()) - 1.0) < 1e-04);

        if constexpr (Derived::RowsAtCompileTime == 2)
        {
            // Just ignore reflection (if any)
            // TODO: Discuss whether this is a good strategy, or not
            return R;
        }
        else if constexpr (Derived::RowsAtCompileTime == 3)
        {
            // Correct reflection (if any)
            return (R.determinant() > 0) ? R : -R;
        }
    }

    /// \brief Calculate the first Lame parameter, $\lambda$.
    ///
    /// \details Reference: [1]
    template <typename Scalar> constexpr Scalar calcFirstLame(const Scalar youngs_modulus, const Scalar poisson_ratio)
    {
        return youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    }

    /// \brief Calculate the second Lame parameter, $\mu$.
    ///
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

    /// \brief Calculate either the "deformed" shape matrix (D_s) or "reference" shape matrix (D_m) of a tetrahedron.
    ///
    /// \param x_0 A 3D vector.
    ///
    /// \param x_1 A 3D vector.
    ///
    /// \param x_2 A 3D vector.
    ///
    /// \param x_3 A 3D vector.
    ///
    /// \details Reference: [1]
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 3> calc3dShapeMatrix(const Eigen::MatrixBase<Derived>& x_0,
                                                                    const Eigen::MatrixBase<Derived>& x_1,
                                                                    const Eigen::MatrixBase<Derived>& x_2,
                                                                    const Eigen::MatrixBase<Derived>& x_3)
    {
        using Mat = Eigen::Matrix<typename Derived::Scalar, 3, 3>;

        Mat shape_matrix;
        shape_matrix.col(0) = x_1 - x_0;
        shape_matrix.col(1) = x_2 - x_0;
        shape_matrix.col(2) = x_3 - x_0;

        return shape_matrix;
    }

    /// \brief Calculate the area of a triangle in 2D
    template <typename Derived>
    typename Derived::Scalar calc2dTriangleArea(const Eigen::MatrixBase<Derived>& x_0,
                                                const Eigen::MatrixBase<Derived>& x_1,
                                                const Eigen::MatrixBase<Derived>& x_2)
    {
        const auto r_1 = x_1 - x_0;
        const auto r_2 = x_2 - x_0;

        return 0.5 * std::abs(r_1(0) * r_2(1) - r_2(0) * r_1(1));
    }

    /// \brief Calculate the volume of a tetrahedron
    template <typename Derived>
    typename Derived::Scalar calcTetrahedronVolume(const Eigen::MatrixBase<Derived>& x_0,
                                                   const Eigen::MatrixBase<Derived>& x_1,
                                                   const Eigen::MatrixBase<Derived>& x_2,
                                                   const Eigen::MatrixBase<Derived>& x_3)
    {
        const auto r_1 = x_1 - x_0;
        const auto r_2 = x_2 - x_0;
        const auto r_3 = x_3 - x_0;

        return std::abs(r_1.dot(r_2.cross(r_3))) / 6.0;
    }

    /// \brief Calculate the diagonal elements of the lumped mass matrix.
    ///
    /// \details This function takes the "barycentric" approach. See https://www.alecjacobson.com/weblog/?p=1146 .
    template <typename DerivedV, typename DerivedF>
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1>
    calcLumpedMasses(const Eigen::MatrixBase<DerivedV>& verts,
                     const Eigen::MatrixBase<DerivedF>& elems,
                     const typename DerivedV::Scalar    total_mass)
    {
        using Scalar = typename DerivedV::Scalar;
        using Vec    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        assert(verts.cols() == 1);
        assert(verts.size() % 2 == 0);
        assert(elems.rows() == 3);

        const auto num_verts = verts.size() / 2;

        Scalar total_area = 0;
        Vec    masses     = Vec::Zero(verts.size());

        for (std::size_t elem_index = 0; elem_index < elems.cols(); ++elem_index)
        {
            const auto& indices = elems.col(elem_index);

            const Scalar area = calc2dTriangleArea(verts.template segment<2>(2 * indices[0]),
                                                   verts.template segment<2>(2 * indices[1]),
                                                   verts.template segment<2>(2 * indices[2]));

            const Scalar one_third_area = (1.0 / 3.0) * area;

            masses(2 * indices[0] + 0) += one_third_area;
            masses(2 * indices[0] + 1) += one_third_area;
            masses(2 * indices[1] + 0) += one_third_area;
            masses(2 * indices[1] + 1) += one_third_area;
            masses(2 * indices[2] + 0) += one_third_area;
            masses(2 * indices[2] + 1) += one_third_area;

            total_area += area;
        }

        assert(total_mass > 0);
        assert(total_area > 0);

        return (total_mass / total_area) * masses;
    }

    /// \brief Calculate the Green strain tensor for a finite element.
    ///
    /// \param deform_grad The deformation gradient matrix, which should be either 2-by-2 (2D element in 2D), 3-by-3 (3D
    /// element in 3D), or 3-by-2 (2D element in 3D).
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, Derived::ColsAtCompileTime, Derived::ColsAtCompileTime>
    calcGreenStrain(const Eigen::MatrixBase<Derived>& deform_grad)
    {
        using Mat = Eigen::Matrix<typename Derived::Scalar, Derived::ColsAtCompileTime, Derived::ColsAtCompileTime>;

        return 0.5 * (deform_grad.transpose() * deform_grad - Mat::Identity());
    }

    /// \details Eq. 3.4 in [1]
    template <typename Derived>
    typename Derived::Scalar calcCoRotationalEnergyDensity(const Eigen::MatrixBase<Derived>& deform_grad,
                                                           const typename Derived::Scalar    first_lame,
                                                           const typename Derived::Scalar    second_lame)
    {
        using Mat = Eigen::Matrix<typename Derived::Scalar, Derived::ColsAtCompileTime, Derived::ColsAtCompileTime>;

        const auto R     = extractRotation(deform_grad);
        const auto S     = R.transpose() * deform_grad;
        const auto I     = Mat::Identity();
        const auto trace = (S - I).trace();

        assert(deform_grad.isApprox(R * S));

        return second_lame * (deform_grad - R).squaredNorm() + 0.5 * first_lame * trace * trace;
    }

    /// \details Eq. 3.5 in [1]
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>
    calcCoRotationalPiolaStress(const Eigen::MatrixBase<Derived>& deform_grad,
                                const typename Derived::Scalar    first_lame,
                                const typename Derived::Scalar    second_lame)
    {
        using Mat = Eigen::Matrix<typename Derived::Scalar, Derived::ColsAtCompileTime, Derived::ColsAtCompileTime>;

        const auto R     = extractRotation(deform_grad);
        const auto S     = R.transpose() * deform_grad;
        const auto I     = Mat::Identity();
        const auto trace = (S - I).trace();

        assert(deform_grad.isApprox(R * S));

        return 2.0 * second_lame * (deform_grad - R) + first_lame * trace * R;
    }

    /// \details The first equation in Sec. 3.3 in [1]
    template <typename Derived>
    typename Derived::Scalar calcStVenantKirchhoffEnergyDensity(const Eigen::MatrixBase<Derived>& deform_grad,
                                                                const typename Derived::Scalar    first_lame,
                                                                const typename Derived::Scalar    second_lame)
    {
        const auto E     = calcGreenStrain(deform_grad);
        const auto trace = E.trace();

        return second_lame * E.squaredNorm() + 0.5 * first_lame * trace * trace;
    }

    /// \details Eq. 3.3 in [1]
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>
    calcStVenantKirchhoffPiolaStress(const Eigen::MatrixBase<Derived>& deform_grad,
                                     const typename Derived::Scalar    first_lame,
                                     const typename Derived::Scalar    second_lame)
    {
        const auto E = calcGreenStrain(deform_grad);

        return 2.0 * second_lame * deform_grad * E + first_lame * E.trace() * deform_grad;
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

    template <typename DerivedVec, typename DerivedMat>
    Eigen::Matrix<typename DerivedVec::Scalar, 3, 3>
    calcTetrahedronDeformGrad(const Eigen::MatrixBase<DerivedVec>& x_0,
                              const Eigen::MatrixBase<DerivedVec>& x_1,
                              const Eigen::MatrixBase<DerivedVec>& x_2,
                              const Eigen::MatrixBase<DerivedVec>& x_3,
                              const Eigen::MatrixBase<DerivedMat>& rest_shape_mat_inv)
    {
        const auto D_s = elasty::fem::calc3dShapeMatrix(x_0, x_1, x_2, x_3);
        const auto F   = D_s * rest_shape_mat_inv;

        return F;
    }

    /// \brief Calculate analytic partial derivatives of the deformation gradient $\mathbf{F}$ with respect to the
    /// vertex positions $\mathbf{x}$ (i.e., $\frac{\partial \mathbf{F}}{\partial \mathbf{x}}$) and return it in the
    /// "flattened" format.
    ///
    /// \details The result is a 2-by-2-by-6 third-order tensor but flattened into a 4-by-6 matrix.
    ///
    /// Reference: [2, Appendix D]
    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 4, 6>
    calcFlattenedPartDeformGradPartPos(const Eigen::MatrixBase<Derived>& rest_shape_mat_inv)
    {
        using Scalar = typename Derived::Scalar;

        // Analytics partial derivatives $\frac{\partial \mathbf{D}_\text{s}{\partial x_{i}}$
        Eigen::Matrix<Scalar, 2, 2> PDsPx[6];

        // x[0] (= x_0)
        PDsPx[0] << -1.0, -1.0, 0.0, 0.0;

        // x[1] (= y_0)
        PDsPx[1] << 0.0, 0.0, -1.0, -1.0;

        // x[2] (= x_1)
        PDsPx[2] << 1.0, 0.0, 0.0, 0.0;

        // x[3] (= y_1)
        PDsPx[3] << 0.0, 0.0, 1.0, 0.0;

        // x[4] (= x_2)
        PDsPx[4] << 0.0, 1.0, 0.0, 0.0;

        // x[5] (= y_2)
        PDsPx[5] << 0.0, 0.0, 0.0, 1.0;

        Eigen::Matrix<Scalar, 4, 6> vec_PFPx;
        for (std::size_t i = 0; i < 6; ++i)
        {
            const Eigen::Matrix<Scalar, 2, 2> PFPx_i = PDsPx[i] * rest_shape_mat_inv;

            vec_PFPx.col(i) = Eigen::Map<const Eigen::Matrix<Scalar, 4, 1>>(PFPx_i.data(), PFPx_i.size());
        }

        return vec_PFPx;
    }
} // namespace elasty::fem

#endif // ELASTY_FEM_HPP
