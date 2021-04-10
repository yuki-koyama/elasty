#include <Eigen/LU>
#include <elasty/fem.hpp>
#include <gtest/gtest.h>

namespace
{
    constexpr double k_youngs_modulus = 800.0;
    constexpr double k_poisson_ratio  = 0.40;

    constexpr double k_first_lame  = elasty::fem::calcFirstLame(k_youngs_modulus, k_poisson_ratio);
    constexpr double k_second_lame = elasty::fem::calcSecondLame(k_youngs_modulus, k_poisson_ratio);
} // namespace

TEST(FemTest, StVenantKirchhoff2d)
{
    auto calcPiolaStress = [](const Eigen::Matrix2d& F) -> Eigen::Matrix2d {
        return elasty::fem::calcStVenantKirchhoffPiolaStress(F, k_first_lame, k_second_lame);
    };
    auto calcEnergyDensity = [](const Eigen::Matrix2d& F) -> double {
        return elasty::fem::calcStVenantKirchhoffEnergyDensity(F, k_first_lame, k_second_lame);
    };

    Eigen::MatrixXd x_rest{2, 3};
    x_rest.col(0) = Eigen::Vector2d{0.0, 0.0};
    x_rest.col(1) = Eigen::Vector2d{1.0, 0.0};
    x_rest.col(2) = Eigen::Vector2d{0.0, 1.0};

    Eigen::MatrixXd x{2, 3};
    x.col(0) = Eigen::Vector2d{0.0, 0.5};
    x.col(1) = Eigen::Vector2d{2.0, 0.0};
    x.col(2) = Eigen::Vector2d{0.0, 1.5};

    // Calculate $\frac{\partial \Psi}{\partial \mathbf{x}}$ analytically
    const auto D_m      = elasty::fem::calc2dShapeMatrix(x_rest.col(0), x_rest.col(1), x_rest.col(2));
    const auto D_m_inv  = D_m.inverse();
    const auto F        = elasty::fem::calc2dTriangleDeformGrad(x.col(0), x.col(1), x.col(2), D_m_inv);
    const auto P        = calcPiolaStress(F);
    const auto vec_P    = Eigen::Map<const Eigen::Vector4d>{P.data(), P.size()};
    const auto vec_PFPx = elasty::fem::calcFlattenedPartDeformGradPartPos(D_m_inv);
    const auto PPsiPx   = vec_PFPx.transpose() * vec_P;

    // Calculate $\frac{\partial \Psi}{\partial \mathbf{x}}$ numerically
    Eigen::MatrixXd x_temp = x;
    Eigen::VectorXd diff{6};
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 2; ++j)
        {
            constexpr double eps = 1e-06;

            x_temp(j, i) += eps;
            const auto F_p =
                elasty::fem::calc2dTriangleDeformGrad(x_temp.col(0), x_temp.col(1), x_temp.col(2), D_m_inv);
            const auto e_p = calcEnergyDensity(F_p);

            x_temp = x;

            x_temp(j, i) -= eps;
            const auto F_m =
                elasty::fem::calc2dTriangleDeformGrad(x_temp.col(0), x_temp.col(1), x_temp.col(2), D_m_inv);
            const auto e_m = calcEnergyDensity(F_m);

            x_temp = x;

            diff[i * 2 + j] = (e_p - e_m) / (2.0 * eps);
        }
    }

    EXPECT_TRUE((PPsiPx - diff).norm() < 1e-04);
}

TEST(FemTest, CoRotational2d)
{
    auto calcPiolaStress = [](const Eigen::Matrix2d& F) -> Eigen::Matrix2d {
        return elasty::fem::calcCoRotationalPiolaStress(F, k_first_lame, k_second_lame);
    };
    auto calcEnergyDensity = [](const Eigen::Matrix2d& F) -> double {
        return elasty::fem::calcCoRotationalEnergyDensity(F, k_first_lame, k_second_lame);
    };

    Eigen::MatrixXd x_rest{2, 3};
    x_rest.col(0) = Eigen::Vector2d{0.0, 0.0};
    x_rest.col(1) = Eigen::Vector2d{1.0, 0.0};
    x_rest.col(2) = Eigen::Vector2d{0.0, 1.0};

    Eigen::MatrixXd x{2, 3};
    x.col(0) = Eigen::Vector2d{0.0, 0.5};
    x.col(1) = Eigen::Vector2d{2.0, 0.0};
    x.col(2) = Eigen::Vector2d{0.0, 1.5};

    // Calculate $\frac{\partial \Psi}{\partial \mathbf{x}}$ analytically
    const auto D_m      = elasty::fem::calc2dShapeMatrix(x_rest.col(0), x_rest.col(1), x_rest.col(2));
    const auto D_m_inv  = D_m.inverse();
    const auto F        = elasty::fem::calc2dTriangleDeformGrad(x.col(0), x.col(1), x.col(2), D_m_inv);
    const auto P        = calcPiolaStress(F);
    const auto vec_P    = Eigen::Map<const Eigen::Vector4d>{P.data(), P.size()};
    const auto vec_PFPx = elasty::fem::calcFlattenedPartDeformGradPartPos(D_m_inv);
    const auto PPsiPx   = vec_PFPx.transpose() * vec_P;

    // Calculate $\frac{\partial \Psi}{\partial \mathbf{x}}$ numerically
    Eigen::MatrixXd x_temp = x;
    Eigen::VectorXd diff{6};
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 2; ++j)
        {
            constexpr double eps = 1e-06;

            x_temp(j, i) += eps;
            const auto F_p =
                elasty::fem::calc2dTriangleDeformGrad(x_temp.col(0), x_temp.col(1), x_temp.col(2), D_m_inv);
            const auto e_p = calcEnergyDensity(F_p);

            x_temp = x;

            x_temp(j, i) -= eps;
            const auto F_m =
                elasty::fem::calc2dTriangleDeformGrad(x_temp.col(0), x_temp.col(1), x_temp.col(2), D_m_inv);
            const auto e_m = calcEnergyDensity(F_m);

            x_temp = x;

            diff[i * 2 + j] = (e_p - e_m) / (2.0 * eps);
        }
    }

    EXPECT_TRUE((PPsiPx - diff).norm() < 1e-04);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
