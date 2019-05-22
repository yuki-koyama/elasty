#include <elasty/constraint.hpp>
#include <elasty/particle.hpp>

int main()
{
    auto p_0 = std::make_shared<elasty::Particle>(
        Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), 1.0);
    auto p_1 = std::make_shared<elasty::Particle>(
        Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d::Zero(), 1.0);
    auto p_2 = std::make_shared<elasty::Particle>(
        Eigen::Vector3d(-0.5, 0.5, 0.0), Eigen::Vector3d::Zero(), 1.0);
    auto p_3 = std::make_shared<elasty::Particle>(
        Eigen::Vector3d(+0.5, 0.5, 0.0), Eigen::Vector3d::Zero(), 1.0);

    const Eigen::Vector3d& x_0 = p_0->x;
    const Eigen::Vector3d& x_1 = p_1->x;
    const Eigen::Vector3d& x_2 = p_2->x;
    const Eigen::Vector3d& x_3 = p_3->x;

    const Eigen::Vector3d p_10 = x_1 - x_0;
    const Eigen::Vector3d p_20 = x_2 - x_0;
    const Eigen::Vector3d p_30 = x_3 - x_0;

    const Eigen::Vector3d n_0 = p_10.cross(p_20).normalized();
    const Eigen::Vector3d n_1 = p_10.cross(p_30).normalized();

    assert(!n_0.hasNaN());
    assert(!n_1.hasNaN());

    const double dihedral_angle =
        std::acos(std::max(-1.0, std::min(+1.0, n_0.dot(n_1))));

    assert(!std::isnan(dihedral_angle));

    const auto constraint = std::make_shared<elasty::BendingConstraint>(
        p_0, p_1, p_2, p_3, 1.0, dihedral_angle);

    const double value = constraint->calculateValue();

    Eigen::Matrix<double, 12, 1> grad;
    constraint->calculateGrad(grad.data());

    constexpr double epsilon = 1e-20;
    if (!(std::abs(value) < epsilon) || !(grad.norm() < epsilon))
    {
        throw std::runtime_error("");
    }

    return 0;
}
