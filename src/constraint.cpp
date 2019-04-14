#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <Eigen/Geometry>

namespace
{
    inline Eigen::Matrix3d convert_vector_to_cross_operator(const Eigen::Vector3d& vec)
    {
        Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
        mat(0, 1) = + vec(2);
        mat(0, 2) = - vec(1);
        mat(1, 0) = - vec(2);
        mat(1, 2) = + vec(0);
        mat(2, 0) = + vec(1);
        mat(2, 1) = - vec(0);
        return mat;
    };
}

double elasty::BendingConstraint::calculateValue()
{
    const Eigen::Vector3d& x_0 = m_engine->m_particles[m_indices[0]].p;
    const Eigen::Vector3d& x_1 = m_engine->m_particles[m_indices[1]].p;
    const Eigen::Vector3d& x_2 = m_engine->m_particles[m_indices[2]].p;
    const Eigen::Vector3d& x_3 = m_engine->m_particles[m_indices[3]].p;

    const Eigen::Vector3d p_10 = x_1 - x_0;
    const Eigen::Vector3d p_20 = x_2 - x_0;
    const Eigen::Vector3d p_30 = x_3 - x_0;

    const Eigen::Vector3d n_0 = p_10.cross(p_20).normalized();
    const Eigen::Vector3d n_1 = p_10.cross(p_30).normalized();

    assert(!n_0.hasNaN());
    assert(!n_1.hasNaN());

    const double current_dihedral_angle = std::acos(std::min(+ 1.0, std::max(- 1.0, n_0.dot(n_1))));

    assert(!std::isnan(current_dihedral_angle));

    return current_dihedral_angle - m_dihedral_angle;
}

Eigen::VectorXd elasty::BendingConstraint::calculateGrad()
{
    const Eigen::Vector3d& x_0 = m_engine->m_particles[m_indices[0]].p;
    const Eigen::Vector3d& x_1 = m_engine->m_particles[m_indices[1]].p;
    const Eigen::Vector3d& x_2 = m_engine->m_particles[m_indices[2]].p;
    const Eigen::Vector3d& x_3 = m_engine->m_particles[m_indices[3]].p;

    // Assuming that p_0 = [ 0, 0, 0 ]^T without loss of generality
    const Eigen::Vector3d p_1 = x_1 - x_0;
    const Eigen::Vector3d p_2 = x_2 - x_0;
    const Eigen::Vector3d p_3 = x_3 - x_0;

    const Eigen::Vector3d p_1_cross_p_2 = p_1.cross(p_2);
    const Eigen::Vector3d p_1_cross_p_3 = p_1.cross(p_3);

    const Eigen::Vector3d n_0 = p_1_cross_p_2.normalized();
    const Eigen::Vector3d n_1 = p_1_cross_p_3.normalized();

    const double d = n_0.dot(n_1);

    constexpr double epsilon = 1e-12;
    if (std::abs(d) - 1.0 < epsilon) { return Eigen::VectorXd::Random(3 * 4); }

    const double common_coeff = - 1.0 / std::sqrt(1.0 - d * d);

    auto calculate_gradient_of_normalized_cross_product_wrt_p_1 = [](const Eigen::Vector3d& p_1,
                                                                     const Eigen::Vector3d& p_2,
                                                                     const Eigen::Vector3d& n)
    -> Eigen::Matrix3d
    {
        return + (1.0 / p_1.cross(p_2).norm()) * (- convert_vector_to_cross_operator(p_2) + n * (n.cross(p_2)).transpose());
    };

    auto calculate_gradient_of_normalized_cross_product_wrt_p_2 = [](const Eigen::Vector3d& p_1,
                                                                     const Eigen::Vector3d& p_2,
                                                                     const Eigen::Vector3d& n)
    -> Eigen::Matrix3d
    {
        return - (1.0 / p_1.cross(p_2).norm()) * (- convert_vector_to_cross_operator(p_1) + n * (n.cross(p_1)).transpose());
    };

    const Eigen::Vector3d grad_C_wrt_p_1 = common_coeff * (calculate_gradient_of_normalized_cross_product_wrt_p_1(n_0, p_1, p_2).transpose() * n_1 + calculate_gradient_of_normalized_cross_product_wrt_p_1(n_1, p_1, p_3).transpose() * n_0);
    const Eigen::Vector3d grad_C_wrt_p_2 = common_coeff * calculate_gradient_of_normalized_cross_product_wrt_p_2(n_0, p_1, p_2).transpose() * n_1;
    const Eigen::Vector3d grad_C_wrt_p_3 = common_coeff * calculate_gradient_of_normalized_cross_product_wrt_p_2(n_1, p_1, p_3).transpose() * n_0;
    const Eigen::Vector3d grad_C_wrt_p_0 = - grad_C_wrt_p_1 - grad_C_wrt_p_2 - grad_C_wrt_p_3;

    Eigen::VectorXd grad_C = Eigen::VectorXd(3 * 4);
    grad_C.segment<3>(3 * 0) = grad_C_wrt_p_0;
    grad_C.segment<3>(3 * 1) = grad_C_wrt_p_1;
    grad_C.segment<3>(3 * 2) = grad_C_wrt_p_2;
    grad_C.segment<3>(3 * 3) = grad_C_wrt_p_3;

    assert(!grad_C.hasNaN());

    return grad_C;
}

double elasty::DistanceConstraint::calculateValue()
{
    const Eigen::Vector3d& x_0 = m_engine->m_particles[m_indices[0]].p;
    const Eigen::Vector3d& x_1 = m_engine->m_particles[m_indices[1]].p;

    return (x_0 - x_1).norm() - m_d;
}

Eigen::VectorXd elasty::DistanceConstraint::calculateGrad()
{
    const Eigen::Vector3d& x_0 = m_engine->m_particles[m_indices[0]].p;
    const Eigen::Vector3d& x_1 = m_engine->m_particles[m_indices[1]].p;

    Eigen::Vector3d n = (x_0 - x_1).normalized();

    if (n.hasNaN()) { n = Eigen::Vector3d::Random(3).normalized(); }

    Eigen::VectorXd grad_C = Eigen::VectorXd(3 * 2);
    grad_C.segment<3>(0) = n;
    grad_C.segment<3>(3) = - n;

    return grad_C;
}

double elasty::EnvironmentalCollisionConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_engine->m_particles[m_indices[0]].p;
    return m_n.transpose() * x - m_d;
}

Eigen::VectorXd elasty::EnvironmentalCollisionConstraint::calculateGrad()
{
    return m_n;
}

elasty::FixedPointConstraint::FixedPointConstraint(const Engine* engine,
                                                   const std::vector<unsigned int>& indices,
                                                   const double stiffness,
                                                   const Eigen::Vector3d& point) :
Constraint(engine, indices, stiffness),
m_point(point)
{
    assert(indices.size() == 1);
}

double elasty::FixedPointConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_engine->m_particles[m_indices[0]].p;
    return (x - m_point).norm();
}

Eigen::VectorXd elasty::FixedPointConstraint::calculateGrad()
{
    const Eigen::Vector3d& x = m_engine->m_particles[m_indices[0]].p;
    const Eigen::Vector3d n = (x - m_point).normalized();

    if (n.hasNaN()) { return Eigen::Vector3d::Zero(); }

    return n;
}
