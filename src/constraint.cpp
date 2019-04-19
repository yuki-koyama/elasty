#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>
#include <elasty/particle.hpp>
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

    inline double calculateCotTheta(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
    {
        const double cos_theta = x.dot(y);
        const double sin_theta = x.cross(y).norm();
        return cos_theta / sin_theta;
    }
}

elasty::BendingConstraint::BendingConstraint(const std::shared_ptr<Particle> p_0,
                                             const std::shared_ptr<Particle> p_1,
                                             const std::shared_ptr<Particle> p_2,
                                             const std::shared_ptr<Particle> p_3,
                                             const double stiffness,
                                             const double dihedral_angle) :
Constraint(std::vector<std::shared_ptr<Particle>>{ p_0, p_1, p_2, p_3 }, stiffness),
m_dihedral_angle(dihedral_angle)
{
}

double elasty::BendingConstraint::calculateValue()
{
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;
    const Eigen::Vector3d& x_2 = m_particles[2]->p;
    const Eigen::Vector3d& x_3 = m_particles[3]->p;

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
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;
    const Eigen::Vector3d& x_2 = m_particles[2]->p;
    const Eigen::Vector3d& x_3 = m_particles[3]->p;

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
    if (std::abs(d) - 1.0 < epsilon) { return Eigen::VectorXd::Zero(3 * 4); }

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

elasty::DistanceConstraint::DistanceConstraint(const std::shared_ptr<Particle> p_0,
                                               const std::shared_ptr<Particle> p_1,
                                               const double stiffness,
                                               const double d) :
Constraint(std::vector<std::shared_ptr<Particle>>{ p_0, p_1 }, stiffness),
m_d(d)
{
    assert(d >= 0.0);
}

double elasty::DistanceConstraint::calculateValue()
{
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;

    return (x_0 - x_1).norm() - m_d;
}

Eigen::VectorXd elasty::DistanceConstraint::calculateGrad()
{
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;

    Eigen::Vector3d n = (x_0 - x_1).normalized();

    if (n.hasNaN()) { n = Eigen::Vector3d::Random(3).normalized(); }

    Eigen::VectorXd grad_C = Eigen::VectorXd(3 * 2);
    grad_C.segment<3>(0) = n;
    grad_C.segment<3>(3) = - n;

    return grad_C;
}

elasty::EnvironmentalCollisionConstraint::EnvironmentalCollisionConstraint(const std::shared_ptr<Particle> p_0,
                                                                           const double stiffness,
                                                                           const Eigen::Vector3d& n,
                                                                           const double d) :
Constraint(std::vector<std::shared_ptr<Particle>>{ p_0 }, stiffness),
m_n(n),
m_d(d)
{
}

double elasty::EnvironmentalCollisionConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_particles[0]->p;
    return m_n.transpose() * x - m_d;
}

Eigen::VectorXd elasty::EnvironmentalCollisionConstraint::calculateGrad()
{
    return m_n;
}

elasty::FixedPointConstraint::FixedPointConstraint(const std::shared_ptr<Particle> p_0,
                                                   const double stiffness,
                                                   const Eigen::Vector3d& point) :
Constraint(std::vector<std::shared_ptr<Particle>>{ p_0 }, stiffness),
m_point(point)
{
}

double elasty::FixedPointConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_particles[0]->p;
    return (x - m_point).norm();
}

Eigen::VectorXd elasty::FixedPointConstraint::calculateGrad()
{
    const Eigen::Vector3d& x = m_particles[0]->p;
    const Eigen::Vector3d n = (x - m_point).normalized();

    if (n.hasNaN()) { return Eigen::Vector3d::Zero(); }

    return n;
}

elasty::IsometricBendingConstraint::IsometricBendingConstraint(const std::shared_ptr<Particle> p_0,
                                                               const std::shared_ptr<Particle> p_1,
                                                               const std::shared_ptr<Particle> p_2,
                                                               const std::shared_ptr<Particle> p_3,
                                                               const double stiffness) :
Constraint(std::vector<std::shared_ptr<Particle>>{ p_0, p_1, p_2, p_3 }, stiffness)
{
    const Eigen::Vector3d& x_0 = p_0->x;
    const Eigen::Vector3d& x_1 = p_1->x;
    const Eigen::Vector3d& x_2 = p_2->x;
    const Eigen::Vector3d& x_3 = p_3->x;

    const Eigen::Vector3d e0 = x_1 - x_0;
    const Eigen::Vector3d e1 = x_2 - x_1;
    const Eigen::Vector3d e2 = x_0 - x_2;
    const Eigen::Vector3d e3 = x_3 - x_0;
    const Eigen::Vector3d e4 = x_1 - x_3;

    const double cot_01 = calculateCotTheta(e0, - e1);
    const double cot_02 = calculateCotTheta(e0, - e2);
    const double cot_03 = calculateCotTheta(e0, e3);
    const double cot_04 = calculateCotTheta(e0, e4);

    const Eigen::Vector4d K = Eigen::Vector4d(cot_01 + cot_04, cot_02 + cot_03, - cot_01 - cot_02, - cot_03 - cot_04);

    const double A_0 = 0.5 * e0.cross(e1).norm();
    const double A_1 = 0.5 * e0.cross(e3).norm();

    m_Q = (3.0 / (A_0 + A_1)) * K * K.transpose();
}

double elasty::IsometricBendingConstraint::calculateValue()
{
    double sum = 0.0;
    for (unsigned int i = 0; i < 4; ++ i)
    {
        for (unsigned int j = 0; j < 4; ++ j)
        {
            sum += m_Q(i, j) * double(m_particles[i]->p.transpose() * m_particles[j]->p);
        }
    }
    return 0.5 * sum;
}

Eigen::VectorXd elasty::IsometricBendingConstraint::calculateGrad()
{
    Eigen::VectorXd grad_C = Eigen::VectorXd(3 * 4);

    for (unsigned int i = 0; i < 4; ++ i)
    {
        Eigen::Vector3d sum = Eigen::Vector3d::Zero();
        for (unsigned int j = 0; j < 4; ++ j)
        {
            sum += m_Q(i, j) * m_particles[j]->p;
        }
        grad_C.segment<3>(3 * i) = sum;
    }

    return grad_C;
}
