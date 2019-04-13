#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>

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
    const Eigen::Vector3d n = (x_0 - x_1).normalized();

    assert(!n.hasNaN());

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
