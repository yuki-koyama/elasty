#include <elasty/constraint.hpp>
#include <elasty/engine.hpp>

double elasty::EnvironmentalCollisionConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_engine->m_particles[m_indices[0]].p;
    return m_n.transpose() * x - m_d;
}
