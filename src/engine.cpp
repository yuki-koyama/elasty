#include <elasty/engine.hpp>
#include <elasty/constraint.hpp>

void elasty::Engine::stepTime()
{
    // Apply external forces
    setExternalForces();
    for (auto& particle : m_particles)
    {
        particle.v = particle.v + m_dt * (1.0 / particle.m) * particle.f;
    }

    // Calculate predicted positions
    for (auto& particle : m_particles)
    {
        particle.p = particle.x + m_dt * particle.v;
    }

    // Generate collision constraints
    generateCollisionConstraints();

    // Solve constraints
    for (unsigned int i = 0; i < m_num_iterations; ++ i)
    {
        for (auto constraint : m_constraints)
        {
            projectConstraint(constraint);
        }

        for (auto constraint : m_instant_constraints)
        {
            projectConstraint(constraint);
        }
    }

    // Apply the results
    for (auto& particle : m_particles)
    {
        particle.v = (particle.p - particle.x) * (1.0 / m_dt);
        particle.x = particle.p;
    }

    // Update velocities
    updateVelocities();

    // Clear instant constraints
    m_instant_constraints.clear();
}

void elasty::Engine::addConstraint(std::shared_ptr<Constraint> constraint)
{
    m_constraints.push_back(constraint);
}

void elasty::Engine::addInstantConstraint(std::shared_ptr<Constraint> constraint)
{
    m_instant_constraints.push_back(constraint);
}

void elasty::Engine::projectConstraint(std::shared_ptr<Constraint> constraint)
{
    // Calculate $C$
    const double C = constraint->calculateValue();

    // Skip if the constraint is unilateral and is satisfied
    if (constraint->getType() == ConstraintType::Unilateral && C >= 0.0) { return; }

    // Calculate $\Nabla C$
    const Eigen::VectorXd grad_C = constraint->calculateGrad();

    // Skip if the gradient is sufficiently small
    if (grad_C.isApprox(Eigen::VectorXd::Zero(grad_C.size()))) { return; }

    // Calculate $\mathbf{M}^{-1}$
    const unsigned int n = constraint->m_indices.size();
    std::vector<double> inverse_mass_raw(3 * n);
    for (unsigned int j = 0; j < n; ++ j)
    {
        inverse_mass_raw[j * 3 + 0] = 1.0 / m_particles[constraint->m_indices[j]].m;
        inverse_mass_raw[j * 3 + 1] = 1.0 / m_particles[constraint->m_indices[j]].m;
        inverse_mass_raw[j * 3 + 2] = 1.0 / m_particles[constraint->m_indices[j]].m;
    }
    const auto inverse_M_diagnal = Eigen::Map<Eigen::VectorXd>(inverse_mass_raw.data(), 3 * n);

    // Calculate $s$
    const double s = C / (grad_C.transpose() * inverse_M_diagnal.asDiagonal() * grad_C);

    // Calculate $\Delta x$
    const Eigen::VectorXd delta_x = - s * inverse_M_diagnal.asDiagonal() * grad_C;
    assert(!delta_x.hasNaN());

    // Update predicted positions
    for (unsigned int j = 0; j < n; ++ j)
    {
        m_particles[constraint->m_indices[j]].p += constraint->m_stiffness * delta_x.segment<3>(3 * j);
    }

    assert(std::abs(C) >= std::abs(constraint->calculateValue()));
}
