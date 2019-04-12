#include <elasty/elasty.hpp>

void elasty::Engine::projectConstraint(std::shared_ptr<Constraint> constraint)
{
    // Calculate $C$
    const double C = constraint->calculateValue();

    // Calculate $\Nabla C$
    const Eigen::VectorXd grad_C = constraint->calculateGrad();

    // Calculate $\mathbf{M}^{-1}$
    const unsigned int n = constraint->m_indices.size();
    std::vector<double> inverse_mass_raw(3 * n);
    for (unsigned int j = 0; j < n; ++ j)
    {
        inverse_mass_raw[j + 0] = 1.0 / m_particles[constraint->m_indices[j]].m;
        inverse_mass_raw[j + 1] = 1.0 / m_particles[constraint->m_indices[j]].m;
        inverse_mass_raw[j + 2] = 1.0 / m_particles[constraint->m_indices[j]].m;
    }
    const auto inverse_M_diagnal = Eigen::Map<Eigen::VectorXd>(inverse_mass_raw.data(), 3 * n);

    // Calculate $s$
    const double s = C / (grad_C.transpose() * inverse_M_diagnal.asDiagonal() * grad_C);

    // Calculate $\Delta x$
    const Eigen::VectorXd delta_x = - s * inverse_M_diagnal.asDiagonal() * grad_C;

    // Update predicted positions
    for (unsigned int j = 0; j < n; ++ j)
    {
        m_particles[constraint->m_indices[j]].p += delta_x;
    }
}

