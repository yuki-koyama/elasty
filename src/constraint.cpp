#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <algorithm>
#include <cstring>
#include <elasty/constraint.hpp>
#include <elasty/fem.hpp>
#include <elasty/particle.hpp>

namespace
{
    inline Eigen::Matrix3d convertVecToCrossOp(const Eigen::Vector3d& vec)
    {
        Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();

        mat(0, 1) = -vec(2);
        mat(0, 2) = +vec(1);
        mat(1, 0) = +vec(2);
        mat(1, 2) = -vec(0);
        mat(2, 0) = -vec(1);
        mat(2, 1) = +vec(0);

        return mat;
    };

    inline double calculateCotTheta(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
    {
        const double scaled_cos_theta = x.dot(y);
        const double scaled_sin_theta = x.cross(y).norm();
        return scaled_cos_theta / scaled_sin_theta;
    }
} // namespace

elasty::BendingConstraint::BendingConstraint(const std::shared_ptr<Particle> p_0,
                                             const std::shared_ptr<Particle> p_1,
                                             const std::shared_ptr<Particle> p_2,
                                             const std::shared_ptr<Particle> p_3,
                                             const double                    stiffness,
                                             const double                    compliance,
                                             const double                    delta_time,
                                             const double                    dihedral_angle)
    : FixedNumAbstractConstraint(std::vector<std::shared_ptr<Particle>>{p_0, p_1, p_2, p_3},
                                 stiffness,
                                 compliance,
                                 delta_time),
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

    const double current_dihedral_angle = std::acos(std::clamp(n_0.dot(n_1), -1.0, 1.0));

    assert(n_0.norm() > 0.0);
    assert(n_1.norm() > 0.0);
    assert(!std::isnan(current_dihedral_angle));

    return current_dihedral_angle - m_dihedral_angle;
}

// See the appendix of the original paper by Muller et al. (2007) for details.
void elasty::BendingConstraint::calculateGrad(double* grad_C)
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

    // If the current dihedral angle is sufficiently small or large (i.e., zero or pi), return zeros.
    // This is only an ad-hoc solution for stability and it needs to be solved in a more theoretically grounded way.
    constexpr double epsilon = 1e-12;
    if (1.0 - d * d < epsilon)
    {
        std::fill(grad_C, grad_C + 12, 0.0);
        return;
    }

    const double common_coeff = -1.0 / std::sqrt(1.0 - d * d);

    auto calc_grad_of_normalized_cross_prod_wrt_p_a =
        [](const Eigen::Vector3d& p_a, const Eigen::Vector3d& p_b, const Eigen::Vector3d& n) -> Eigen::Matrix3d {
        return +(1.0 / p_a.cross(p_b).norm()) * (-convertVecToCrossOp(p_b) + n * (n.cross(p_b)).transpose());
    };

    auto calc_grad_of_normalized_cross_prod_wrt_p_b =
        [](const Eigen::Vector3d& p_a, const Eigen::Vector3d& p_b, const Eigen::Vector3d& n) -> Eigen::Matrix3d {
        return -(1.0 / p_a.cross(p_b).norm()) * (-convertVecToCrossOp(p_a) + n * (n.cross(p_a)).transpose());
    };

    const Eigen::Matrix3d partial_n_0_per_partial_p_1 = calc_grad_of_normalized_cross_prod_wrt_p_a(p_1, p_2, n_0);
    const Eigen::Matrix3d partial_n_1_per_partial_p_1 = calc_grad_of_normalized_cross_prod_wrt_p_a(p_1, p_3, n_1);
    const Eigen::Matrix3d partial_n_0_per_partial_p_2 = calc_grad_of_normalized_cross_prod_wrt_p_b(p_1, p_2, n_0);
    const Eigen::Matrix3d partial_n_1_per_partial_p_3 = calc_grad_of_normalized_cross_prod_wrt_p_b(p_1, p_3, n_1);

    const Eigen::Vector3d grad_C_wrt_p_1 =
        common_coeff * (partial_n_0_per_partial_p_1.transpose() * n_1 + partial_n_1_per_partial_p_1.transpose() * n_0);
    const Eigen::Vector3d grad_C_wrt_p_2 = common_coeff * partial_n_0_per_partial_p_2.transpose() * n_1;
    const Eigen::Vector3d grad_C_wrt_p_3 = common_coeff * partial_n_1_per_partial_p_3.transpose() * n_0;
    const Eigen::Vector3d grad_C_wrt_p_0 = -grad_C_wrt_p_1 - grad_C_wrt_p_2 - grad_C_wrt_p_3;

    std::memcpy(grad_C + (3 * 0), grad_C_wrt_p_0.data(), sizeof(double) * 3);
    std::memcpy(grad_C + (3 * 1), grad_C_wrt_p_1.data(), sizeof(double) * 3);
    std::memcpy(grad_C + (3 * 2), grad_C_wrt_p_2.data(), sizeof(double) * 3);
    std::memcpy(grad_C + (3 * 3), grad_C_wrt_p_3.data(), sizeof(double) * 3);
}

elasty::ContinuumTriangleConstraint::ContinuumTriangleConstraint(const std::shared_ptr<Particle> p_0,
                                                                 const std::shared_ptr<Particle> p_1,
                                                                 const std::shared_ptr<Particle> p_2,
                                                                 const double                    stiffness,
                                                                 const double                    compliance,
                                                                 const double                    delta_time,
                                                                 const double                    youngs_modulus,
                                                                 const double                    poisson_ratio)
    : FixedNumAbstractConstraint(std::vector<std::shared_ptr<Particle>>{p_0, p_1, p_2},
                                 stiffness,
                                 compliance,
                                 delta_time),
      m_first_lame(fem::calcFirstLame(youngs_modulus, poisson_ratio)),
      m_second_lame(fem::calcSecondLame(youngs_modulus, poisson_ratio))
{
    const Eigen::Vector3d& x_0 = m_particles[0]->x;
    const Eigen::Vector3d& x_1 = m_particles[1]->x;
    const Eigen::Vector3d& x_2 = m_particles[2]->x;

    // Calculate the two axes for defining material coordinates
    const Eigen::Vector3d r_1    = x_1 - x_0;
    const Eigen::Vector3d r_2    = x_2 - x_0;
    const Eigen::Vector3d cross  = r_1.cross(r_2);
    const Eigen::Vector3d axis_1 = r_1.normalized();
    const Eigen::Vector3d axis_2 = cross.cross(axis_1).normalized();

    // Calculate the rest positions in the material coordinates
    const Eigen::Vector2d mat_x_0(axis_1.dot(x_0), axis_2.dot(x_0));
    const Eigen::Vector2d mat_x_1(axis_1.dot(x_1), axis_2.dot(x_1));
    const Eigen::Vector2d mat_x_2(axis_1.dot(x_2), axis_2.dot(x_2));

    // Calculate the rest shape matrix
    Eigen::Matrix2d rest_D;
    rest_D.col(0) = mat_x_1 - mat_x_0;
    rest_D.col(1) = mat_x_2 - mat_x_0;

    // Calculate the inverse of the rest shape matrix
    assert(rest_D.determinant() > 0.0);
    m_rest_D_inv = rest_D.inverse();

    // Calculate the area of the rest configuration
    m_rest_area = 0.5 * cross.norm();
}

double elasty::ContinuumTriangleConstraint::calculateValue()
{
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;
    const Eigen::Vector3d& x_2 = m_particles[2]->p;

    // Calculate the shape matrix
    Eigen::Matrix<double, 3, 2> D;
    D.col(0) = x_1 - x_0;
    D.col(1) = x_2 - x_0;

    // Calculate the deformation gradient
    const Eigen::Matrix<double, 3, 2> F = D * m_rest_D_inv;

    // Calculate the Green strain tensor
    const Eigen::Matrix2d epsilon = fem::calcGreenStrain(F);

    // Calculate the strain tensor based on the Saint Venant–Kirchhoff model
    const Eigen::Matrix2d S =
        m_first_lame * epsilon.trace() * Eigen::Matrix2d::Identity() + 2.0 * m_second_lame * epsilon;

    // Calculate the strain energy density
    const double psi = 0.5 * (epsilon.transpose() * S).trace();

    // Return the constraint value
    return m_rest_area * psi;
}

void elasty::ContinuumTriangleConstraint::calculateGrad(double* grad_C)
{
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;
    const Eigen::Vector3d& x_2 = m_particles[2]->p;

    // Calculate the shape matrix
    Eigen::Matrix<double, 3, 2> D;
    D.col(0) = x_1 - x_0;
    D.col(1) = x_2 - x_0;

    // Calculate the deformation gradient
    const Eigen::Matrix<double, 3, 2> F = D * m_rest_D_inv;

    // Calculate the Green strain tensor
    const Eigen::Matrix2d epsilon = fem::calcGreenStrain(F);

    // Calculate the strain tensor based on the Saint Venant–Kirchhoff model
    const Eigen::Matrix2d S =
        m_first_lame * epsilon.trace() * Eigen::Matrix2d::Identity() + 2.0 * m_second_lame * epsilon;

    // Calculate the first Piola-Kirchhoff stress tensor
    const Eigen::Matrix<double, 3, 2> P = F * S;

    // Calculate the gradient of the constraint
    const Eigen::Matrix<double, 3, 2> grad_12 = m_rest_area * P * m_rest_D_inv.transpose();
    const Eigen::Vector3d             grad_0  = -grad_12.col(0) - grad_12.col(1);

    // Copy the results
    std::memcpy(grad_C + 0, grad_0.data(), sizeof(double) * 3);
    std::memcpy(grad_C + 3, grad_12.data(), sizeof(double) * 6);
}

elasty::DistanceConstraint::DistanceConstraint(const std::shared_ptr<Particle> p_0,
                                               const std::shared_ptr<Particle> p_1,
                                               const double                    stiffness,
                                               const double                    compliance,
                                               const double                    delta_time,
                                               const double                    d)
    : FixedNumAbstractConstraint(std::vector<std::shared_ptr<Particle>>{p_0, p_1}, stiffness, compliance, delta_time),
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

void elasty::DistanceConstraint::calculateGrad(double* grad_C)
{
    const Eigen::Vector3d& x_0 = m_particles[0]->p;
    const Eigen::Vector3d& x_1 = m_particles[1]->p;

    const Eigen::Vector3d r = x_0 - x_1;

    const double dist = r.norm();

    constexpr double epsilon = 1e-24;

    // Calculate a normalized vector, where a random direction is selected when the points are degenerated
    const Eigen::Vector3d n = (dist < epsilon) ? Eigen::Vector3d::Random().normalized() : (1.0 / dist) * r;

    grad_C[0] = +n(0);
    grad_C[1] = +n(1);
    grad_C[2] = +n(2);
    grad_C[3] = -n(0);
    grad_C[4] = -n(1);
    grad_C[5] = -n(2);
}

elasty::EnvironmentalCollisionConstraint::EnvironmentalCollisionConstraint(const std::shared_ptr<Particle> p_0,
                                                                           const double                    stiffness,
                                                                           const double                    compliance,
                                                                           const double                    delta_time,
                                                                           const Eigen::Vector3d&          n,
                                                                           const double                    d)
    : FixedNumAbstractConstraint(std::vector<std::shared_ptr<Particle>>{p_0}, stiffness, compliance, delta_time),
      m_n(n),
      m_d(d)
{
}

double elasty::EnvironmentalCollisionConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_particles[0]->p;
    return m_n.transpose() * x - m_d;
}

void elasty::EnvironmentalCollisionConstraint::calculateGrad(double* grad_C)
{
    std::memcpy(grad_C, m_n.data(), sizeof(double) * 3);
}

elasty::FixedPointConstraint::FixedPointConstraint(const std::shared_ptr<Particle> p_0,
                                                   const double                    stiffness,
                                                   const double                    compliance,
                                                   const double                    delta_time,
                                                   const Eigen::Vector3d&          point)
    : FixedNumAbstractConstraint(std::vector<std::shared_ptr<Particle>>{p_0}, stiffness, compliance, delta_time),
      m_point(point)
{
}

double elasty::FixedPointConstraint::calculateValue()
{
    const Eigen::Vector3d& x = m_particles[0]->p;
    return (x - m_point).norm();
}

void elasty::FixedPointConstraint::calculateGrad(double* grad_C)
{
    const Eigen::Vector3d& x    = m_particles[0]->p;
    const Eigen::Vector3d  r    = x - m_point;
    const double           dist = r.norm();

    constexpr double epsilon = 1e-24;

    // Calculate a normalized vector, where a random direction is selected when the points are degenerated
    const Eigen::Vector3d n = (dist < epsilon) ? Eigen::Vector3d::Random().normalized() : (1.0 / dist) * r;

    std::memcpy(grad_C, n.data(), sizeof(double) * 3);
}

elasty::IsometricBendingConstraint::IsometricBendingConstraint(const std::shared_ptr<Particle> p_0,
                                                               const std::shared_ptr<Particle> p_1,
                                                               const std::shared_ptr<Particle> p_2,
                                                               const std::shared_ptr<Particle> p_3,
                                                               const double                    stiffness,
                                                               const double                    compliance,
                                                               const double                    delta_time)
    : FixedNumAbstractConstraint(std::vector<std::shared_ptr<Particle>>{p_0, p_1, p_2, p_3},
                                 stiffness,
                                 compliance,
                                 delta_time)
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

    const double cot_01 = calculateCotTheta(e0, -e1);
    const double cot_02 = calculateCotTheta(e0, -e2);
    const double cot_03 = calculateCotTheta(e0, e3);
    const double cot_04 = calculateCotTheta(e0, e4);

    const Eigen::Vector4d K = Eigen::Vector4d(cot_01 + cot_04, cot_02 + cot_03, -cot_01 - cot_02, -cot_03 - cot_04);

    const double A_0 = 0.5 * e0.cross(e1).norm();
    const double A_1 = 0.5 * e0.cross(e3).norm();

    m_Q = (3.0 / (A_0 + A_1)) * K * K.transpose();
}

double elasty::IsometricBendingConstraint::calculateValue()
{
    double sum = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
        for (unsigned int j = 0; j < 4; ++j)
        {
            sum += m_Q(i, j) * double(m_particles[i]->p.transpose() * m_particles[j]->p);
        }
    }
    return 0.5 * sum;
}

void elasty::IsometricBendingConstraint::calculateGrad(double* grad_C)
{
    for (unsigned int i = 0; i < 4; ++i)
    {
        Eigen::Vector3d sum = Eigen::Vector3d::Zero();
        for (unsigned int j = 0; j < 4; ++j)
        {
            sum += m_Q(i, j) * m_particles[j]->p;
        }
        std::memcpy(grad_C + (3 * i), sum.data(), sizeof(double) * 3);
    }
}

elasty::ShapeMatchingConstraint::ShapeMatchingConstraint(const std::vector<std::shared_ptr<Particle>>& particles,
                                                         const double                                  stiffness,
                                                         const double                                  compliance,
                                                         const double                                  delta_time)
    : VariableNumConstraint(particles, stiffness, compliance, delta_time)
{
    // Calculate the initial center of mass and the total mass
    Eigen::Vector3d x_0_cm = Eigen::Vector3d::Zero();
    m_total_mass           = 0.0;
    for (int i = 0; i < m_particles.size(); ++i)
    {
        m_total_mass += m_particles[i]->m;
        x_0_cm += m_particles[i]->m * m_particles[i]->p;
    }
    x_0_cm /= m_total_mass;

    // Calculate q
    m_q.resize(m_particles.size());
    for (int i = 0; i < m_particles.size(); ++i)
    {
        m_q[i] = m_particles[i]->x - x_0_cm;
    }
}

double elasty::ShapeMatchingConstraint::calculateValue()
{
    throw std::runtime_error("ShapeMatchingConstraint does not directly provide its cost value or the gradient.");
}

void elasty::ShapeMatchingConstraint::calculateGrad(double* grad_C)
{
    throw std::runtime_error("ShapeMatchingConstraint does not directly provide its cost value or the gradient.");
}

void elasty::ShapeMatchingConstraint::projectParticles(const AlgorithmType type)
{
    // Calculate the current center of mass
    Eigen::Vector3d x_cm = Eigen::Vector3d::Zero();
    for (int i = 0; i < m_particles.size(); ++i)
    {
        x_cm += m_particles[i]->m * m_particles[i]->p;
    }
    x_cm /= m_total_mass;

    // Calculate A_pq
    Eigen::Matrix3d A_pq = Eigen::Matrix3d::Zero();
    for (int i = 0; i < m_particles.size(); ++i)
    {
        A_pq += m_particles[i]->m * (m_particles[i]->p - x_cm) * m_q[i].transpose();
    }

    // Calculate the rotation matrix
    const auto            ATA          = A_pq.transpose() * A_pq;
    const auto            eigen_solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(ATA);
    const auto            S_inv        = eigen_solver.operatorInverseSqrt();
    const Eigen::Matrix3d R            = A_pq * S_inv;

    assert(R.determinant() > 0);

    // Update the particle positions
    for (int i = 0; i < m_particles.size(); ++i)
    {
        // Calculate the goal position
        const Eigen::Vector3d g = R * m_q[i] + x_cm;

        // Move the particle
        m_particles[i]->p = m_stiffness * g + (1.0 - m_stiffness) * m_particles[i]->p;
    }
}
