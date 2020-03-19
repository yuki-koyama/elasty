#ifndef constraint_hpp
#define constraint_hpp

#include <Eigen/Core>
#include <elasty/particle.hpp>
#include <memory>
#include <vector>

namespace elasty
{
    enum class ConstraintType
    {
        Bilateral,
        Unilateral
    };

    enum class AlgorithmType
    {
        Pbd,
        Xpbd
    };

    class AbstractConstraint
    {
    public:
        AbstractConstraint(const std::vector<std::shared_ptr<Particle>>& particles,
                           const double                                  stiffness,  ///< for PBD
                           const double                                  compliance, ///< for XPBD
                           const double                                  dt          ///< for XPBD
                           )
            : m_stiffness(stiffness),
              m_lagrange_multiplier(0.0),
              m_compliance(compliance),
              m_dt(dt),
              m_particles(particles)
        {
        }

        /// \brief Calculates the constraint function value C(x).
        virtual double calculateValue() = 0;

        /// \brief Calculate the derivative of the constraint function grad C(x).
        ///
        /// \details As constraints can have different vector sizes, it will store the result to the passed raw buffer
        /// that should be allocated in the caller, rather than returning a dynamically allocated variable-length
        /// vector. This method does not check whether the buffer is adequately allocated, or not.
        virtual void calculateGrad(double* grad_C) = 0;

        /// \brief Manipulate the associated particles by projecting them to the constraint manifold.
        ///
        /// \details This method should be called by the core engine. As this method directly updates the predicted
        /// positions of the associated particles, it is intended to be used in a Gauss-Seidel-style solver.
        virtual void projectParticles(const AlgorithmType type = AlgorithmType::Pbd) = 0;

        /// \brief Return the constraint type (i.e., either unilateral or bilateral).
        virtual ConstraintType getType() = 0;

        /// \brief Stiffness of this constraint, which should be in [0, 1].
        ///
        /// \details This value will not be used in XPBD; instead, the compliance value will be used.
        double m_stiffness;

        /// \brief Lagrange multipler for this constraint, used in XPBD only.
        double m_lagrange_multiplier;

        /// \brief Compliance (inverse stiffness), used in XPBD only.
        ///
        /// \details Should be set to zero if the constraint is hard (e.g., collision constraints).
        double m_compliance;

        /// \brief Time step, used in XPBD only.
        double m_dt;

    protected:
        /// \brief Associated particles.
        ///
        /// \details The number of particles is (in most cases) solely determined in each constraint. For example, a
        /// distance constraint need to have exactly two particles. Some special constraints (e.g., shape-matching
        /// constraint) could have a variable number of particles.
        std::vector<std::shared_ptr<Particle>> m_particles;
    };

    template <int Num> class FixedNumAbstractConstraint : public AbstractConstraint
    {
    public:
        FixedNumAbstractConstraint(const std::vector<std::shared_ptr<Particle>>& particles,
                                   const double                                  stiffness,  ///< for PBD
                                   const double                                  compliance, ///< for XPBD
                                   const double                                  dt          ///< for XPBD
                                   )
            : AbstractConstraint(particles, stiffness, compliance, dt), m_inv_M(constructInverseMassMatrix(m_particles))
        {
        }

        virtual double calculateValue()              = 0;
        virtual void   calculateGrad(double* grad_C) = 0;

        void projectParticles(const AlgorithmType type = AlgorithmType::Pbd) final
        {
            // Calculate the constraint function value
            const double C = calculateValue();

            // Skip if it is a unilateral constraint and is satisfied
            if (getType() == ConstraintType::Unilateral && C >= 0.0)
            {
                return;
            }

            // Calculate the derivative of the constraint function
            const Eigen::Matrix<double, Num * 3, 1> grad_C = calculateGrad();

            // Skip if the gradient is sufficiently small
#if 1
            constexpr double very_small_value = 1e-10;
            if (grad_C.norm() < very_small_value)
#else
            if (grad_C.isApprox(Eigen::Matrix<double, Num * 3, 1>::Zero()))
#endif
            {
                return;
            }

            switch (type)
            {
                case AlgorithmType::Pbd:
                {
                    // Calculate s
                    const double s = -C / (grad_C.transpose() * m_inv_M.asDiagonal() * grad_C);

                    // Calculate \Delta x
                    const Eigen::Matrix<double, Num * 3, 1> delta_x = m_stiffness * s * m_inv_M.asDiagonal() * grad_C;
                    assert(!delta_x.hasNaN());

                    // Update predicted positions
                    for (unsigned int j = 0; j < Num; ++j)
                    {
                        m_particles[j]->p += delta_x.segment(3 * j, 3);
                    }

                    break;
                }
                case AlgorithmType::Xpbd:
                {
                    // Calculate time-scaled compliance
                    const double alpha_tilde = m_compliance / (m_dt * m_dt);

                    // Calculate \Delta lagrange multiplier
                    const double delta_lagrange_multiplier =
                        (-C - alpha_tilde * m_lagrange_multiplier) /
                        (grad_C.transpose() * m_inv_M.asDiagonal() * grad_C + alpha_tilde);

                    // Calculate \Delta x
                    const Eigen::Matrix<double, Num * 3, 1> delta_x =
                        delta_lagrange_multiplier * m_inv_M.asDiagonal() * grad_C;
                    assert(!delta_x.hasNaN());

                    // Update predicted positions
                    for (unsigned int j = 0; j < Num; ++j)
                    {
                        m_particles[j]->p += delta_x.segment(3 * j, 3);
                    }

                    // Update the lagrange multiplier
                    m_lagrange_multiplier += delta_lagrange_multiplier;

                    break;
                }
            }
        }

    private:
        const Eigen::Matrix<double, Num * 3, 1> m_inv_M;

        Eigen::Matrix<double, Num * 3, 1> calculateGrad()
        {
            Eigen::Matrix<double, Num * 3, 1> grad_C;
            calculateGrad(grad_C.data());
            return grad_C;
        }

        static Eigen::Matrix<double, Num * 3, 1>
        constructInverseMassMatrix(const std::vector<std::shared_ptr<elasty::Particle>>& particles)
        {
            Eigen::Matrix<double, Num * 3, 1> inv_M;
            for (unsigned int j = 0; j < Num; ++j)
            {
                inv_M(j * 3 + 0) = particles[j]->w;
                inv_M(j * 3 + 1) = particles[j]->w;
                inv_M(j * 3 + 2) = particles[j]->w;
            }
            return inv_M;
        }
    };

    class VariableNumConstraint : public AbstractConstraint
    {
    public:
        VariableNumConstraint(const std::vector<std::shared_ptr<Particle>>& particles,
                              const double                                  stiffness,  ///< for PBD
                              const double                                  compliance, ///< for XPBD
                              const double                                  dt          ///< for XPBD
                              )
            : AbstractConstraint(particles, stiffness, compliance, dt), m_inv_M(constructInverseMassMatrix(m_particles))
        {
        }

        virtual double calculateValue()              = 0;
        virtual void   calculateGrad(double* grad_C) = 0;

        virtual void projectParticles(const AlgorithmType type = AlgorithmType::Pbd)
        {
            // Calculate the constraint function value
            const double C = calculateValue();

            // Skip if it is a unilateral constraint and is satisfied
            if (getType() == ConstraintType::Unilateral && C >= 0.0)
            {
                return;
            }

            // Calculate the derivative of the constraint function
            const Eigen::VectorXd grad_C = calculateGrad();

            // Skip if the gradient is sufficiently small
#if 1
            constexpr double very_small_value = 1e-10;
            if (grad_C.norm() < very_small_value)
#else
            if (grad_C.isApprox(Eigen::VectorXd::Zero(grad_C.size())))
#endif
            {
                return;
            }

            switch (type)
            {
                case AlgorithmType::Pbd:
                {
                    // Calculate s
                    const double s = -C / (grad_C.transpose() * m_inv_M.asDiagonal() * grad_C);

                    // Calculate \Delta x
                    const Eigen::VectorXd delta_x = m_stiffness * s * m_inv_M.asDiagonal() * grad_C;
                    assert(!delta_x.hasNaN());

                    // Update predicted positions
                    for (unsigned int j = 0; j < m_particles.size(); ++j)
                    {
                        m_particles[j]->p += delta_x.segment(3 * j, 3);
                    }

                    break;
                }
                case AlgorithmType::Xpbd:
                {
                    // Calculate time-scaled compliance
                    const double alpha_tilde = m_compliance / (m_dt * m_dt);

                    // Calculate \Delta lagrange multiplier
                    const double delta_lagrange_multiplier =
                        (-C - alpha_tilde * m_lagrange_multiplier) /
                        (grad_C.transpose() * m_inv_M.asDiagonal() * grad_C + alpha_tilde);

                    // Calculate \Delta x
                    const Eigen::VectorXd delta_x = delta_lagrange_multiplier * m_inv_M.asDiagonal() * grad_C;
                    assert(!delta_x.hasNaN());

                    // Update predicted positions
                    for (unsigned int j = 0; j < m_particles.size(); ++j)
                    {
                        m_particles[j]->p += delta_x.segment(3 * j, 3);
                    }

                    // Update the lagrange multiplier
                    m_lagrange_multiplier += delta_lagrange_multiplier;

                    break;
                }
            }
        }

    private:
        const Eigen::VectorXd m_inv_M;

        Eigen::VectorXd calculateGrad()
        {
            Eigen::VectorXd grad_C(m_particles.size() * 3);
            calculateGrad(grad_C.data());
            return grad_C;
        }

        static Eigen::VectorXd
        constructInverseMassMatrix(const std::vector<std::shared_ptr<elasty::Particle>>& particles)
        {
            Eigen::VectorXd inv_M(particles.size() * 3);
            for (unsigned int j = 0; j < particles.size(); ++j)
            {
                inv_M(j * 3 + 0) = particles[j]->w;
                inv_M(j * 3 + 1) = particles[j]->w;
                inv_M(j * 3 + 2) = particles[j]->w;
            }
            return inv_M;
        }
    };

    class BendingConstraint final : public FixedNumAbstractConstraint<4>
    {
    public:
        BendingConstraint(const std::shared_ptr<Particle> p_0,
                          const std::shared_ptr<Particle> p_1,
                          const std::shared_ptr<Particle> p_2,
                          const std::shared_ptr<Particle> p_3,
                          const double                    stiffness,  ///< for PBD
                          const double                    compliance, ///< for XPBD
                          const double                    dt,         ///< for XPBD
                          const double                    dihedral_angle);

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:
        const double m_dihedral_angle;
    };

    class ContinuumTriangleConstraint final : public FixedNumAbstractConstraint<3>
    {
    public:
        ContinuumTriangleConstraint(const std::shared_ptr<Particle> p_0,
                                    const std::shared_ptr<Particle> p_1,
                                    const std::shared_ptr<Particle> p_2,
                                    const double                    stiffness,  ///< for PBD
                                    const double                    compliance, ///< for XPBD
                                    const double                    dt,         ///< for XPBD
                                    const double                    youngs_modulus,
                                    const double                    poisson_ratio);

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:
        const double m_first_lame;
        const double m_second_lame;

        double          m_rest_area;
        Eigen::Matrix2d m_rest_D_inv;
    };

    class DistanceConstraint final : public FixedNumAbstractConstraint<2>
    {
    public:
        DistanceConstraint(const std::shared_ptr<Particle> p_0,
                           const std::shared_ptr<Particle> p_1,
                           const double                    stiffness,  ///< for PBD
                           const double                    compliance, ///< for XPBD
                           const double                    dt,         ///< for XPBD
                           const double                    d);

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:
        const double m_d;
    };

    class EnvironmentalCollisionConstraint final : public FixedNumAbstractConstraint<1>
    {
    public:
        EnvironmentalCollisionConstraint(const std::shared_ptr<Particle> p_0,
                                         const double                    stiffness,  ///< for PBD
                                         const double                    compliance, ///< for XPBD
                                         const double                    dt,         ///< for XPBD
                                         const Eigen::Vector3d&          n,
                                         const double                    d);

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Unilateral; }

    private:
        const Eigen::Vector3d m_n;
        const double          m_d;
    };

    class FixedPointConstraint final : public FixedNumAbstractConstraint<1>
    {
    public:
        FixedPointConstraint(const std::shared_ptr<Particle> p_0,
                             const double                    stiffness,  ///< for PBD
                             const double                    compliance, ///< for XPBD
                             const double                    dt,         ///< for XPBD
                             const Eigen::Vector3d&          point);

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:
        const Eigen::Vector3d m_point;
    };

    class IsometricBendingConstraint final : public FixedNumAbstractConstraint<4>
    {
    public:
        IsometricBendingConstraint(const std::shared_ptr<Particle> p_0,
                                   const std::shared_ptr<Particle> p_1,
                                   const std::shared_ptr<Particle> p_2,
                                   const std::shared_ptr<Particle> p_3,
                                   const double                    stiffness,  ///< for PBD
                                   const double                    compliance, ///< for XPBD
                                   const double                    dt          ///< for XPBD
        );

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:
        Eigen::Matrix4d m_Q;
    };

    class ShapeMatchingConstraint final : public VariableNumConstraint
    {
    public:
        ShapeMatchingConstraint(const std::vector<std::shared_ptr<Particle>>& particles,
                                const double                                  stiffness,  ///< for PBD
                                const double                                  compliance, ///< for XPBD
                                const double                                  dt          ///< for XPBD
        );

        double         calculateValue() override;
        void           calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

        /// \details Unlike other PBD constraints, ShapeMatchingConstraint does not offer a closed-form constraint
        /// function, but does offer a particle projection procedure directly. Thus, this class needs to override the
        /// particle projection method.
        void projectParticles(const AlgorithmType type = AlgorithmType::Pbd) override;

    private:
        double                       m_total_mass;
        std::vector<Eigen::Vector3d> m_q;
    };
} // namespace elasty

#endif /* constraint_hpp */
