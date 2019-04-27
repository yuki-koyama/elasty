#ifndef constraint_hpp
#define constraint_hpp

#include <elasty/particle.hpp>
#include <memory>
#include <vector>
#include <Eigen/Core>

namespace elasty
{
    enum class ConstraintType { Bilateral, Unilateral };

    class Constraint
    {
    public:

        Constraint(const std::vector<std::shared_ptr<Particle>>& particles,
                   const double stiffness) :
        m_stiffness(stiffness),
        m_particles(particles)
        {
        }

        /// \brief Calculates the constraint function value C(x).
        virtual double calculateValue() = 0;

        /// \brief Calculate the derivative of the constraint function
        /// grad C(x).
        /// \details As constraints can have different vector sizes, it will
        /// store the result to the passed raw buffer that should be allocated
        /// in the caller, rather than returning a dynamically allocated
        /// variable-length vector. This method does not check whether the
        /// buffer is adequately allocated, or not.
        virtual void calculateGrad(double* grad_C) = 0;

        /// \brief Manipulate the associated particles by projecting them to the
        /// constraint manifold.
        /// \details This method should be called by the core engine. As this
        /// method directly updates the predicted positions of the associated
        /// particles, it is intended to be used in a Gauss-Seidel-style solver.
        virtual void projectParticles() = 0;

        /// \brief Return the constraint type (i.e., either unilateral or
        /// bilateral).
        virtual ConstraintType getType() = 0;

        /// \brief Stiffness of this constraint, which should be in [0, 1].
        double m_stiffness;

    protected:

        /// \brief Associated particles.
        /// \details The number of particles is (in most cases) solely
        /// determined in each constraint. For example, a distance constraint
        /// need to have exactly two particles. Some special constraints
        /// (e.g., shape-matching constraint) could have a variable number of
        /// particles.
        std::vector<std::shared_ptr<Particle>> m_particles;
    };

    template <int Num>
    class FixedNumConstraint : public Constraint
    {
    public:

        FixedNumConstraint(const std::vector<std::shared_ptr<Particle>>& particles,
                           const double stiffness) :
        Constraint(particles, stiffness),
        m_inv_M(constructInverseMassMatrix(m_particles))
        {
        }

        virtual double calculateValue() = 0;
        virtual void calculateGrad(double* grad_C) = 0;

        void projectParticles() final
        {
            // Calculate the constraint function value
            const double C = calculateValue();

            // Skip if it is a unilateral constraint and is satisfied
            if (getType() == ConstraintType::Unilateral && C >= 0.0) { return; }

            // Calculate the derivative of the constraint function
            const Eigen::Matrix<double, Num * 3, 1> grad_C = calculateGrad();

            // Skip if the gradient is sufficiently small
            if (grad_C.isApprox(Eigen::Matrix<double, Num * 3, 1>::Zero())) { return; }

            // Calculate s
            const double s = C / (grad_C.transpose() * m_inv_M.asDiagonal() * grad_C);

            // Calculate \Delta x
            const Eigen::Matrix<double, Num * 3, 1> delta_x = - s * m_inv_M.asDiagonal() * grad_C;
            assert(!delta_x.hasNaN());

            // Update predicted positions
            for (unsigned int j = 0; j < Num; ++ j)
            {
                m_particles[j]->p += m_stiffness * delta_x.segment(3 * j, 3);
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

        static Eigen::Matrix<double, Num * 3, 1> constructInverseMassMatrix(const std::vector<std::shared_ptr<elasty::Particle>>& particles)
        {
            Eigen::Matrix<double, Num * 3, 1> inv_M;
            for (unsigned int j = 0; j < Num; ++ j)
            {
                inv_M(j * 3 + 0) = particles[j]->w;
                inv_M(j * 3 + 1) = particles[j]->w;
                inv_M(j * 3 + 2) = particles[j]->w;
            }
            return inv_M;
        }
    };

    class BendingConstraint final : public FixedNumConstraint<4>
    {
    public:

        BendingConstraint(const std::shared_ptr<Particle> p_0,
                          const std::shared_ptr<Particle> p_1,
                          const std::shared_ptr<Particle> p_2,
                          const std::shared_ptr<Particle> p_3,
                          const double stiffness,
                          const double dihedral_angle);

        double calculateValue() override;
        void calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:

        const double m_dihedral_angle;
    };

    class DistanceConstraint final : public FixedNumConstraint<2>
    {
    public:

        DistanceConstraint(const std::shared_ptr<Particle> p_0,
                           const std::shared_ptr<Particle> p_1,
                           const double stiffness,
                           const double d);
        
        double calculateValue() override;
        void calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:

        const double m_d;
    };

    class EnvironmentalCollisionConstraint final : public FixedNumConstraint<1>
    {
    public:

        EnvironmentalCollisionConstraint(const std::shared_ptr<Particle> p_0,
                                         const double stiffness,
                                         const Eigen::Vector3d& n,
                                         const double d);

        double calculateValue() override;
        void calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Unilateral; }

    private:

        const Eigen::Vector3d m_n;
        const double m_d;
    };

    class FixedPointConstraint final : public FixedNumConstraint<1>
    {
    public:

        FixedPointConstraint(const std::shared_ptr<Particle> p_0,
                             const double stiffness,
                             const Eigen::Vector3d& point);

        double calculateValue() override;
        void calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:

        const Eigen::Vector3d m_point;
    };

    class IsometricBendingConstraint final : public FixedNumConstraint<4>
    {
    public:

        IsometricBendingConstraint(const std::shared_ptr<Particle> p_0,
                                   const std::shared_ptr<Particle> p_1,
                                   const std::shared_ptr<Particle> p_2,
                                   const std::shared_ptr<Particle> p_3,
                                   const double stiffness);

        double calculateValue() override;
        void calculateGrad(double* grad_C) override;
        ConstraintType getType() override { return ConstraintType::Bilateral; }

    private:

        Eigen::Matrix4d m_Q;
    };
}

#endif /* constraint_hpp */
