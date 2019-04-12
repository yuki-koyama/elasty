#ifndef ELASTY_HPP
#define ELASTY_HPP

#include <cassert>
#include <memory>
#include <vector>
#include <Eigen/Core>

namespace elasty
{
    class Constraint;

    struct Particle
    {
        Eigen::Vector3d x;
        Eigen::Vector3d v;
        Eigen::Vector3d p;
        double m;
        unsigned int i;
    };

    class Engine
    {
    public:

        virtual void initializeScene() = 0;
        virtual void stepTime() = 0;

        virtual void addConstraint(std::shared_ptr<Constraint> constraint)
        {
            m_constraints.push_back(constraint);
        }

        virtual void addInstantConstraint(std::shared_ptr<Constraint> constraint)
        {
            m_instant_constraints.push_back(constraint);
        }

        std::vector<Particle> m_particles;

        std::vector<std::shared_ptr<Constraint>> m_constraints;
        std::vector<std::shared_ptr<Constraint>> m_instant_constraints;

    protected:

        void projectConstraint(std::shared_ptr<Constraint> constraint);
    };

    class Constraint
    {
    public:
        Constraint(const Engine* engine, const std::vector<unsigned int>& indices) :
        m_engine(engine),
        m_indices(indices)
        {
        }

        virtual double calculateValue() = 0;
        virtual Eigen::VectorXd calculateGrad() = 0;

        const Engine* m_engine;
        const std::vector<unsigned int> m_indices;
    };

    class EnvironmentalCollisionConstraint final : public Constraint
    {
    public:

        EnvironmentalCollisionConstraint(const Engine* engine, const std::vector<unsigned int>& indices, const Eigen::Vector3d& n, const double d) :
        Constraint(engine, indices),
        m_n(n),
        m_d(d)
        {
            assert(m_indices.size() == 1);
        }

        double calculateValue() override
        {
            const Eigen::Vector3d& x = m_engine->m_particles[m_indices[0]].p;
            return m_n.transpose() * x - m_d;
        }

        Eigen::VectorXd calculateGrad() override
        {
            return m_n;
        }

    private:

        const Eigen::Vector3d m_n;
        const double m_d;
    };
}

#endif // ELASTY_HPP
