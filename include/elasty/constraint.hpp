#ifndef constraint_hpp
#define constraint_hpp

#include <vector>
#include <Eigen/Core>

namespace elasty
{
    class Engine;

    class Constraint
    {
    public:

        Constraint(const Engine* engine, const std::vector<unsigned int>& indices, const double stiffness) :
        m_engine(engine),
        m_indices(indices),
        m_stiffness(stiffness)
        {
        }

        virtual double calculateValue() = 0;
        virtual Eigen::VectorXd calculateGrad() = 0;

        const Engine* m_engine;
        const std::vector<unsigned int> m_indices;
        const double m_stiffness;
    };

    class DistanceConstraint final : public Constraint
    {
    public:

        DistanceConstraint(const Engine* engine, const std::vector<unsigned int>& indices, const double stiffness, const double d) :
        Constraint(engine, indices, stiffness),
        m_d(d)
        {
            assert(indices.size() == 2);
            assert(d >= 0.0);
        }

        double calculateValue() override;
        Eigen::VectorXd calculateGrad() override;

    private:

        const double m_d;
    };

    class EnvironmentalCollisionConstraint final : public Constraint
    {
    public:

        EnvironmentalCollisionConstraint(const Engine* engine,
                                         const std::vector<unsigned int>& indices,
                                         const double stiffness,
                                         const Eigen::Vector3d& n,
                                         const double d) :
        Constraint(engine, indices, stiffness),
        m_n(n),
        m_d(d)
        {
            assert(m_indices.size() == 1);
        }

        double calculateValue() override;
        Eigen::VectorXd calculateGrad() override;

    private:

        const Eigen::Vector3d m_n;
        const double m_d;
    };

    class FixedPointConstraint final : public Constraint
    {
    public:

        FixedPointConstraint(const Engine* engine,
                             const std::vector<unsigned int>& indices,
                             const double stiffness,
                             const Eigen::Vector3d& point);

        double calculateValue() override;
        Eigen::VectorXd calculateGrad() override;

    private:

        const Eigen::Vector3d m_point;
    };
}

#endif /* constraint_hpp */
