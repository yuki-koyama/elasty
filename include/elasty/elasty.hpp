#ifndef ELASTY_HPP
#define ELASTY_HPP

#include <memory>
#include <vector>
#include <Eigen/Core>

namespace elasty
{
    struct Vertex
    {
        Eigen::Vector3d x;
        Eigen::Vector3d v;
    };

    struct Constraint
    {
        virtual void calculateValue(double& C) = 0;
        virtual void calculateGrad(Eigen::Matrix3Xd& grad_C) = 0;
    };

    class Engine
    {
    public:

        virtual void initializeScene() = 0;
        virtual void stepTime() = 0;
        virtual void addConstraint(std::shared_ptr<Constraint> constraint) = 0;

    private:

        std::vector<std::shared_ptr<Constraint>> m_constraints;
    };
}

#endif // ELASTY_HPP
