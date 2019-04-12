#ifndef engine_hpp
#define engine_hpp

#include <memory>
#include <vector>
#include <elasty/particle.hpp>

namespace elasty
{
    class Constraint;

    class Engine
    {
    public:

        void stepTime();

        virtual void initializeScene() = 0;

        virtual void setExternalForces() = 0;
        virtual void generateCollisionConstraints() = 0;
        virtual void updateVelocities() = 0;

        void addConstraint(std::shared_ptr<Constraint> constraint)
        {
            m_constraints.push_back(constraint);
        }

        void addInstantConstraint(std::shared_ptr<Constraint> constraint)
        {
            m_instant_constraints.push_back(constraint);
        }

        std::vector<Particle> m_particles;

        std::vector<std::shared_ptr<Constraint>> m_constraints;
        std::vector<std::shared_ptr<Constraint>> m_instant_constraints;

    protected:

        void projectConstraint(std::shared_ptr<Constraint> constraint);

        double m_dt = 1.0 / 60.0;
        unsigned int m_num_iterations = 10;
    };
}

#endif /* engine_hpp */
