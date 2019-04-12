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
}

#endif /* engine_hpp */
