#ifndef ELASTY_ENGINE_HPP
#define ELASTY_ENGINE_HPP

#include <memory>
#include <vector>

namespace elasty
{
    class AbstractConstraint;
    struct Particle;

    class AbstractEngine
    {
    public:
        void stepTime();

        virtual void initializeScene() = 0;

        virtual void setExternalForces()            = 0;
        virtual void generateCollisionConstraints() = 0;
        virtual void updateVelocities()             = 0;

        void clearScene();

        std::vector<std::shared_ptr<Particle>> m_particles;

        std::vector<std::shared_ptr<AbstractConstraint>> m_constraints;
        std::vector<std::shared_ptr<AbstractConstraint>> m_instant_constraints;

        double getDeltaTime() const { return m_delta_time; }

        void setNumIters(unsigned int num_iters) { m_num_iterations = num_iters; }

    protected:
        void addConstraint(std::shared_ptr<AbstractConstraint> constraint);
        void addInstantConstraint(std::shared_ptr<AbstractConstraint> constraint);

    private:
        double       m_delta_time     = 1.0 / 60.0;
        unsigned int m_num_iterations = 10;
    };
} // namespace elasty

#endif // ELASTY_ENGINE_HPP
