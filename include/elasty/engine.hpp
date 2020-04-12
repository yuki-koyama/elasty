#ifndef ELASTY_ENGINE_HPP
#define ELASTY_ENGINE_HPP

#include <elasty/algorithm-type.hpp>
#include <memory>
#include <vector>

namespace elasty
{
    class AbstractConstraint;
    struct Particle;

    class AbstractEngine
    {
    public:
        AbstractEngine(const double        delta_time           = 1.0 / 60.0,
                       const unsigned int  num_constraint_iters = 10,
                       const unsigned int  num_substeps         = 1,
                       const AlgorithmType algorithm_type       = AlgorithmType::Xpbd);

        void stepTime();

        virtual void initializeScene() = 0;

        virtual void setExternalForces()            = 0;
        virtual void generateCollisionConstraints() = 0;
        virtual void updateVelocities()             = 0;

        void clearScene();

        /// \brief Getter of the delta time value
        double getDeltaTime() const { return m_delta_time; }

        /// \brief Getter of the particles this engine manages
        const std::vector<std::shared_ptr<Particle>>& getParticles() const { return m_particles; }

    protected:
        void addConstraint(std::shared_ptr<AbstractConstraint> constraint);
        void addInstantConstraint(std::shared_ptr<AbstractConstraint> constraint);

        std::vector<std::shared_ptr<Particle>> m_particles;

        std::vector<std::shared_ptr<AbstractConstraint>> m_constraints;
        std::vector<std::shared_ptr<AbstractConstraint>> m_instant_constraints;

    private:
        const double        m_delta_time;
        const unsigned int  m_num_constraint_iters;
        const unsigned int  m_num_substeps;
        const AlgorithmType m_algorithm_type;
    };
} // namespace elasty

#endif // ELASTY_ENGINE_HPP
