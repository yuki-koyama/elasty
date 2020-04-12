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
        AbstractEngine(const double        delta_frame_time     = 1.0 / 60.0,
                       const unsigned int  num_constraint_iters = 1,
                       const unsigned int  num_substeps         = 8,
                       const AlgorithmType algorithm_type       = AlgorithmType::Xpbd);

        void proceedFrame();

        virtual void initializeScene() = 0;

        virtual void setExternalForces()            = 0;
        virtual void generateCollisionConstraints() = 0;
        virtual void updateVelocities()             = 0;

        void clearScene();

        /// \brief Getter of the delta physics time
        ///
        /// \details The value equals to the delta frame time devided by the number of substeps
        double getDeltaPhysicsTime() const { return m_delta_physics_time; }

        /// \brief Getter of the delta frame time
        ///
        /// \details The value equals to the delta physics time multipled by the number of substeps
        double getDeltaFrameTime() const { return m_delta_frame_time; }

        /// \brief Getter of the particles this engine manages
        const std::vector<std::shared_ptr<Particle>>& getParticles() const { return m_particles; }

    protected:
        void addConstraint(std::shared_ptr<AbstractConstraint> constraint);
        void addInstantConstraint(std::shared_ptr<AbstractConstraint> constraint);

        std::vector<std::shared_ptr<Particle>> m_particles;

        std::vector<std::shared_ptr<AbstractConstraint>> m_constraints;
        std::vector<std::shared_ptr<AbstractConstraint>> m_instant_constraints;

    private:
        const double        m_delta_physics_time;
        const double        m_delta_frame_time;
        const unsigned int  m_num_constraint_iters;
        const unsigned int  m_num_substeps;
        const AlgorithmType m_algorithm_type;

        void proceedSubstep();
    };
} // namespace elasty

#endif // ELASTY_ENGINE_HPP
