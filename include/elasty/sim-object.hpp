#ifndef ELASTY_SIM_OBJECT_HPP
#define ELASTY_SIM_OBJECT_HPP

#include <memory>
#include <vector>

namespace elasty
{
    struct Particle;
    class AbstractConstraint;

    class SimObject
    {
    public:
        std::vector<std::shared_ptr<Particle>>           m_particles;
        std::vector<std::shared_ptr<AbstractConstraint>> m_constraints;
    };
} // namespace elasty

#endif // ELASTY_SIM_OBJECT_HPP
