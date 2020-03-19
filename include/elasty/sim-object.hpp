#ifndef sim_object_hpp
#define sim_object_hpp

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

#endif /* sim_object_hpp */
