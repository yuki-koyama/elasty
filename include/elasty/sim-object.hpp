#ifndef sim_object_hpp
#define sim_object_hpp

#include <memory>

namespace elasty
{
    class Particle;
    class Constraint;

    class SimObject
    {
    public:

        std::vector<std::shared_ptr<Particle>> m_particles;
        std::vector<std::shared_ptr<Constraint>> m_constraints;
    };
}

#endif /* sim_object_hpp */
