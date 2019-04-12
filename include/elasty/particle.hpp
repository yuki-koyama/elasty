#ifndef particle_hpp
#define particle_hpp

#include <cassert>
#include <memory>
#include <vector>
#include <Eigen/Core>

namespace elasty
{
    struct Particle
    {
        Eigen::Vector3d x;
        Eigen::Vector3d v;
        Eigen::Vector3d p;
        double m;
        unsigned int i;
    };
}

#endif /* particle_hpp */
