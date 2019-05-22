#ifndef particle_hpp
#define particle_hpp

#include <Eigen/Core>

namespace elasty
{
    struct Particle
    {
        Particle(const Eigen::Vector3d& x,
                 const Eigen::Vector3d& v,
                 const double           m)
            : x(x), v(v), p(x), f(Eigen::Vector3d::Zero()), m(m), w(1.0 / m)
        {
        }

        Eigen::Vector3d x;
        Eigen::Vector3d v;
        Eigen::Vector3d p;
        Eigen::Vector3d f;
        double          m;
        double          w;
    };
} // namespace elasty

#endif /* particle_hpp */
