#ifndef ELASTY_FEM_HPP
#define ELASTY_FEM_HPP

#include <Eigen/Core>

// References:
//
// [1] Eftychios Sifakis and Jernej Barbic. 2012. FEM simulation of 3D deformable solids: a practitioner's guide to
// theory, discretization and model reduction. In ACM SIGGRAPH 2012 Courses (SIGGRAPH '12). Association for Computing
// Machinery, New York, NY, USA, Article 20, 1–50. DOI:https://doi.org/10.1145/2343483.2343501
//
// [2] Theodore Kim and David Eberle. 2020. Dynamic deformables: implementation and production practicalities. In ACM
// SIGGRAPH 2020 Courses (SIGGRAPH '20). Association for Computing Machinery, New York, NY, USA, Article 23, 1–182.
// DOI:https://doi.org/10.1145/3388769.3407490

namespace elasty::fem
{
    // Reference: [1]
    template <typename Scalar> constexpr Scalar calcFirstLame(const Scalar youngs_modulus, const Scalar poisson_ratio)
    {
        return youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    }

    // Reference: [1]
    template <typename Scalar> constexpr Scalar calcSecondLame(const Scalar youngs_modulus, const Scalar poisson_ratio)
    {
        return youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    }
} // namespace elasty::fem

#endif // ELASTY_FEM_HPP
