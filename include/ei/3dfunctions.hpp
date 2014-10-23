#pragma once

#include "3dtypes.hpp"

namespace ei
{
    /// \brief Get the volume of any object.
    float volume( const Sphere& _sphere);
    float volume( const Box& _box);

    /// \brief Get the surface area of any object.
    float surface( const Sphere& _sphere);
    float surface( const Box& _box);

    // Include inline implementations
#   include "details/3dfunctions.inl"
}