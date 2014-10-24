#pragma once

#include "3dtypes.hpp"

namespace ei
{
    /// \brief Get the volume of any object.
    float volume( const Sphere& _sphere);                                      // TESTED
    float volume( const Box& _box);                                            // TESTED
    float volume( const Triangle& _triangle);                                  // TESTED

    /// \brief Get the surface area of any object.
    float surface( const Sphere& _sphere);                                     // TESTED
    float surface( const Box& _box);                                           // TESTED
    float surface( const Triangle& _triangle);                                 // TESTED

    // Include inline implementations
#   include "details/3dfunctions.inl"
}