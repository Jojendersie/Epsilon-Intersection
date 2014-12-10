#pragma once

#include "3dtypes.hpp"

namespace ei
{
    /// \brief Get the volume of any object.
    float volume(const Sphere& _sphere);                                       // TESTED
    float volume(const Box& _box);                                             // TESTED
    float volume(const Triangle& _triangle);                                   // TESTED
    float volume(const Disc& _disc);
    float volume(const Plane& _plane);                                         // TESTED
    float volume(const DOP& _dop);
    float volume(const Ellipsoid& _ellipsoid);                                 // TESTED

    /// \brief Get the surface area of any object.
    float surface(const Sphere& _sphere);                                      // TESTED
    float surface(const Box& _box);                                            // TESTED
    float surface(const Triangle& _triangle);                                  // TESTED
    float surface(const Disc& _disc);
    float surface(const Plane& _plane);                                        // TESTED
    float surface(const DOP& _dop);
    float surface(const Ellipsoid& _ellipsoid);                                // TESTED

    // Include inline implementations
#   include "details/3dfunctions.inl"
}