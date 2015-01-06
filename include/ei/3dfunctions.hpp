#pragma once

#include "3dtypes.hpp"

namespace ei
{
    /// \brief Get the volume of any object.
    float volume(const Sphere& _sphere);                                       // TESTED
    float volume(const Box& _box);                                             // TESTED
    float volume(const Thetrahedron& _thetrahedron);                           // TESTED
    float volume(const Triangle& _triangle);                                   // TESTED
    float volume(const Disc& _disc);                                           // TESTED
    float volume(const Plane& _plane);                                         // TESTED
    float volume(const DOP& _dop);                                             // TESTED
    float volume(const Ellipsoid& _ellipsoid);                                 // TESTED
    float volume(const Ray& _ray);                                             // TESTED
    float volume(const Line& _line);                                           // TESTED
    float volume(const Capsule& _capsule);                                     // TESTED
    float volume(const Frustum& _frustum);                                     // TESTED

    /// \brief Get the surface area of any object.
    float surface(const Sphere& _sphere);                                      // TESTED
    float surface(const Box& _box);                                            // TESTED
    float surface(const Thetrahedron& _thetrahedron);                          // TESTED
    float surface(const Triangle& _triangle);                                  // TESTED
    float surface(const Disc& _disc);                                          // TESTED
    float surface(const Plane& _plane);                                        // TESTED
    float surface(const DOP& _dop);                                            // TESTED
    float surface(const Ellipsoid& _ellipsoid);                                // TESTED
    float surface(const Ray& _ray);                                            // TESTED
    float surface(const Line& _line);                                          // TESTED
    float surface(const Capsule& _capsule);                                    // TESTED
    float surface(const Frustum& _frustum);                                    // TESTED

    // Include inline implementations
#   include "details/3dfunctions.inl"
}