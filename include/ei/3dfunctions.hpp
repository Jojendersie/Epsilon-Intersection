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
    float volume(const Segment& _line);                                        // TESTED
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
    float surface(const Segment& _line);                                       // TESTED
    float surface(const Capsule& _capsule);                                    // TESTED
    float surface(const Frustum& _frustum);                                    // TESTED

    /// \brief Get the euclidean distance between two objects.
    /// \details The distance for point-solid queries can be negative. All
    ///     other geometries return 0 if they intersect.
    float distance(const Vec3& _point0, const Vec3& _point1);                  // TESTED
    float distance(const Vec3& _point, const Segment& _line);                  // TESTED
    float distance(const Vec3& _point, const Triangle& _triangle);             // TESTED
    float distance(const Vec3& _point, const Sphere& _sphere);                 // TESTED
    float distance(const Vec3& _point, const Capsule& _capsule);               // TESTED
    float distance(const Vec3& _point, const Box& _box);                       // TESTED
    float distance(const Vec3& _point, const Plane& _plane);                   // TESTED
    /// \returns -x if dot(n,_point) larger then both, x if smaller and 0 if it
    ///     is between both planes
    float distance(const Vec3& _point, const DOP& _dop);                   
    float distance(const Sphere& _sphere, const Segment& _segment);            // TESTED
    float distance(const Sphere& _sphere, const Capsule& _capsule);            // TESTED
    float distance(const Sphere& _sphere, const Box& _box);                    // TESTED
    float distance(const Sphere& _sphere, const Plane& _plane);                // TESTED
    float distance(const Segment& _line0, const Segment& _line1);              // TESTED
    float distance(const Capsule& _capsule0, const Capsule& _capsule1);        // TESTED
    inline float distance(const Segment& _line, const Vec3& _point)            { return distance(_point, _line); }
    inline float distance(const Triangle& _triangle, const Vec3& _point)       { return distance(_point, _triangle); }
    inline float distance(const Sphere& _sphere, const Vec3& _point)           { return distance(_point, _sphere); }
    inline float distance(const Capsule& _capsule, const Vec3& _point)         { return distance(_point, _capsule); }
    inline float distance(const Capsule& _capsule, const Sphere& _sphere)      { return distance(_sphere, _capsule); }
    inline float distance(const Box& _box, const Vec3& _point)                 { return distance(_point, _box); }
    inline float distance(const Box& _box, const Sphere& _sphere)              { return distance(_sphere, _box); }
    inline float distance(const Plane& _plane, const Vec3& _point)             { return distance(_point, _plane); }
    inline float distance(const Plane& _plane, const Sphere& _sphere)          { return distance(_sphere, _plane); }
    inline float distance(const DOP& _dop, const Vec3& _point)                 { return distance(_point, _dop); }

    // Include inline implementations
#   include "details/3dfunctions.inl"
}