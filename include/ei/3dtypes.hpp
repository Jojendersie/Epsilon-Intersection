#pragma once

#include "matrix.hpp"

namespace ei {

    // Predeclarations for all types to enable conversion operators.
    struct Sphere;
    struct Plane;
    struct Box;
    struct OBox;
    struct Disc;
    struct Triangle;
    struct Thetrahedron;
    struct Ray;
    struct Line;
    struct Frustum;
    struct PyramidFrustum;
    struct Ellipsoid;
    struct OEllipsoid;
    struct Capsule;

    /// \brief A sphere in 3D space.
    struct Sphere
    {
        Vec3 center;
        float radius;

        /// \brief Create uninitialized sphere.
        Sphere() {}

        /// \brief Create sphere from center and radius
        Sphere( const Vec3& _center, float _radius );
    };

    /// \brief Axis aligned box.
    struct Box
    {
        Vec3 min;
        Vec3 max;

        /// \brief Create uninitialized box.
        Box() {}

        /// \brief Create from minimal and maximal coordinates
        Box( const Vec3& _min, const Vec3& _max );
    };

    // Include inline implementations
    #include "details/3dtypes.inl"
}
