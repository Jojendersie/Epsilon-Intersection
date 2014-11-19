#pragma once

#include "matrix.hpp"

namespace ei {

    // Predeclarations for all types to enable conversion operators.
    struct Sphere;
    struct Plane;
	struct DOP;
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

    /// \brief A list of all supported 3d types
    enum struct Types3D
    {
        SPHERE,
        PLANE,
        DOUBLE_PLANE,
        BOX,
        OBOX,
        DISC,
        TRIANGLE,
        THETRAHEDRON,
        RAY,
        LINE,
        FRUSTUM,
        PYRAMIDFRUSTUM,
        ELLIPSOID,
        OELLIPSOID,
        CAPSULE,

        NUM_TYPES
    };

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

        /// \brief Get the smallest box containing two boxes.
        Box( const Box& _box0, const Box& _box1 );

        /// \brief Create the bounding box for a sphere.
        explicit Box( const Sphere& _sphere );

        /// \brief Create the bounding box for a triangle.
        explicit Box( const Triangle& _triangle );
    };

    /// \brief A triangle in 3D space.
    struct Triangle
    {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;

        /// \brief Indexed access to the 3 vertices
        Vec3& v(int _index);
        const Vec3& v(int _index) const;

        /// \brief Create uninitialized Triangle.
        Triangle() {}

        /// \brief Create from three vertex coordinates
        Triangle( const Vec3& _v0, const Vec3& _v1, const Vec3& _v2 );
    };

    /// \brief A plane in 3D. If you want to use 2 parallel planes use DOPs
    ///     instead.
    struct Plane
    {
        Vec3 n;     ///< The normal on the plane
        float d;    ///< The distance to the origin

        /// \brief Create uninitialized Plane.
        Plane() {}

        /// \brief Create a plane from direct parameters
        Plane(const Vec3& _normal, float _d);                                  // TESTED

        /// \brief Create a plane from a support vector and a direction vector.
        Plane(const Vec3& _normal, const Vec3& _support);                      // TESTED

        /// \brief Create a plane from three points.
        /// \details Creates the RHS normal for courter-clock-wise sorted
        ///     vertices. With other words: the normal is that of the
        ///     triangle.
        Plane(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2);              // TESTED
    };

    struct Ellipsoid
    {
        Vec3 center;
        Vec3 radii;

        /// \brief Create uninitialized Ellipsoid.
        Ellipsoid() {}

        /// \brief Create an Ellipsoid from center and radii.
        Ellipsoid(const Vec3& _center, const Vec3& _radii);                    // TESTED

        /// \brief Create bounding Ellipsoid from axis aligned bounding box.
        Ellipsoid(const Box& _box);                                            // TESTED
    };

    // Include inline implementations
    #include "details/3dtypes.inl"
}
