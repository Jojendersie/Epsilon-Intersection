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
        Sphere( const Vec3& _center, float _radius );                          // TESTED

        /// \brief Create the bounding sphere of a box
        explicit Sphere( const Box& _box );                                    // TESTED
    };

    /// \brief Axis aligned box.
    struct Box
    {
        Vec3 min;
        Vec3 max;

        /// \brief Create uninitialized box.
        Box() {}

        /// \brief Create from minimal and maximal coordinates
        Box( const Vec3& _min, const Vec3& _max );                             // TESTED

        /// \brief Get the smallest box containing two boxes.
        Box( const Box& _box0, const Box& _box1 );                             // TESTED

        /// \brief Create the bounding box for a sphere.
        explicit Box( const Sphere& _sphere );                                 // TESTED

        /// \brief Create the bounding box for a triangle.
        explicit Box( const Triangle& _triangle );                             // TESTED

        /// \brief Create the bounding box for an ellipsoid
        explicit Box( const Ellipsoid& _ellipsoid );                           // TESTED
    };

    /// \brief Oriented bounding box.
    /// \details If you are going to use oriented bounding boxes you might want
    ///     to use multiple double oriented planes (k-DOPs) instead. The
    ///     parametrization of the OBox is more natuaral wheras k-DOPs are faster????
    /*struct OBox
    {
        Vec3 center;
        Vec3 sides;                 ///< Side lengths of the box
        Quaternion orientation;     ///< Orientation of the box
    };*/

    struct Thetrahedron
    {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;
        Vec3 v3;

        /// \brief Indexed access to the 4 vertices
        Vec3& v(int _index);
        const Vec3& v(int _index) const;

        /// \brief Create uninitialized thetrahedron.
        Thetrahedron() {}

        /// \brief Create from four vertices
        Thetrahedron(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, const Vec3& _v3);
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
        Triangle(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2);
    };

    /// \brief A 2D cicular element in 3D space
    struct Disc
    {
        Vec3 center;        ///< Center/Position of the disc
        Vec3 normal;        ///< Disc normal
        float radius;       ///< Disc radius [0,INF) allowed

        /// \brief Create uninitialized Disc.
        Disc() {}

        /// \brief Create from parameters
        Disc(const Vec3& _center, const Vec3& _normal, float _radius);
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

    /// \brief A double oriented plane (i.e. two parallel planes).
    /// \details DOPs are often used for generalized bounding volumes (k-DOP).
    ///     E.g. an axis aligned 3-DOP is the same as an axis aligned bounding
    ///     box.
    struct DOP
    {
        Vec3 n;     ///< The normal on the first plane
        float d0;   ///< The distance to the origin of the first plane
        float d1;   ///< The distance to the origin of the second plane

        /// \brief Create uninitialized DOP.
        DOP() {}

        /// \brief Create a DOP from direct parameters
        /// \param [in] _d0 Negated distance from the origin to the first plane
        ///     (-dot(_normal, _support0)).
        /// \param [in] _d0 Negated distance from the origin to the second plane
        ///     (-dot(_normal, _support1)).
        DOP(const Vec3& _normal, float _d0, float _d1);

        /// \brief Create a DOP from a direction (normal) and two support
        ///     vectors.
        DOP(const Vec3& _normal, const Vec3& _support0, const Vec3& _support1);
    };

    /// \brief An axis alignd ellipsoid.
    struct Ellipsoid
    {
        Vec3 center;
        Vec3 radii;         ///< 3 radii greater 0

        /// \brief Create uninitialized Ellipsoid.
        Ellipsoid() {}

        /// \brief Create an Ellipsoid from center and radii.
        /// \param [in] _radii The scaling radii. If a radius is <= 1e-30f the
        ///     constructor replaces it with 1e-30f for reasons of stability.
        Ellipsoid(const Vec3& _center, const Vec3& _radii);                    // TESTED

        /// \brief Create bounding Ellipsoid from axis aligned bounding box.
        explicit Ellipsoid(const Box& _box);                                   // TESTED
    };

    /// \brief A ray starts in one point and extends to infinity
    struct Ray
    {
        Vec3 origin;        ///< Origin of the ray
        Vec3 direction;     ///< Normlized direction vector

        /// \brief Create uninitialized Ray.
        Ray() {}

        /// \brief Create Ray from origin and direction.
        /// \param [in] _direction A normalized direction vector.
        ///     The method does no normalization because it could
        ///     be a redundant operation.
        Ray(const Vec3& _origin, const Vec3& _direction);
    };

    /// \brief A line is the connection between two points
    struct Line
    {
        Vec3 a;             ///< Start of the line
        Vec3 b;             ///< End of the line

        /// \brief Create uninitialized Line.
        Line() {}

        /// \brief Create from two points
        Line(const Vec3& _a, const Vec3& _b);

        /// \brief Create from bounded ray
        /// \param [in] _distance Length of the ray to define the end point of
        ///     the line.
        Line(const Ray& _ray, float _distance);
    };

    /// \brief A cylinder with hemispherical ends.
    struct Capsule
    {
        Vec3 a;             ///< Start of the inner line (cylinder center bottom)
        Vec3 b;             ///< End of the inner line (cylinder center top)
        float radius;       ///< Size of the boundary (radius of cylinder and hemisphere)

        /// \brief Create uninitialized Capsule.
        Capsule() {}

        /// \brief Direct create from parameters
        Capsule(const Vec3& _a, const Vec3& _b, float _radius);

        /// \brief Create from line and add boundary
        Capsule(const Line& _line, float _radius);
    };

    /// \brief
    struct Frustum
    {
        DOP nf;         ///< Parallel near and far planes
        Plane l;        ///< left plane
        Plane r;        ///< right plane
        Plane t;        ///< top plane
        Plane b;        ///< bottom plane

        /// \brief Create uninitialized Frustum.
        Frustum() {}

        /// \brief Create from camera like parametrization (rhs)
        /// \param [in] _direction Normalized direction vector.
        /// \param [in] _l Distance to the left plane from center to border on the near plane.
        /// \param [in] _r Distance to the right plane from center to border on the near plane.
        /// \param [in] _t Distance to the top plane from center to border on the near plane.
        /// \param [in] _b Distance to the bottom plane from center to border on the near plane.
        Frustum(const Vec3& _origin, const Vec3& _direction, const Vec3& _up, float _n, float _f, float _l, float _r, float _t, float _b);
    };

    // Include inline implementations
    #include "details/3dtypes.inl"
}
