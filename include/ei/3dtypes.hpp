#pragma once

#include "vector.hpp"

namespace ei {

    // Declarations for all types to enable conversion operators.
    struct Sphere;
    struct Plane;
    struct DOP;
    struct Box;
    struct OBox;
    struct Disc;
    struct Triangle;
    struct Tetrahedron;
    struct Ray;
    struct Segment;
    struct Cone;
    struct Frustum;
    struct Ellipsoid;
    struct OEllipsoid;
    struct Capsule;

    // Fast types are designed if a primitive should be tested against many
    // others. They may take more memory for precomputed and redundant
    // information. For that reason they contain only const members to avoid
    // inconsistent states.
    struct FastFrustum;

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
        Sphere() noexcept {}

        /// \brief Create sphere from center and radius
        Sphere( const Vec3& _center, float _radius ) noexcept;                 // TESTED

        /// \brief Create the bounding sphere of a box
        explicit Sphere( const Box& _box ) noexcept;                           // TESTED

        /// \brief Create the bounding sphere for two points
        Sphere( const Vec3& _p0, const Vec3& _p1 ) noexcept;

        /// \brief Create the bounding sphere for three points
        Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2 ) noexcept;

        /// \brief Create the bounding sphere for four points
        Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const Vec3& _p3 ) noexcept;

        /// \brief Create the bounding sphere for n points using Welzl's
        ///     algorithm.
        /// \details The algorithm has expected linear run time.
        Sphere( const Vec3* _points, uint32 _numPoints ) noexcept;
    };

    /// \brief Axis aligned box.
    /// \details Box invariant: the max must always be larger or equal than min.
    struct Box
    {
        Vec3 min;
        Vec3 max;

        /// \brief Create uninitialized box.
        Box() noexcept {}

        /// \brief Create from minimal and maximal coordinates
        //Box( const Vec3& _min, const Vec3& _max );                             // TESTED

        /// \brief Create a box for a single point.
        /// \details This is also the recursion and for the point list constructor.
        explicit Box( const Vec3& _point ) noexcept;

        template<typename... Args>
        Box( const Vec3& _point, Args... _morePoints ) noexcept;

        /// \brief Get the smallest box containing two boxes.
        Box( const Box& _box0, const Box& _box1 ) noexcept;                    // TESTED

        /// \brief Create the bounding box for a sphere.
        explicit Box( const Sphere& _sphere ) noexcept;                        // TESTED

        /// \brief Create the bounding box for a triangle.
        explicit Box( const Triangle& _triangle ) noexcept;                    // TESTED

        /// \brief Create the bounding box of a tetrahedron
        explicit Box( const Tetrahedron& _tetrahedron ) noexcept;

        /// \brief Create the bounding box for an ellipsoid
        explicit Box( const Ellipsoid& _ellipsoid ) noexcept;                  // TESTED

        /// \brief Create the bounding box for an oriented box.
        explicit Box( const OBox& _box ) noexcept;                             // TESTED

        /// \brief Create an optimal box for a set of points
        Box( const Vec3* _points, uint32 _numPoints ) noexcept;                // TESTED
    };

    /// \brief Oriented bounding box.
    /// \details If you are going to use oriented bounding boxes you might want
    ///     to use multiple double oriented planes (k-DOPs) instead. The
    ///     parametrization of the OBox is more natural whereas k-DOPs are faster????
    struct OBox
    {
        Vec3 center;
        Vec3 halfSides;             ///< Half side lengths of the box.
        Quaternion orientation;     ///< Orientation of the box. This gives the rotation from an AABox to the rotated box.

        /// \brief Create uninitialized box.
        OBox() noexcept {}

        /// \brief Create from parametrization
        OBox( const Vec3& _center, const Vec3& _halfSides, const Quaternion& _orientation ) noexcept;

        /// \brief Create an oriented box from a simple box
        explicit OBox( const Box& _box ) noexcept;

        /// \brief Create an oriented box from a disc.
        /// \details Since the disc is isotropic there is one degree of freedom.
        ///     This ambiguity is solved by using the smallest rotation of the
        ///     z-axis towards the disc normal.
        explicit OBox( const Disc& _disc ) noexcept;

        /// \brief Create an oriented box which contains an aabox
        OBox( const Quaternion& _orientation, const Box& _box ) noexcept;
        OBox( const Mat3x3& _orientation, const Box& _box ) noexcept;

        /// \brief Create an oriented box which contains a set of points
        OBox( const Quaternion& _orientation, const Vec3* _points, uint32 _numPoints ) noexcept;

        /// \brief Find the best oriented box by brute force.
        /// \details Uses O(n^4) brute force algorithm. The exact runtime is
        ///     T(n * binomial(n,3)) = T((n^4 - 3n^3 + 2n^2)/6).
        /// \param [in] _points The point set for which the box is searched.
        /// \param [in] _numPoints Size of the point array.
        /// \param [in] _tries Number of random tries for the orientation.
        OBox( const Vec3* _points, uint32 _numPoints ) noexcept;               // TESTED
    };

    struct Tetrahedron
    {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;
        Vec3 v3;

        /// \brief Indexed access to the 4 vertices
        Vec3& v(int _index) noexcept;
        const Vec3& v(int _index) const noexcept;

        /// \brief Create uninitialized tetrahedron.
        Tetrahedron() noexcept {}

        /// \brief Create from four vertices
        Tetrahedron(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, const Vec3& _v3) noexcept;
    };

    /// \brief A triangle in 3D space.
    struct Triangle
    {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;

        /// \brief Indexed access to the 3 vertices
        Vec3& v(int _index) noexcept;
        const Vec3& v(int _index) const noexcept;

        /// \brief Create uninitialized Triangle.
        Triangle() noexcept {}

        /// \brief Create from three vertex coordinates
        Triangle(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2) noexcept;
    };

    /// \brief A 2D circular element in 3D space
    struct Disc
    {
        Vec3 center;        ///< Center/Position of the disc
        Vec3 normal;        ///< Disc normal
        float radius;       ///< Disc radius [0,INF) allowed

        /// \brief Create uninitialized Disc.
        Disc() noexcept {}

        /// \brief Create from parameters
        Disc(const Vec3& _center, const Vec3& _normal, float _radius) noexcept;
    };

    /// \brief A plane in 3D. If you want to use 2 parallel planes use DOPs
    ///     instead.
    struct Plane
    {
        Vec3 n;     ///< The normal on the plane
        float d;    ///< The distance to the origin such that dot(n,x) + d = 0 for all point in the plane

        /// \brief Create uninitialized Plane.
        Plane() noexcept {}

        /// \brief Create a plane from direct parameters
        Plane(const Vec3& _normal, float _d) noexcept;                         // TESTED

        /// \brief Create a plane from a support vector and a direction vector.
        Plane(const Vec3& _normal, const Vec3& _support) noexcept;             // TESTED

        /// \brief Create a plane from three points.
        /// \details Creates the RHS normal for courter-clock-wise sorted
        ///     vertices. With other words: the normal is that of the
        ///     triangle.
        Plane(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2) noexcept;     // TESTED
    };

    /// \brief A double oriented plane (i.e. two parallel planes).
    /// \details DOPs are often used for generalized bounding volumes (k-DOP).
    ///     E.g. an axis aligned 3-DOP is the same as an axis aligned bounding
    ///     box.
    struct DOP
    {
        Vec3 n;     ///< The normal on the first plane
        float d0;   ///< The negative distance to the origin of the first plane
        float d1;   ///< The negative distance to the origin of the second plane. d0 >= d1

        /// \brief Create uninitialized DOP.
        DOP() noexcept {}

        /// \brief Create a DOP from direct parameters
        /// \param [in] _d0 Negated distance from the origin to the first plane
        ///     (-dot(_normal, _support0)).
        /// \param [in] _d0 Negated distance from the origin to the second plane
        ///     (-dot(_normal, _support1)).
        DOP(const Vec3& _normal, float _d0, float _d1) noexcept;

        /// \brief Create a DOP from a direction (normal) and two support
        ///     vectors.
        DOP(const Vec3& _normal, const Vec3& _support0, const Vec3& _support1) noexcept;
    };

    /// \brief An axis aligned ellipsoid.
    struct Ellipsoid
    {
        Vec3 center;
        Vec3 radii;         ///< 3 radii greater 0

        /// \brief Create uninitialized Ellipsoid.
        Ellipsoid() noexcept {}

        /// \brief Create an Ellipsoid from center and radii.
        /// \param [in] _radii The scaling radii. If a radius is <= 1e-30f the
        ///     constructor replaces it with 1e-30f for reasons of stability.
        Ellipsoid(const Vec3& _center, const Vec3& _radii) noexcept;           // TESTED

        /// \brief Create bounding Ellipsoid from axis aligned bounding box.
        explicit Ellipsoid(const Box& _box) noexcept;                          // TESTED
    };

    /// \brief An oriented ellipsoid.
    struct OEllipsoid
    {
        Vec3 center;
        Vec3 radii;         ///< 3 radii greater 0
        Quaternion orientation;

        /// \brief Create uninitialized Ellipsoid.
        OEllipsoid() noexcept {}

        /// \brief Create an Ellipsoid from parametrization.
        /// \param [in] _radii The scaling radii. If a radius is <= 1e-30f the
        ///     constructor replaces it with 1e-30f for reasons of stability.
        OEllipsoid(const Vec3& _center, const Vec3& _radii, const Quaternion& _orientation) noexcept;

        /// \brief Create bounding Ellipsoid from axis aligned box.
        explicit OEllipsoid(const Box& _box) noexcept;

        /// \brief Create bounding Ellipsoid from oriented box.
        explicit OEllipsoid(const OBox& _box) noexcept;
    };

    /// \brief A ray starts in one point and extends to infinity
    struct Ray
    {
        Vec3 origin;        ///< Origin of the ray
        Vec3 direction;     ///< Normalized direction vector

        /// \brief Create uninitialized Ray.
        Ray() noexcept {}

        /// \brief Create Ray from origin and direction.
        /// \param [in] _direction A normalized direction vector.
        ///     The method does no normalization because it could
        ///     be a redundant operation.
        Ray(const Vec3& _origin, const Vec3& _direction) noexcept;
    };

    /// \brief A line segment is the connection between two points
    struct Segment
    {
        Vec3 a;             ///< Start of the line
        Vec3 b;             ///< End of the line

        /// \brief Create uninitialized Line.
        Segment() noexcept {}

        /// \brief Create from two points
        Segment(const Vec3& _a, const Vec3& _b) noexcept;

        /// \brief Create from bounded ray
        /// \param [in] _distance Length of the ray to define the end point of
        ///     the line.
        Segment(const Ray& _ray, float _distance) noexcept;
    };

    /// \brief A cone starts in one point and extends to a perpendicular disk as
    ///     base. The maximum half opening must be smaller than 90°.
    struct Cone
    {
        Ray centralRay;
        float tanTheta;     ///< Tangents of the half opening angle.
        float height;       ///< Distance from origin to the base.

        /// \brief Create an uninitialized Cone.
        Cone() noexcept {}

        /// \brief Create from intuitive parametrization.
        /// \param [in] _tanHalfOpeningAngle Tangents of the angle from the
        ///     central ray to the hull.
        Cone(const Vec3& _origin, const Vec3& _direction, float _tanHalfOpeningAngle, float _height) noexcept
        {
            eiAssert(approx(lensq(_direction), 1.0f), "Expected a normalized direction!");
            tanTheta = _tanHalfOpeningAngle;
            height = _height;
            centralRay.origin = _origin;
            centralRay.direction = _direction;
        }

        /// \brief Create from direct parametrization.
        /// \param [in] _tanHalfOpeningAngle Tangents of the angle from the
        ///     central ray to the hull.
        Cone(const Ray& _ray, float _tanHalfOpeningAngle, float _height) noexcept :
            centralRay(_ray),
            tanTheta(_tanHalfOpeningAngle),
            height(_height)
        {}
    };

    /// \brief A cylinder with hemispherical ends.
    struct Capsule
    {
        Segment seg;        ///< Start and end of the inner line (cylinder center bottom/top)
        float radius;       ///< Size of the boundary (radius of cylinder and hemisphere)

        /// \brief Create uninitialized Capsule.
        Capsule() noexcept {}

        /// \brief Direct create from parameters
        Capsule(const Vec3& _a, const Vec3& _b, float _radius) noexcept;

        /// \brief Create from line and add boundary
        Capsule(const Segment& _line, float _radius) noexcept;
    };

    /// \brief A pyramid frustum with four planes which intersect in one point.
    /// \details There is also a FastFrustum type which should be used for camera
    ///     frustums etc. This type is intended for calculations with a frustum
    ///     not for fast intersection tests.
    struct Frustum
    {
        Vec3 apex;          ///< The origin / tip of the pyramid
        Vec3 up;
        Vec3 direction;
        float l, r;         ///< Left and right distances on the far plane. Assumption: l < r
        float b, t;         ///< Bottom and top distances on the far plane. Assumption: b < t
        float n, f;         ///< Near and far distance. Assumptions: 0 <= n < f

        /// \brief Create from camera like parametrization (LHS)
        /// \param [in] _direction Normalized direction vector.
        /// \param [in] _l Distance to the left plane from center to border on
        ///     the far plane.
        /// \param [in] _r Distance to the right plane from center to border
        ///     on the far plane.
        /// \param [in] _b Distance to the bottom plane from center to border
        ///     on the far plane.
        /// \param [in] _t Distance to the top plane from center to border on
        ///     the far plane.
        Frustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) noexcept;
    };


    // ********************************************************************* //
    //                           FAST VARIANTS                               //
    // ********************************************************************* //

    struct FastFrustum
    {
        const DOP nf;         ///< Parallel near and far planes
        const Plane l;        ///< Left plane (normal points inward)
        const Plane r;        ///< Right plane (normal points inward)
        const Plane b;        ///< Bottom plane (normal points inward)
        const Plane t;        ///< Top plane (normal points inward)
        const Vec3 vertices[8]; ///< All vertices in the orderd: nlb, nlt, nrb, nlt, flb, flt, frb, frt

        /// \brief Construction from dynamic variant
        FastFrustum(const Frustum& _frustum) noexcept;                         // TESTED

        /// \brief Create from standard frustum parametrization
        FastFrustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) noexcept;

        /// \brief Overwrite the current data (auto generation not possible because of const members)
        FastFrustum& operator = (const FastFrustum& _frustum) noexcept;
    };

    struct FastCone
    {
        const Ray centralRay;
        const float cosThetaSq;
        const float height;       ///< Distance from origin to the base.

        /// \brief Construction from dynamic variant.
        FastCone(const Cone & _cone) noexcept :
            centralRay(_cone.centralRay),
            cosThetaSq(1.0f / (1.0f + _cone.tanTheta * _cone.tanTheta)), // cos(atan(x))^2 == 1/(x^2+1)
            height(_cone.height)
        {}

        /// \brief Overwrite the current data (auto generation not possible because of const members)
        FastCone& operator = (const FastCone& _cone) noexcept;
    };

    /// \brief Get the volume of any object.
    float volume(const Sphere& _sphere) noexcept;                              // TESTED
    float volume(const Box& _box) noexcept;                                    // TESTED
    float volume(const OBox& _obox) noexcept;                                  // TESTED
    float volume(const Tetrahedron& _thetrahedron) noexcept;                   // TESTED
    float volume(const Triangle& _triangle) noexcept;                          // TESTED
    float volume(const Disc& _disc) noexcept;                                  // TESTED
    float volume(const Plane& _plane) noexcept;                                // TESTED
    float volume(const DOP& _dop) noexcept;                                    // TESTED
    float volume(const Ellipsoid& _ellipsoid) noexcept;                        // TESTED
    float volume(const OEllipsoid& _oellipsoid) noexcept;                      // TESTED
    float volume(const Ray& _ray) noexcept;                                    // TESTED
    float volume(const Segment& _line) noexcept;                               // TESTED
    float volume(const Cone& _cone) noexcept;                                  // TESTED
    float volume(const Capsule& _capsule) noexcept;                            // TESTED
    float volume(const Frustum& _frustum) noexcept;                            // TESTED

    /// \brief Get the surface area of any object.
    float surface(const Sphere& _sphere) noexcept;                             // TESTED
    float surface(const Box& _box) noexcept;                                   // TESTED
    float surface(const OBox& _obox) noexcept;                                 // TESTED
    float surface(const Tetrahedron& _thetrahedron) noexcept;                  // TESTED
    float surface(const Triangle& _triangle) noexcept;                         // TESTED
    float surface(const Disc& _disc) noexcept;                                 // TESTED
    float surface(const Plane& _plane) noexcept;                               // TESTED
    float surface(const DOP& _dop) noexcept;                                   // TESTED
    float surface(const Ellipsoid& _ellipsoid) noexcept;                       // TESTED
    float surface(const OEllipsoid& _oellipsoid) noexcept;                     // TESTED
    float surface(const Ray& _ray) noexcept;                                   // TESTED
    float surface(const Segment& _line) noexcept;                              // TESTED
    float surface(const Cone& _cone) noexcept;                                 // TESTED
    float surface(const Capsule& _capsule) noexcept;                           // TESTED
    float surface(const Frustum& _frustum) noexcept;                           // TESTED

    /// \brief Transform a box (rotation).
    OBox transform(const Box& _box, const Quaternion& _rotation) noexcept;
    OBox transform(const OBox& _box, const Quaternion& _rotation) noexcept;
    /// \brief Transform a box (translation).
    Box transform(const Box& _box, const Vec3& _translation) noexcept;
    OBox transform(const OBox& _box, const Vec3& _translation) noexcept;
    /// \brief Transform a box (first rotate then translation).
    OBox transform(const Box& _box, const Quaternion& _rotation, const Vec3& _translation) noexcept;
    OBox transform(const OBox& _box, const Quaternion& _rotation, const Vec3& _translation) noexcept;

    /// \brief Compute the centroid for all kinds of bounded geometry.
    /// \details The geometric centroid coincides with the center of mass if
    ///    the mass is uniformly distributed.
    ///
    ///    A geometric composition of simple geometry can be combined via
    ///    sum(center(G_i)*volume(G_i)) / sum(volume(G_i)). Where negative
    ///    volumes can be used to model holes or to subtract the overlapping
    ///    regions again.
    Vec3 center(const Sphere& _sphere) noexcept;                               // TESTED
    Vec3 center(const Box& _box) noexcept;                                     // TESTED
    Vec3 center(const OBox& _obox) noexcept;                                   // TESTED
    Vec3 center(const Tetrahedron& _thetrahedron) noexcept;                    // TESTED
    Vec3 center(const Triangle& _triangle) noexcept;                           // TESTED
    Vec3 center(const Disc& _disc) noexcept;                                   // TESTED
    Vec3 center(const Ellipsoid& _ellipsoid) noexcept;                         // TESTED
    Vec3 center(const Segment& _line) noexcept;                                // TESTED
    Vec3 center(const Cone& _cone) noexcept;                                   // TESTED
    Vec3 center(const Capsule& _capsule) noexcept;                             // TESTED
    Vec3 center(const Frustum& _frustum) noexcept;                             // TESTED
    //Vec3 center(const FastFrustum& _frustum);

    /// \brief Remove all points from the array, which are not part of the
    ///     convex hull.
    /// \details The algorithm moves the points on the convex hull to the front
    ///     of the _points array. The other points are overwritten.
    /// \param [in] _threshold Discard vertices which are closer to the
    ///     previous convex hull than this threshold. This also includes
    ///     duplicates of vertices.
    /// \return Number of points in the convex set (these are the first elements
    ///     in _points after call).
    uint32 convexSet(Vec3* _points, uint32 _numPoints, float _threshold = 0.0f);

    // Include inline implementations
#   include "details/3dtypes.inl"
}
