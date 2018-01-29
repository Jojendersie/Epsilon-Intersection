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
        Sphere( const Vec3& _center, float _radius ) noexcept :                // TESTED
            center(_center),
            radius(_radius)
        {}

        /// \brief Create the bounding sphere of a box
        explicit Sphere( const Box& _box ) noexcept;                           // TESTED

        /// \brief Create the bounding sphere for two points
        Sphere( const Vec3& _p0, const Vec3& _p1 ) noexcept :
            center((_p0 + _p1) * 0.5f),
            radius(len(_p0 - _p1) * 0.5f)
        {}

        /// \brief Create the bounding sphere for three points
        Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2 ) noexcept;

        /// \brief Create the bounding sphere for four points
        Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const Vec3& _p3 ) noexcept;

        /// \brief Create the bounding sphere for n points using Welzl's
        ///     algorithm.
        /// \details The algorithm has expected linear run time.
        Sphere( const Vec3* _points, uint32 _numPoints ) noexcept;
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
        Disc(const Vec3& _center, const Vec3& _normal, float _radius) noexcept :
        center(_center),
            normal(_normal),
            radius(_radius)
        {
            eiAssert(_radius >= 0.0f, "Expected a positive radius!");
        }
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
        explicit Box( const Vec3& _point ) noexcept :
            min(_point),
            max(_point)
        {}

        template<typename... Args>
        Box( const Vec3& _point, Args... _morePoints ) noexcept :
            Box(_morePoints...)
        {
            min = ei::min(min, _point);
            max = ei::max(max, _point);
        }

        /// \brief Get the smallest box containing two boxes.
        Box( const Box& _box0, const Box& _box1 ) noexcept :                   // TESTED
            min(ei::min(_box0.min, _box1.min)),
            max(ei::max(_box0.max, _box1.max))
        {
            eiAssert( max >= min,
                "Minimum coordinates must be smaller or equal the maximum." );
        }

        /// \brief Create the bounding box for a sphere.
        explicit Box( const Sphere& _sphere ) noexcept :                       // TESTED
            min(_sphere.center - _sphere.radius),
            max(_sphere.center + _sphere.radius)
        {
            eiAssertWeak( max >= min,
                "Subtraction or addition of a scalar failed or sphere had negative radius!" );
        }

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
        OBox( const Vec3& _center, const Vec3& _halfSides, const Quaternion& _orientation ) noexcept :
            center(_center),
            halfSides(_halfSides),
            orientation(_orientation)
        {}

        /// \brief Create an oriented box from a simple box
        explicit OBox( const Box& _box ) noexcept :
            center((_box.max + _box.min) * 0.5f),
            halfSides((_box.max - _box.min) * 0.5f),
            orientation(qidentity())
        {}

        /// \brief Create an oriented box from a disc.
        /// \details Since the disc is isotropic there is one degree of freedom.
        ///     This ambiguity is solved by using the smallest rotation of the
        ///     z-axis towards the disc normal.
        explicit OBox( const Disc& _disc ) noexcept :
            center(_disc.center),
            halfSides(_disc.radius, _disc.radius, 0.0f),
            orientation( Vec3(0.0f, 0.0f, 1.0f), _disc.normal)
        {}

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
        Vec3& v(int _index) noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 4 vertices!");
            return reinterpret_cast<Vec3*>(this)[_index];
        }

        const Vec3& v(int _index) const noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 4 vertices!");
            return reinterpret_cast<const Vec3*>(this)[_index];
        }

        /// \brief Create uninitialized tetrahedron.
        Tetrahedron() noexcept {}

        /// \brief Create from four vertices
        Tetrahedron(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, const Vec3& _v3) noexcept :
            v0(_v0),
            v1(_v1),
            v2(_v2),
            v3(_v3)
        {}
    };

    /// \brief A triangle in 3D space.
    struct Triangle
    {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;

        /// \brief Indexed access to the 3 vertices
        Vec3& v(int _index) noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 3, "A triangle only has 3 vertices!");
            return reinterpret_cast<Vec3*>(this)[_index];
        }

        const Vec3& v(int _index) const noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 3, "A triangle only has 3 vertices!");
            return reinterpret_cast<const Vec3*>(this)[_index];
        }

        /// \brief Create uninitialized Triangle.
        Triangle() noexcept {}

        /// \brief Create from three vertex coordinates
        Triangle(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2) noexcept :
            v0(_v0),
            v1(_v1),
            v2(_v2)
        {}
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
        Plane(const Vec3& _normal, float _d) noexcept :                        // TESTED
            n(_normal),
            d(_d)
        {
            eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
        }

        /// \brief Create a plane from a support vector and a direction vector.
        Plane(const Vec3& _normal, const Vec3& _support) noexcept :            // TESTED
            n(_normal),
            d(-dot(_normal, _support))
        {
            eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
        }

        /// \brief Create a plane from three points.
        /// \details Creates the RHS normal for courter-clock-wise sorted
        ///     vertices. With other words: the normal is that of the
        ///     triangle.
        Plane(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2) noexcept      // TESTED
        {
            n = normalize(cross(_v1 - _v0, _v2 - _v0));
            d = -dot(n, _v0);
        }
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
        DOP(const Vec3& _normal, float _d0, float _d1) noexcept :
            n(_normal),
            d0(max(_d0, _d1)),
            d1(min(_d0, _d1))
        {
            eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
        }

        /// \brief Create a DOP from a direction (normal) and two support
        ///     vectors.
        DOP(const Vec3& _normal, const Vec3& _support0, const Vec3& _support1) noexcept :
            DOP(_normal, -dot(_normal, _support0), -dot(_normal, _support1))
        {}
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
        Ellipsoid(const Vec3& _center, const Vec3& _radii) noexcept :          // TESTED
            center(_center),
            radii(max(_radii, Vec3(1e-16f)))
        {}

        /// \brief Create bounding Ellipsoid from axis aligned bounding box.
        explicit Ellipsoid(const Box& _box) noexcept                           // TESTED
        {
            eiAssert( _box.max >= _box.min, "Invalid bounding box." );
            center = (_box.max + _box.min) * 0.5f;
            /// sqrt(n) * side length / 2, where n is the number of dimensions with
            /// an extension (side length 0 allows to generate ellipses or rays)
            Vec3 sideLen = _box.max - _box.min;
            radii = (sqrt((float)sum(neq(sideLen, 0.0f))) * 0.5f) * sideLen;
            radii = max(radii, Vec3(1e-16f));
        }

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
        OEllipsoid(const Vec3& _center, const Vec3& _radii, const Quaternion& _orientation) noexcept :
            center(_center),
            radii(_radii),
            orientation(_orientation)
        {}

        /// \brief Create bounding Ellipsoid from axis aligned box.
        explicit OEllipsoid(const Box& _box) noexcept
        {
            eiAssert( _box.max >= _box.min, "Invalid box." );
            center = (_box.max + _box.min) * 0.5f;
            /// sqrt(n) * side length / 2, where n is the number of dimensions with
            /// an extension (side length 0 allows to generate ellipses or rays)
            Vec3 sideLen = _box.max - _box.min;
            radii = (sqrt((float)sum(neq(sideLen, 0.0f))) * 0.5f) * sideLen;
            radii = max(radii, Vec3(1e-16f));
            orientation = qidentity();
        }

        /// \brief Create bounding Ellipsoid from oriented box.
        explicit OEllipsoid(const OBox& _box) noexcept
        {
            eiAssert( _box.halfSides >= 0.0f, "Invalid box." );
            center = _box.center;
            /// sqrt(n) * side length / 2, where n is the number of dimensions with
            /// an extension (side length 0 allows to generate ellipses or rays)
            radii = sqrt((float)sum(neq(_box.halfSides, 0.0f))) * _box.halfSides;
            radii = max(radii, Vec3(1e-16f));
            orientation = _box.orientation;
        }
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
        Ray(const Vec3& _origin, const Vec3& _direction) noexcept :
            origin(_origin),
            direction(_direction)
        {
            eiAssert(approx(lensq(_direction), 1.0f), "Insert a normalized normal!");
        }
    };

    /// \brief A line segment is the connection between two points
    struct Segment
    {
        Vec3 a;             ///< Start of the line
        Vec3 b;             ///< End of the line

        /// \brief Create uninitialized Line.
        Segment() noexcept {}

        /// \brief Create from two points
        Segment(const Vec3& _a, const Vec3& _b) noexcept :
            a(_a),
            b(_b)
        {}

        /// \brief Create from bounded ray
        /// \param [in] _distance Length of the ray to define the end point of
        ///     the line.
        Segment(const Ray& _ray, float _distance) noexcept :
            a(_ray.origin),
            b(_ray.origin + _ray.direction * _distance)
        {
            eiAssertWeak(approx(lensq(_ray.direction), 1.0f), "The input ray is not normalized!");
        }
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
        Capsule(const Vec3& _a, const Vec3& _b, float _radius) noexcept :
            seg(_a, _b),
            radius(_radius)
        {
            eiAssertWeak(_radius >= 0.0f, "Radius must be positive!");
        }

        /// \brief Create from line and add boundary
        Capsule(const Segment& _line, float _radius) noexcept :
            seg(_line),
            radius(_radius)
        {
            eiAssertWeak(_radius >= 0.0f, "Radius must be positive!");
        }
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
        Frustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) noexcept :
            l(_l), r(_r), b(_b), t(_t), n(_n), f(_f),
            apex(_apex),
            up(_up),
            direction(_direction)
        {
            eiAssert(approx(lensq(_direction), 1.0f), "Insert a normalized direction!");
            eiAssert(approx(lensq(_up), 1.0f), "Insert a normalized up vector!");
            eiAssert(_n < _f && 0 <= _n, "Near and far frustum planes are sorted wrongly.");
            eiAssert(_l < _r, "Left and right frustum planes are sorted wrongly.");
            eiAssert(_b < _t, "Top and bottom frustum planes are sorted wrongly.");
        }
    };


    // ********************************************************************* //
    //                           FAST VARIANTS                               //
    // ********************************************************************* //

    struct FastRay
    {
        const Vec3 origin;
        const Vec3 direction;
        const Vec3 invDirection;        ///< 1/direction

        /// \brief Construction from dynamic ray struct
        FastRay(const Ray & _ray) noexcept :
            origin(_ray.origin),
            direction(_ray.direction),
            invDirection(1.0f / _ray.direction)
        {}
    };

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
        FastFrustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) noexcept :
            FastFrustum(Frustum(_apex, _direction, _up, _l, _r, _b, _t, _n, _f))
        {}

        /// \brief Overwrite the current data (auto generation not possible because of const members)
        FastFrustum& operator = (const FastFrustum& _frustum) noexcept
        {
            const_cast<DOP&>(nf) = _frustum.nf;
            const_cast<Plane&>(l) = _frustum.l;
            const_cast<Plane&>(r) = _frustum.r;
            const_cast<Plane&>(b) = _frustum.b;
            const_cast<Plane&>(t) = _frustum.t;
            return *this;
        }
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
        FastCone& operator = (const FastCone& _cone) noexcept
        {
            const_cast<Ray&>(centralRay) = _cone.centralRay;
            const_cast<float&>(cosThetaSq) = _cone.cosThetaSq;
            const_cast<float&>(height) = _cone.height;
            return *this;
        }
    };

    struct FastTriangle
    {
        const Vec3 v0;
        const Vec3 e01;     ///< Edge from v0 to v1
        const Vec3 e02;     ///< Edge from v0 to v2
        const Vec3 normal;  ///< Normalized normal
        const float area;

        /// \brief Construction from dynamic variant.
        explicit FastTriangle(const Triangle & _triangle) noexcept :
            v0(_triangle.v0),
            e01(_triangle.v1 - _triangle.v0),
            e02(_triangle.v2 - _triangle.v0),
            normal(0.0f), area(0.0f)
        {
            const_cast<Vec3&>(normal) = cross(e01, e02);
            float a2 = len(normal);
            const_cast<Vec3&>(normal) /= a2;
            const_cast<float&>(area) = a2 * 0.5f;
        }

        /// \brief Overwrite the current data (auto generation not possible because of const members)
        FastTriangle& operator = (const FastTriangle& _triangle) noexcept
        {
            const_cast<Vec3&>(v0) = _triangle.v0;
            const_cast<Vec3&>(e01) = _triangle.e01;
            const_cast<Vec3&>(e02) = _triangle.e02;
            const_cast<Vec3&>(normal) = _triangle.normal;
            const_cast<float&>(area) = _triangle.area;
            return *this;
        }
    };


    // ************************************************************************* //
    // VOLUME AND SURFACE METHODS                                                //
    // ************************************************************************* //
    /// \brief Get the volume of any object.
    inline float volume(const Sphere& _sphere) noexcept                        // TESTED
    {
        return 4.0f / 3.0f * PI * _sphere.radius * _sphere.radius * _sphere.radius;
    }

    inline float volume(const Box& _box) noexcept                              // TESTED
    {
        Vec3 size = _box.max - _box.min;
        return size.x * size.y * size.z;
    }

    inline float volume(const OBox& _obox) noexcept                            // TESTED
    {
        return 8.0f * _obox.halfSides.x * _obox.halfSides.y * _obox.halfSides.z;
    }

    inline float volume(const Tetrahedron& _thetrahedron) noexcept             // TESTED
    {
        return dot(_thetrahedron.v3-_thetrahedron.v0, cross(_thetrahedron.v2-_thetrahedron.v0, _thetrahedron.v1-_thetrahedron.v0)) / 6.0f;
    }

    inline float volume(const Triangle&) noexcept                              // TESTED
    {
        return 0.0f;
    }

    inline float volume(const Disc&) noexcept                                  // TESTED
    {
        return 0.0f;
    }

    inline float volume(const Plane&) noexcept                                 // TESTED
    {
        return 0.0f;
    }

    inline float volume(const DOP&) noexcept                                   // TESTED
    {
        return 0.0f;
    }

    inline float volume(const Ellipsoid& _ellipsoid) noexcept                  // TESTED
    {
        return 4.0f / 3.0f * PI * _ellipsoid.radii.x * _ellipsoid.radii.y * _ellipsoid.radii.z;
    }

    inline float volume(const OEllipsoid& _oellipsoid) noexcept                // TESTED
    {
        return 4.0f / 3.0f * PI * _oellipsoid.radii.x * _oellipsoid.radii.y * _oellipsoid.radii.z;
    }

    inline float volume(const Ray&) noexcept                                   // TESTED
    {
        return 0.0f;
    }

    inline float volume(const Segment&) noexcept                               // TESTED
    {
        return 0.0f;
    }

    inline float volume(const Cone& _cone) noexcept                            // TESTED
    {
        float r = _cone.height * _cone.tanTheta;
        return PI / 3.0f * r * r * _cone.height;
    }

    inline float volume(const Capsule& _capsule) noexcept                      // TESTED
    {
        return PI * sq(_capsule.radius) * (_capsule.radius*4.0f/3.0f + len(_capsule.seg.b-_capsule.seg.a));
    }

    float volume(const Frustum& _frustum) noexcept;                            // TESTED

    /// \brief Get the surface area of any object.
    inline float surface(const Sphere& _sphere) noexcept                       // TESTED
    {
        return 4.0f * PI * sq(_sphere.radius);
    }

    inline float surface(const Box& _box) noexcept                             // TESTED
    {
        Vec3 size = _box.max - _box.min;
        return 2.0f * (size.x * size.y + size.x * size.z + size.y * size.z);
    }

    inline float surface(const OBox& _obox) noexcept                           // TESTED
    {
        return 8.0f * (_obox.halfSides.x * _obox.halfSides.y + _obox.halfSides.x * _obox.halfSides.z + _obox.halfSides.y * _obox.halfSides.z);
    }

    inline float surface(const Tetrahedron& _thetrahedron) noexcept            // TESTED
    {
        // Analogous to a triangle (repeated four times)
        Vec3 a = _thetrahedron.v1 - _thetrahedron.v0;
        Vec3 b = _thetrahedron.v2 - _thetrahedron.v0;
        Vec3 c = _thetrahedron.v3 - _thetrahedron.v0;
        Vec3 d = _thetrahedron.v2 - _thetrahedron.v1;
        Vec3 e = _thetrahedron.v3 - _thetrahedron.v1;
        return 0.5f * (len( cross(a, b) )
            + len( cross(a, c) )
            + len( cross(b, c) )
            + len( cross(d, e) ));
    }

    inline float surface(const Triangle& _triangle) noexcept                   // TESTED
    {
        // Heron's formula is much more expensive than cross product because
        // the 3 side lengths must be computed first.
        return len( cross(_triangle.v1 - _triangle.v0, _triangle.v2 - _triangle.v0) ) * 0.5f;
    }

    inline float surface(const Disc& _disc) noexcept                           // TESTED
    {
        return PI * _disc.radius;
    }

    inline float surface(const Plane&) noexcept                                // TESTED
    {
        return INF;
    }

    inline float surface(const DOP&) noexcept                                  // TESTED
    {
        return INF;
    }

    inline float surface(const Ellipsoid& _ellipsoid) noexcept                 // TESTED
    {
        // Use approximation (Knud Thomsen's formula) only! Everything else is a
        // lot larger.
        Vec3 pr( pow(_ellipsoid.radii.x, 1.6075f),
            pow(_ellipsoid.radii.y, 1.6075f),
            pow(_ellipsoid.radii.z, 1.6075f) );
        return 4.0f * PI * pow((pr.x * pr.y + pr.x * pr.z + pr.y * pr.z) / 3.0f, 0.622083981f);
    }

    inline float surface(const OEllipsoid& _oellipsoid) noexcept               // TESTED
    {
        // Use approximation (Knud Thomsen's formula) only! Everything else is a
        // lot larger.
        Vec3 pr( pow(_oellipsoid.radii.x, 1.6075f),
            pow(_oellipsoid.radii.y, 1.6075f),
            pow(_oellipsoid.radii.z, 1.6075f) );
        return 4.0f * PI * pow((pr.x * pr.y + pr.x * pr.z + pr.y * pr.z) / 3.0f, 0.622083981f);
    }

    inline float surface(const Ray&) noexcept                                  // TESTED
    {
        return 0.0f;
    }

    inline float surface(const Segment&) noexcept                              // TESTED
    {
        return 0.0f;
    }

    inline float surface(const Cone& _cone) noexcept                           // TESTED
    {
        float r = _cone.height * _cone.tanTheta;
        return PI * r * (r + sqrt(r*r + _cone.height*_cone.height));
    }

    inline float surface(const Capsule& _capsule) noexcept                     // TESTED
    {
        return 2 * PI * _capsule.radius * (2 * _capsule.radius + len(_capsule.seg.b-_capsule.seg.a));
    }

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

    // ************************************************************************* //
    // CENTROID METHODS                                                          //
    // ************************************************************************* //
    /// \brief Compute the centroid for all kinds of bounded geometry.
    /// \details The geometric centroid coincides with the center of mass if
    ///    the mass is uniformly distributed.
    ///
    ///    A geometric composition of simple geometry can be combined via
    ///    sum(center(G_i)*volume(G_i)) / sum(volume(G_i)). Where negative
    ///    volumes can be used to model holes or to subtract the overlapping
    ///    regions again.
    inline Vec3 center(const Sphere& _sphere) noexcept                         // TESTED
    {
        return _sphere.center;
    }

    inline Vec3 center(const Box& _box) noexcept                               // TESTED
    {
        return (_box.min + _box.max) * 0.5f;
    }

    inline Vec3 center(const OBox& _obox) noexcept                             // TESTED
    {
        return _obox.center;
    }

    inline Vec3 center(const Tetrahedron& _thetrahedron) noexcept              // TESTED
    {
        return (_thetrahedron.v0 + _thetrahedron.v1 + _thetrahedron.v2 + _thetrahedron.v3) / 4.0f;
    }

    inline Vec3 center(const Triangle& _triangle) noexcept                     // TESTED
    {
        return (_triangle.v0 + _triangle.v1 + _triangle.v2) / 3.0f;
    }

    inline Vec3 center(const Disc& _disc) noexcept                             // TESTED
    {
        return _disc.center;
    }

    inline Vec3 center(const Ellipsoid& _ellipsoid) noexcept                   // TESTED
    {
        return _ellipsoid.center;
    }

    inline Vec3 center(const Segment& _line) noexcept                          // TESTED
    {
        return (_line.a + _line.b) * 0.5f;
    }

    inline Vec3 center(const Cone& _cone) noexcept                             // TESTED
    {
        return _cone.centralRay.origin + _cone.centralRay.direction * 0.75f;
    }

    inline Vec3 center(const Capsule& _capsule) noexcept                       // TESTED
    {
        return (_capsule.seg.a + _capsule.seg.b) * 0.5f;
    }

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

    // ********************************************************************* //
    // Dependent implementations (of types which are not knows at            //
    // declaration point)                                                    //
    // ********************************************************************* //
    inline Sphere::Sphere( const Box& _box ) noexcept :
        center((_box.min + _box.max) * 0.5f),
        radius(len(_box.max - _box.min) * 0.5f)
    {
        eiAssert( _box.max >= _box.min, "Invalid bounding box." );
    }


    inline Box::Box( const Triangle& _triangle ) noexcept :
        min(ei::min(_triangle.v0, _triangle.v1, _triangle.v2)),
        max(ei::max(_triangle.v0, _triangle.v1, _triangle.v2))
    {
        eiAssertWeak( max >= min,
            "min() or max() failed for a vector!" );
    }

    inline Box::Box( const Tetrahedron& _tetrahedron ) noexcept :
        min(ei::min(_tetrahedron.v0, _tetrahedron.v1, _tetrahedron.v2, _tetrahedron.v3)),
        max(ei::max(_tetrahedron.v0, _tetrahedron.v1, _tetrahedron.v2, _tetrahedron.v3))
    {
        eiAssertWeak( max >= min,
            "min() or max() failed for a vector!" );
    }

    inline Box::Box( const Ellipsoid& _ellipsoid ) noexcept :
        min(_ellipsoid.center - _ellipsoid.radii),
        max(_ellipsoid.center + _ellipsoid.radii)
    {
        eiAssertWeak(_ellipsoid.radii >= 0.0f, "Invalid ellipsoid!");
    }
}
