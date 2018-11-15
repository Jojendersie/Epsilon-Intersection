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
        EIAPI Sphere() noexcept {}

        /// \brief Create sphere from center and radius
        EIAPI Sphere( const Vec3& _center, float _radius ) noexcept :                // TESTED
            center(_center),
            radius(_radius)
        {}

        /// \brief Create the bounding sphere of a box
        EIAPI explicit Sphere( const Box& _box ) noexcept;                           // TESTED

        /// \brief Create the bounding sphere for two points
        EIAPI Sphere( const Vec3& _p0, const Vec3& _p1 ) noexcept :
            center((_p0 + _p1) * 0.5f),
            radius(len(_p0 - _p1) * 0.5f)
        {}

        /// \brief Create the bounding sphere for three points
        EIAPI Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2 ) noexcept
        {
            // The center of the circumscribed circle is at (barycentric coordinates)
            // v0*sin(2 alpha) + v1*sin(2 beta) + v2*sin(2 gamma) and has the radius
            // abc/4A.
            Vec3 c = _p0 - _p1;    float csq = lensq(c);
            Vec3 a = _p1 - _p2;    float asq = lensq(a);
            Vec3 b = _p2 - _p0;    float bsq = lensq(b);

            // One of the sides could be the longest side - the minimum sphere is
            // defined through only two points.
            // This can also handle the coplanar case.
            if( csq + bsq <= asq ) *this = Sphere(_p1, _p2);
            else if( asq + bsq <= csq ) *this = Sphere(_p1, _p0);
            else if( asq + csq <= bsq ) *this = Sphere(_p2, _p0);
            else {
                float area2Sq = 2*lensq(cross(a, c));
                center =
                         _p0 * (-dot(c,b)*asq/area2Sq)
                       + _p1 * (-dot(c,a)*bsq/area2Sq)
                       + _p2 * (-dot(b,a)*csq/area2Sq);
                radius = sqrt(asq*bsq*csq/(2*area2Sq));
            }
        }

        /// \brief Create the bounding sphere for four points
        EIAPI Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const Vec3& _p3 ) noexcept
        {
            // It is possible that not all 4 points lie on the surface of the sphere.
            // Just two of them could already define a sphere enclosing all other.
            // So we need to compute any combination of possible spheres (14), but
            // luckily we know a direct solution for any combination of 3 points.
            // This reduces the work to 4 cases: build a bounding sphere for 3 points
            // and have a look if the fourth point is inside.
            *this = Sphere( _p1, _p2, _p3 );
            float rsq = radius*radius;
            if( lensq(_p3 - center) > rsq ) {
                *this = Sphere( _p0, _p1, _p3 );
                if( lensq(_p2 - center) > rsq ) {
                    *this = Sphere( _p0, _p2, _p3 );
                    if( lensq(_p1 - center) > rsq ) {
                        *this = Sphere( _p1, _p2, _p3 );
                        if( lensq(_p0 - center) > rsq ) {
                            // All 4 points are on the boundary -> construct sphere
                            // from 4 points.
                            Vec3 a = _p1 - _p0;
                            Vec3 b = _p2 - _p0;
                            Vec3 c = _p3 - _p0;

                            Mat3x3 m = Mat3x3( _p1[0], _p1[1], _p1[2],
                                               _p2[0], _p2[1], _p2[2],
                                               _p3[0], _p3[1], _p3[2] );

                            float denominator = 0.5f / determinant( m );

                            Vec3 o = (lensq(c) * cross(a,b) +
                                      lensq(b) * cross(c,a) +
                                      lensq(a) * cross(b,c)) * denominator;

                            center = _p0 + o;
                            radius = len(o);
                        }
                    }
                }
            }
        }

        /// \brief Create the bounding sphere for n points using Welzl's
        ///     algorithm.
        /// \details The algorithm has expected linear run time.
        EIAPI Sphere( const Vec3* _points, uint32 _numPoints ) noexcept
        {
            eiAssert( _numPoints > 0, "The point list must have at least one point." );
            // Create a single linked list for the move to front heuristic
            SingleLinkedPointList* list = new SingleLinkedPointList[_numPoints];
            for(uint32 i = 0; i < _numPoints; ++i) {
                list[i].p = _points[i];
                list[i].next = i+1;
            }
            list[_numPoints-1].next = -1;
            *this = minimalBoundingSphere(list, 0, _numPoints, 1);
            delete[] list;
        }
    private:
        struct SingleLinkedPointList
        {
            Vec3 p;
            int next;   ///< -1 for the last element
        };

        EIAPI static Sphere minimalBoundingSphere( SingleLinkedPointList* _points, uint32 _first, uint32 _n, uint32 _boundarySet )
        {
            Sphere mbs;
            eiAssertWeak(_boundarySet > 0, "Expected at least one point.");

            // If the boundary list is full or all points where added stop
            switch(_boundarySet)
            {
            case 1: mbs = Sphere(_points[_first].p, 0.0f);
                break;
            case 2: {
                Vec3 v0 = _points[_first].p; uint32 next = _points[_first].next;
                mbs = Sphere(v0, _points[next].p);
                break;
            }
            case 3: {
                Vec3 v0 = _points[_first].p; uint32 next = _points[_first].next;
                Vec3 v1 = _points[next].p; next = _points[next].next;
                mbs = Sphere(v0, v1, _points[next].p);
                break;
            }
            case 4: {
                Vec3 v0 = _points[_first].p; uint32 next = _points[_first].next;
                Vec3 v1 = _points[next].p; next = _points[next].next;
                Vec3 v2 = _points[next].p; next = _points[next].next;
                return Sphere(v0, v1, v2, _points[next].p);
            }
            }

            uint32 it = _first;
            uint32 last = _first; // At the beginning this is wrong but does not cause damage, afterwards it will be ovewritten by real numbers
            for(uint32 i = 0; i < _boundarySet; ++i) {last = it; it = _points[it].next;}
            for(uint32 i = _boundarySet; i < _n; ++i)
            {
                // Save next pointer to advance from this point even if the list is changed
                uint32 next = _points[it].next;
                eiAssert(it != 0xffffffff, "Iteration should not have stopped.");
                if(lensq(mbs.center - _points[it].p) > mbs.radius * mbs.radius)
                {
                    // Move to front
                    _points[last].next = _points[it].next;
                    _points[it].next = _first;
                    _first = it;
                    // Rebuild the first i elements
                    mbs = minimalBoundingSphere(_points, it, i, _boundarySet+1);
                }
                last = it;
                it = next;
            }

            return mbs;
        }
    };

    /// \brief A 2D circular element in 3D space
    struct Disc
    {
        Vec3 center;        ///< Center/Position of the disc
        Vec3 normal;        ///< Disc normal
        float radius;       ///< Disc radius [0,INF) allowed

                            /// \brief Create uninitialized Disc.
        EIAPI Disc() noexcept {}

        /// \brief Create from parameters
        EIAPI Disc(const Vec3& _center, const Vec3& _normal, float _radius) noexcept :
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
        EIAPI Box() noexcept {}

        /// \brief Create from minimal and maximal coordinates
        //Box( const Vec3& _min, const Vec3& _max );                             // TESTED

        /// \brief Create a box for a single point.
        /// \details This is also the recursion and for the point list constructor.
        EIAPI explicit Box( const Vec3& _point ) noexcept :
            min(_point),
            max(_point)
        {}

        template<typename... Args>
        EIAPI Box( const Vec3& _point, Args... _morePoints ) noexcept :
            Box(_morePoints...)
        {
            min = ei::min(min, _point);
            max = ei::max(max, _point);
        }

        /// \brief Get the smallest box containing two boxes.
        EIAPI Box( const Box& _box0, const Box& _box1 ) noexcept :                   // TESTED
            min(ei::min(_box0.min, _box1.min)),
            max(ei::max(_box0.max, _box1.max))
        {
            eiAssert( max >= min,
                "Minimum coordinates must be smaller or equal the maximum." );
        }

        /// \brief Create the bounding box for a sphere.
        EIAPI explicit Box( const Sphere& _sphere ) noexcept :                       // TESTED
            min(_sphere.center - _sphere.radius),
            max(_sphere.center + _sphere.radius)
        {
            eiAssertWeak( max >= min,
                "Subtraction or addition of a scalar failed or sphere had negative radius!" );
        }

        /// \brief Create the bounding box for a triangle.
        EIAPI explicit Box( const Triangle& _triangle ) noexcept;                    // TESTED

        /// \brief Create the bounding box of a tetrahedron
        EIAPI explicit Box( const Tetrahedron& _tetrahedron ) noexcept;

        /// \brief Create the bounding box for an ellipsoid
        EIAPI explicit Box( const Ellipsoid& _ellipsoid ) noexcept;                  // TESTED

        /// \brief Create the bounding box for an oriented box.
        EIAPI explicit Box( const OBox& _box ) noexcept;                             // TESTED

        /// \brief Create an optimal box for a set of points
        EIAPI Box( const Vec3* _points, uint32 _numPoints ) noexcept                 // TESTED
        {
            eiAssert( _points && _numPoints > 0, "The point list must have at least one point." );
            min = max = *_points++;
            for( uint32 i = 1; i < _numPoints; ++i, ++_points )
            {
                min = ei::min(min, *_points);
                max = ei::max(max, *_points);
            }
        }
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
        EIAPI OBox() noexcept {}

        /// \brief Create from parametrization
        EIAPI OBox( const Vec3& _center, const Vec3& _halfSides, const Quaternion& _orientation ) noexcept :
            center(_center),
            halfSides(_halfSides),
            orientation(_orientation)
        {}

        /// \brief Create an oriented box from a simple box
        EIAPI explicit OBox( const Box& _box ) noexcept :
            center((_box.max + _box.min) * 0.5f),
            halfSides((_box.max - _box.min) * 0.5f),
            orientation(qidentity())
        {}

        /// \brief Create an oriented box from a disc.
        /// \details Since the disc is isotropic there is one degree of freedom.
        ///     This ambiguity is solved by using the smallest rotation of the
        ///     z-axis towards the disc normal.
        EIAPI explicit OBox( const Disc& _disc ) noexcept :
            center(_disc.center),
            halfSides(_disc.radius, _disc.radius, 0.0f),
            orientation( Vec3(0.0f, 0.0f, 1.0f), _disc.normal)
        {}

        /// \brief Create an oriented box which contains an aabox
        EIAPI OBox( const Quaternion& _orientation, const Box& _box ) noexcept :
            center((_box.min + _box.max) * 0.5f),
            orientation(_orientation)
        {
            // Use the same farthest plane search like in Box(OBox), but using
            // columns instead rows.
            Mat3x3 rotation(~_orientation);
            Vec3 bmax = _box.max - center;
            halfSides = abs(rotation) * bmax;
        }
        EIAPI OBox( const Mat3x3& _orientation, const Box& _box ) noexcept :
            center((_box.min + _box.max) * 0.5f),
            orientation(_orientation)
        {
            // Use the same farthest plane search like in Box(OBox), but using
            // columns instead rows.
            Vec3 bmax = _box.max - center;
            halfSides = abs(transpose(_orientation)) * bmax;
        }

        /// \brief Create an oriented box which contains a set of points
        EIAPI OBox( const Quaternion& _orientation, const Vec3* _points, uint32 _numPoints ) noexcept :
            orientation(_orientation)
        {
            eiAssert( _numPoints > 0, "The point list must have at least one point." );

            // Project all points to the cube sides by transforming them into local
            // space, such that the box is axis aligned again.
            Mat3x3 rotation(conjugate(_orientation));
            Vec3 min, max;
            min = max = rotation * _points[0];
            for(uint32 i = 1; i < _numPoints; ++i)
            {
                Vec3 p = rotation * _points[i];
                min = ei::min(min, p);
                max = ei::max(max, p);
            }

            // Center known with respect to local rotation, go back to world space.
            center = transform((min + max) * 0.5f, conjugate(_orientation));
            halfSides = (max - min) * 0.5f;
        }

        /// \brief Find the best oriented box by brute force.
        /// \details Uses O(n^4) brute force algorithm. The exact runtime is
        ///     T(n * binomial(n,3)) = T((n^4 - 3n^3 + 2n^2)/6).
        /// \param [in] _points The point set for which the box is searched.
        /// \param [in] _numPoints Size of the point array.
        /// \param [in] _tries Number of random tries for the orientation.
        EIAPI OBox( const Vec3* _points, uint32 _numPoints ) noexcept                // TESTED
        {
            if(_numPoints == 1)
            {
                halfSides = Vec3(0.0f);
                center = *_points;
                orientation = qidentity();
            } else if(_numPoints == 2) {
                Vec3 connection = _points[1] - _points[0];
                halfSides = Vec3(len(connection) * 0.5f, 0.0f, 0.0f);
                orientation = Quaternion(connection/halfSides.x, Vec3(1.0f, 0.0f, 0.0f));
                center = _points[0] + 0.5f * connection;
            } else {
                halfSides = Vec3(1e12f);
                float volume = 1e36f;
                Mat3x3 finalRotation;
                // Try each combination of three vertices to setup an orientation
                for(uint32 i = 0; i < _numPoints-2; ++i)
                {
                    for(uint32 j = i+1; j < _numPoints-1; ++j)
                    {
                        Vec3 xAxis = normalize(_points[i] - _points[j]);
                        for(uint32 k = j+1; k < _numPoints; ++k)
                        {
                            Mat3x3 rotation;
                            Vec3 yAxis = cross(xAxis, _points[i] - _points[k]);
                            float l = len(yAxis);
                            if( l < 1e-6f ) // Colinear points
                                rotation = Mat3x3(Quaternion(xAxis, Vec3(1.0f, 0.0f, 0.0f)));
                            else {
                                yAxis /= l;
                                rotation = Mat3x3(xAxis, yAxis, cross(xAxis, yAxis));
                            }
                            // Refit a box with the current rotation. Since the rotation of
                            // all points is the most expensive part try to early out.
                            Vec3 min, max;
                            min = max = rotation * _points[0];
                            for(uint32 a = 1; a < _numPoints; ++a)
                            {
                                Vec3 p = rotation * _points[a];
                                if(p.x < min.x) min.x = p.x;
                                if(p.y < min.y) min.y = p.y;
                                if(p.z < min.z) min.z = p.z;
                                if(p.x > max.x) max.x = p.x;
                                if(p.y > max.y) max.y = p.y;
                                if(p.z > max.z) max.z = p.z;
                                if(prod(max-min) >= volume) goto NextRotation;
                            }
                            // The new box is better than the current, otherwise
                            // the goto would have skipped this section before.
                            //center = transform((min + max) * 0.5f, conjugate(testOrientation));
                            center = transpose(rotation) * ((min + max) * 0.5f);
                            halfSides = max - min;
                            finalRotation = rotation;
                            //orientation = conjugate(Quaternion(rotation));
                            volume = prod(halfSides);
                            NextRotation:;
                        }
                    }
                }
                orientation = conjugate(Quaternion(finalRotation));
                halfSides *= 0.5f;
            }
        }
    };

    struct Tetrahedron
    {
        Vec3 v0;
        Vec3 v1;
        Vec3 v2;
        Vec3 v3;

        /// \brief Indexed access to the 4 vertices
        EIAPI Vec3& v(int _index) noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 4 vertices!");
            return reinterpret_cast<Vec3*>(this)[_index];
        }

        EIAPI const Vec3& v(int _index) const noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 4, "A thetrahedron only has 4 vertices!");
            return reinterpret_cast<const Vec3*>(this)[_index];
        }

        /// \brief Create uninitialized tetrahedron.
        EIAPI Tetrahedron() noexcept {}

        /// \brief Create from four vertices
        EIAPI Tetrahedron(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, const Vec3& _v3) noexcept :
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
        EIAPI Vec3& v(int _index) noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 3, "A triangle only has 3 vertices!");
            return reinterpret_cast<Vec3*>(this)[_index];
        }

        EIAPI const Vec3& v(int _index) const noexcept
        {
            eiAssertWeak(_index >= 0 && _index < 3, "A triangle only has 3 vertices!");
            return reinterpret_cast<const Vec3*>(this)[_index];
        }

        /// \brief Create uninitialized Triangle.
        EIAPI Triangle() noexcept {}

        /// \brief Create from three vertex coordinates
        EIAPI Triangle(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2) noexcept :
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
        EIAPI Plane() noexcept {}

        /// \brief Create a plane from direct parameters
        EIAPI Plane(const Vec3& _normal, float _d) noexcept :                        // TESTED
            n(_normal),
            d(_d)
        {
            eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
        }

        /// \brief Create a plane from a support vector and a direction vector.
        EIAPI Plane(const Vec3& _normal, const Vec3& _support) noexcept :            // TESTED
            n(_normal),
            d(-dot(_normal, _support))
        {
            eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
        }

        /// \brief Create a plane from three points.
        /// \details Creates the RHS normal for courter-clock-wise sorted
        ///     vertices. With other words: the normal is that of the
        ///     triangle.
        EIAPI Plane(const Vec3& _v0, const Vec3& _v1, const Vec3& _v2) noexcept      // TESTED
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
        EIAPI DOP() noexcept {}

        /// \brief Create a DOP from direct parameters
        /// \param [in] _d0 Negated distance from the origin to the first plane
        ///     (-dot(_normal, _support0)).
        /// \param [in] _d0 Negated distance from the origin to the second plane
        ///     (-dot(_normal, _support1)).
        EIAPI DOP(const Vec3& _normal, float _d0, float _d1) noexcept :
            n(_normal),
            d0(max(_d0, _d1)),
            d1(min(_d0, _d1))
        {
            eiAssert(approx(len(_normal), 1.0f), "Expected normalized vector for the normal!");
        }

        /// \brief Create a DOP from a direction (normal) and two support
        ///     vectors.
        EIAPI DOP(const Vec3& _normal, const Vec3& _support0, const Vec3& _support1) noexcept :
            DOP(_normal, -dot(_normal, _support0), -dot(_normal, _support1))
        {}
    };

    /// \brief An axis aligned ellipsoid.
    struct Ellipsoid
    {
        Vec3 center;
        Vec3 radii;         ///< 3 radii greater 0

        /// \brief Create uninitialized Ellipsoid.
        EIAPI Ellipsoid() noexcept {}

        /// \brief Create an Ellipsoid from center and radii.
        /// \param [in] _radii The scaling radii. If a radius is <= 1e-30f the
        ///     constructor replaces it with 1e-30f for reasons of stability.
        EIAPI Ellipsoid(const Vec3& _center, const Vec3& _radii) noexcept :          // TESTED
            center(_center),
            radii(max(_radii, Vec3(1e-16f)))
        {}

        /// \brief Create bounding Ellipsoid from axis aligned bounding box.
        EIAPI explicit Ellipsoid(const Box& _box) noexcept                           // TESTED
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
        EIAPI OEllipsoid() noexcept {}

        /// \brief Create an Ellipsoid from parametrization.
        /// \param [in] _radii The scaling radii. If a radius is <= 1e-30f the
        ///     constructor replaces it with 1e-30f for reasons of stability.
        EIAPI OEllipsoid(const Vec3& _center, const Vec3& _radii, const Quaternion& _orientation) noexcept :
            center(_center),
            radii(_radii),
            orientation(_orientation)
        {}

        /// \brief Create bounding Ellipsoid from axis aligned box.
        EIAPI explicit OEllipsoid(const Box& _box) noexcept
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
        EIAPI explicit OEllipsoid(const OBox& _box) noexcept
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
        EIAPI Ray() noexcept {}

        /// \brief Create Ray from origin and direction.
        /// \param [in] _direction A normalized direction vector.
        ///     The method does no normalization because it could
        ///     be a redundant operation.
        EIAPI Ray(const Vec3& _origin, const Vec3& _direction) noexcept :
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
        EIAPI Segment() noexcept {}

        /// \brief Create from two points
        EIAPI Segment(const Vec3& _a, const Vec3& _b) noexcept :
            a(_a),
            b(_b)
        {}

        /// \brief Create from bounded ray
        /// \param [in] _distance Length of the ray to define the end point of
        ///     the line.
        EIAPI Segment(const Ray& _ray, float _distance) noexcept :
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
        EIAPI Cone() noexcept {}

        /// \brief Create from intuitive parametrization.
        /// \param [in] _tanHalfOpeningAngle Tangents of the angle from the
        ///     central ray to the hull.
        EIAPI Cone(const Vec3& _origin, const Vec3& _direction, float _tanHalfOpeningAngle, float _height) noexcept
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
        EIAPI Cone(const Ray& _ray, float _tanHalfOpeningAngle, float _height) noexcept :
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
        EIAPI Capsule() noexcept {}

        /// \brief Direct create from parameters
        EIAPI Capsule(const Vec3& _a, const Vec3& _b, float _radius) noexcept :
            seg(_a, _b),
            radius(_radius)
        {
            eiAssertWeak(_radius >= 0.0f, "Radius must be positive!");
        }

        /// \brief Create from line and add boundary
        EIAPI Capsule(const Segment& _line, float _radius) noexcept :
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
        EIAPI Frustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) noexcept :
            apex(_apex),
            up(_up),
            direction(_direction),
            l(_l), r(_r), b(_b), t(_t), n(_n), f(_f)
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
        EIAPI FastRay(const Ray & _ray) noexcept :
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
        EIAPI FastFrustum(const Frustum& _frustum) noexcept :                          // TESTED
             vertices{} // Necessary, because const Vec3 is expected to remain uninitialized.
        {
            // Initialization of planes is difficult in the list, so the const-cast
            // only defers the initialization a bit.

            // Compute all 8 vertices (first get more help vectors).
            Vec3* v = const_cast<Vec3*>(vertices);
            Vec3 far = _frustum.f * _frustum.direction + _frustum.apex;
            Vec3 near = _frustum.n * _frustum.direction + _frustum.apex;
            // Get third axis and central off point
            Vec3 xAxis = cross(_frustum.up, _frustum.direction);
            Vec3 bottom = _frustum.b*_frustum.up;
            Vec3 top = _frustum.t*_frustum.up;
            Vec3 left =  _frustum.l * xAxis;
            Vec3 right =  _frustum.r * xAxis;
            float fton = _frustum.n / _frustum.f;
            v[0] = near + (left  + bottom) * fton;
            v[1] = near + (left  + top) * fton;
            v[2] = near + (right + bottom) * fton;
            v[3] = near + (right + top) * fton;
            v[4] = far  + left  + bottom;
            v[5] = far  + left  + top;
            v[6] = far  + right + bottom;
            v[7] = far  + right + top;

            // Create planes
            const_cast<DOP&>(nf) = DOP(_frustum.direction, near, far);
            // Use two vectors in the planes to derive the normal.
            const_cast<Plane&>(l) = Plane(normalize(cross(_frustum.up, v[4]-_frustum.apex)), v[4]);
            const_cast<Plane&>(r) = Plane(normalize(cross(v[7]-_frustum.apex, _frustum.up)), v[7]);
            const_cast<Plane&>(b) = Plane(normalize(cross(v[4]-_frustum.apex, xAxis)), v[4]);
            const_cast<Plane&>(t) = Plane(normalize(cross(xAxis, v[7]-_frustum.apex)), v[7]);
        }

        /// \brief Create from standard frustum parametrization
        EIAPI FastFrustum(const Vec3& _apex, const Vec3& _direction, const Vec3& _up, float _l, float _r, float _b, float _t, float _n, float _f) noexcept :
            FastFrustum(Frustum(_apex, _direction, _up, _l, _r, _b, _t, _n, _f))
        {}

        /// \brief Overwrite the current data (auto generation not possible because of const members)
        EIAPI FastFrustum& operator = (const FastFrustum& _frustum) noexcept
        {
            const_cast<DOP&>(nf) = _frustum.nf;
            const_cast<Plane&>(l) = _frustum.l;
            const_cast<Plane&>(r) = _frustum.r;
            const_cast<Plane&>(b) = _frustum.b;
            const_cast<Plane&>(t) = _frustum.t;
            for(int i = 0; i < 8; ++i) const_cast<Vec3&>(vertices[i]) = _frustum.vertices[i];
            return *this;
        }
    };

    struct FastCone
    {
        const Ray centralRay;
        const float cosThetaSq;
        const float height;       ///< Distance from origin to the base.

        /// \brief Construction from dynamic variant.
        EIAPI FastCone(const Cone & _cone) noexcept :
            centralRay(_cone.centralRay),
            cosThetaSq(1.0f / (1.0f + _cone.tanTheta * _cone.tanTheta)), // cos(atan(x))^2 == 1/(x^2+1)
            height(_cone.height)
        {}

        /// \brief Overwrite the current data (auto generation not possible because of const members)
        EIAPI FastCone& operator = (const FastCone& _cone) noexcept
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
        EIAPI explicit FastTriangle(const Triangle & _triangle) noexcept :
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
        EIAPI FastTriangle& operator = (const FastTriangle& _triangle) noexcept
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
    EIAPI inline float volume(const Sphere& _sphere) noexcept                        // TESTED
    {
        return 4.0f / 3.0f * PI * _sphere.radius * _sphere.radius * _sphere.radius;
    }

    EIAPI inline float volume(const Box& _box) noexcept                              // TESTED
    {
        Vec3 size = _box.max - _box.min;
        return size.x * size.y * size.z;
    }

    EIAPI inline float volume(const OBox& _obox) noexcept                            // TESTED
    {
        return 8.0f * _obox.halfSides.x * _obox.halfSides.y * _obox.halfSides.z;
    }

    EIAPI inline float volume(const Tetrahedron& _thetrahedron) noexcept             // TESTED
    {
        return dot(_thetrahedron.v3-_thetrahedron.v0, cross(_thetrahedron.v2-_thetrahedron.v0, _thetrahedron.v1-_thetrahedron.v0)) / 6.0f;
    }

    EIAPI inline float volume(const Triangle&) noexcept                              // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float volume(const Disc&) noexcept                                  // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float volume(const Plane&) noexcept                                 // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float volume(const DOP&) noexcept                                   // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float volume(const Ellipsoid& _ellipsoid) noexcept                  // TESTED
    {
        return 4.0f / 3.0f * PI * _ellipsoid.radii.x * _ellipsoid.radii.y * _ellipsoid.radii.z;
    }

    EIAPI inline float volume(const OEllipsoid& _oellipsoid) noexcept                // TESTED
    {
        return 4.0f / 3.0f * PI * _oellipsoid.radii.x * _oellipsoid.radii.y * _oellipsoid.radii.z;
    }

    EIAPI inline float volume(const Ray&) noexcept                                   // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float volume(const Segment&) noexcept                               // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float volume(const Cone& _cone) noexcept                            // TESTED
    {
        float r = _cone.height * _cone.tanTheta;
        return PI / 3.0f * r * r * _cone.height;
    }

    EIAPI inline float volume(const Capsule& _capsule) noexcept                      // TESTED
    {
        return PI * sq(_capsule.radius) * (_capsule.radius*4.0f/3.0f + len(_capsule.seg.b-_capsule.seg.a));
    }

    EIAPI inline float volume(const Frustum& _frustum) noexcept                      // TESTED
    {
        eiAssert( _frustum.r > _frustum.l, "Bad parametrization!" );
        eiAssert( _frustum.t > _frustum.b, "Bad parametrization!" );
        eiAssert( _frustum.f > _frustum.n, "Bad parametrization!" );
        // Use formular: V = (f-n)/3 * (a1 + sqrt(a1*a2) + a2) where a1 is the area on the
        // far plane and a2 that on the near plane.
        // Since a2 = sq(n/f) * a1 the formular simplifies to:
        // (f-n)/3 * (a1 + a1 * (n/f) + a1 * sq(n/f))
        float a1 = (_frustum.r - _frustum.l) * (_frustum.t - _frustum.b);
        float ratio = (_frustum.n / _frustum.f);
        return (_frustum.f - _frustum.n) / 3.0f * (a1 + (a1 + a1 * ratio) * ratio);
    }

    /// \brief Get the surface area of any object.
    EIAPI inline float surface(const Sphere& _sphere) noexcept                       // TESTED
    {
        return 4.0f * PI * sq(_sphere.radius);
    }

    EIAPI inline float surface(const Box& _box) noexcept                             // TESTED
    {
        Vec3 size = _box.max - _box.min;
        return 2.0f * (size.x * size.y + size.x * size.z + size.y * size.z);
    }

    EIAPI inline float surface(const OBox& _obox) noexcept                           // TESTED
    {
        return 8.0f * (_obox.halfSides.x * _obox.halfSides.y + _obox.halfSides.x * _obox.halfSides.z + _obox.halfSides.y * _obox.halfSides.z);
    }

    EIAPI inline float surface(const Tetrahedron& _thetrahedron) noexcept            // TESTED
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

    EIAPI inline float surface(const Triangle& _triangle) noexcept                   // TESTED
    {
        // Heron's formula is much more expensive than cross product because
        // the 3 side lengths must be computed first.
        return len( cross(_triangle.v1 - _triangle.v0, _triangle.v2 - _triangle.v0) ) * 0.5f;
    }

    EIAPI inline float surface(const Disc& _disc) noexcept                           // TESTED
    {
        return PI * _disc.radius;
    }

    EIAPI inline float surface(const Plane&) noexcept                                // TESTED
    {
        return INF;
    }

    EIAPI inline float surface(const DOP&) noexcept                                  // TESTED
    {
        return INF;
    }

    EIAPI inline float surface(const Ellipsoid& _ellipsoid) noexcept                 // TESTED
    {
        // Use approximation (Knud Thomsen's formula) only! Everything else is a
        // lot larger.
        Vec3 pr( pow(_ellipsoid.radii.x, 1.6075f),
            pow(_ellipsoid.radii.y, 1.6075f),
            pow(_ellipsoid.radii.z, 1.6075f) );
        return 4.0f * PI * pow((pr.x * pr.y + pr.x * pr.z + pr.y * pr.z) / 3.0f, 0.622083981f);
    }

    EIAPI inline float surface(const OEllipsoid& _oellipsoid) noexcept               // TESTED
    {
        // Use approximation (Knud Thomsen's formula) only! Everything else is a
        // lot larger.
        Vec3 pr( pow(_oellipsoid.radii.x, 1.6075f),
            pow(_oellipsoid.radii.y, 1.6075f),
            pow(_oellipsoid.radii.z, 1.6075f) );
        return 4.0f * PI * pow((pr.x * pr.y + pr.x * pr.z + pr.y * pr.z) / 3.0f, 0.622083981f);
    }

    EIAPI inline float surface(const Ray&) noexcept                                  // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float surface(const Segment&) noexcept                              // TESTED
    {
        return 0.0f;
    }

    EIAPI inline float surface(const Cone& _cone) noexcept                           // TESTED
    {
        float r = _cone.height * _cone.tanTheta;
        return PI * r * (r + sqrt(r*r + _cone.height*_cone.height));
    }

    EIAPI inline float surface(const Capsule& _capsule) noexcept                     // TESTED
    {
        return 2 * PI * _capsule.radius * (2 * _capsule.radius + len(_capsule.seg.b-_capsule.seg.a));
    }

    EIAPI inline float surface(const Frustum& _frustum) noexcept                     // TESTED
    {
        // Calculate area of all 6 sides and add them
        float hratio = (_frustum.f - _frustum.n) / _frustum.f;
        float ratio = _frustum.n / _frustum.f;
        // Far side
        float a1 = (_frustum.r - _frustum.l) * (_frustum.t - _frustum.b);
        // Near side
        float a2 = sq(ratio) * a1;
        // Top and bottom side
        float a3_a5 = (sqrt(sq(_frustum.f)+sq(_frustum.b)) + sqrt(sq(_frustum.f)+sq(_frustum.t))) * hratio * 0.5f * ((_frustum.r - _frustum.l) * (1.0f + ratio));
        // Left and right side
        float a4_a6 = (sqrt(sq(_frustum.f)+sq(_frustum.l)) + sqrt(sq(_frustum.f)+sq(_frustum.r))) * hratio * 0.5f * ((_frustum.t - _frustum.b) * (1.0f + ratio));
        return a1 + a2 + a3_a5 + a4_a6;
    }

    /// \brief Transform a box (rotation).
    EIAPI inline OBox transform(const Box& _box, const Quaternion& _rotation) noexcept
    {
        OBox box(_box);
        box.orientation *= _rotation;
        box.center = transform(box.center, _rotation);
        return box;
    }
    EIAPI inline OBox transform(const OBox& _box, const Quaternion& _rotation) noexcept
    {
        return OBox(transform(_box.center, _rotation), _box.halfSides, _box.orientation * _rotation);
    }
    /// \brief Transform a box (translation).
    EIAPI inline Box transform(const Box& _box, const Vec3& _translation) noexcept
    {
        return Box(_box.min + _translation, _box.max + _translation);
    }
    EIAPI inline OBox transform(const OBox& _box, const Vec3& _translation) noexcept
    {
        return OBox(_box.center + _translation, _box.halfSides, _box.orientation);
    }
    /// \brief Transform a box (first rotate then translate).
    EIAPI inline OBox transform(const Box& _box, const Quaternion& _rotation, const Vec3& _translation) noexcept
    {
        OBox box(_box);
        box.orientation *= _rotation;
        box.center = transform(box.center, _rotation) + _translation;
        return box;
    }
    EIAPI inline OBox transform(const OBox& _box, const Quaternion& _rotation, const Vec3& _translation) noexcept
    {
        return OBox(transform(_box.center, _rotation) + _translation,
            _box.halfSides, _box.orientation * _rotation);
    }

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
    EIAPI inline Vec3 center(const Sphere& _sphere) noexcept                         // TESTED
    {
        return _sphere.center;
    }

    EIAPI inline Vec3 center(const Box& _box) noexcept                               // TESTED
    {
        return (_box.min + _box.max) * 0.5f;
    }

    EIAPI inline Vec3 center(const OBox& _obox) noexcept                             // TESTED
    {
        return _obox.center;
    }

    EIAPI inline Vec3 center(const Tetrahedron& _thetrahedron) noexcept              // TESTED
    {
        return (_thetrahedron.v0 + _thetrahedron.v1 + _thetrahedron.v2 + _thetrahedron.v3) / 4.0f;
    }

    EIAPI inline Vec3 center(const Triangle& _triangle) noexcept                     // TESTED
    {
        return (_triangle.v0 + _triangle.v1 + _triangle.v2) / 3.0f;
    }

    EIAPI inline Vec3 center(const Disc& _disc) noexcept                             // TESTED
    {
        return _disc.center;
    }

    EIAPI inline Vec3 center(const Ellipsoid& _ellipsoid) noexcept                   // TESTED
    {
        return _ellipsoid.center;
    }

    EIAPI inline Vec3 center(const Segment& _line) noexcept                          // TESTED
    {
        return (_line.a + _line.b) * 0.5f;
    }

    EIAPI inline Vec3 center(const Cone& _cone) noexcept                             // TESTED
    {
        return _cone.centralRay.origin + _cone.centralRay.direction * 0.75f;
    }

    EIAPI inline Vec3 center(const Capsule& _capsule) noexcept                       // TESTED
    {
        return (_capsule.seg.a + _capsule.seg.b) * 0.5f;
    }

    EIAPI inline Vec3 center(const Frustum& _frustum) noexcept                       // TESTED
    {
        // For a full pyramid the centroid is located on the line segment from apex
        // to base centroid at 3/4 of that segment.
        Vec3 baseCentroidDir = (_frustum.l + _frustum.r) * 0.5f * cross(_frustum.up, _frustum.direction)
            + (_frustum.t + _frustum.b) * 0.5f * _frustum.up
            + _frustum.f * _frustum.direction;
        // Since this is a trunctad pyramid the height must be determined:
        // a1 = (r-l) * (t-b)    ground side area
        // a2 = a1 * sq(n/f)     intercept theorem
        // Relative distance of centroid to base centroid
        // d = (a1 + 2*sqrt(a1*a2) + 3*a2) / 4 * (a1 + sqrt(a1*a2) + a2)  http://mathworld.wolfram.com/PyramidalFrustum.html
        //   = (3*n*n + 2*f*n + f*f) / (n*n + f*n + f*f)
        // Distance from apex da = 1-d
        float da = 1 - (3*sq(_frustum.n) + 2*_frustum.n*_frustum.f + sq(_frustum.f))
                       / (4 * (sq(_frustum.n) + _frustum.n*_frustum.f + sq(_frustum.f)));
        // baseCentroidDir is the vector from apex to base and ´da´ the relative
        // distance from near to centroid plane. We need relative distance from
        // apex to centroid, which is ´(n+da*(f-n))/f´;
        return _frustum.apex + ((_frustum.n + da*(_frustum.f-_frustum.n))/_frustum.f) * baseCentroidDir;
    }

    namespace details {
        // Helper struct for faces on the convex hull.
        struct CHFace
        {
            UVec3 indices;
            UVec3 neighbors; // triangle neighbor indices for the edges (0,1), (1,2), (2,0).
            Plane p;
            bool isDeleted;
            CHFace() = default;
            CHFace(uint32 _i0, uint32 _i1, uint32 _i2, const UVec3& _neighbors, const Vec3* _points/*, uint32 _insidePoint*/) :
                indices(_i0, _i1, _i2),
                neighbors(_neighbors),
                p(_points[_i0], _points[_i1], _points[_i2]),
                isDeleted(false)
            {
            }
        };

        inline uint32 quickHull2D(Vec3* _points, uint32& _numPoints, const Vec3& _normal, Vec3 _a, Vec3 _b)
        {
            if(_numPoints == 0) return 0;
            // Find the furthest point on the right side relative to line a-b
            uint32 idx = _numPoints;
            Vec3 e0 = _b - _a;
            Vec3 edgeNormal = cross(e0, _normal);
            float dmax = 0.0f;
            for(uint32 i = 0; i < _numPoints; ++i)
            {
                float d = dot(_points[i] - _a, edgeNormal);
                if(d > dmax) { dmax = d; idx = i; }
                // TODO: list could be sorted after sign. Only positive side must be investigated by recursive calls.
            }

            if(idx != _numPoints)
            {
                if(idx != 0) std::swap(_points[0], _points[idx]);
                // Discard all points within the triangle a,b,p[i].
                Vec3 e1 = _points[0] - _a;
                float d00 = dot(e0, e0);
                float d01 = dot(e0, e1);
                float d11 = dot(e1, e1);
                float denom = d00 * d11 - d01 * d01;
                for(uint32 i = 1; i < _numPoints; ++i)
                {
                    float d20 = dot(_points[i] - _a, e0);
                    float d21 = dot(_points[i] - _a, e1);
                    float bary0 = (d11 * d20 - d01 * d21) / denom;
                    float bary1 = (d00 * d21 - d01 * d20) / denom;
                    if((bary0 >= 0.0f) && (bary1 >= 0.0f) && (bary0 + bary1 <= 1.0f))
                        std::swap(_points[i--], _points[--_numPoints]);
                }

                // Call recursive for the two new edges
                uint32 nconvex = 1;
                eiAssert(_numPoints > 0, "At least the new extreme point must exist!");
                _numPoints--; // Decrement because offsetting requires smaller bounds
                nconvex += quickHull2D(_points + nconvex, _numPoints, _normal, _a, _points[idx]);
                nconvex += quickHull2D(_points + nconvex, _numPoints, _normal, _points[idx], _b);
                return nconvex;
            }
            return 0;
        }
    } // namespace details

    float distance(const Vec3& _point, const Segment& _line); // Forward declaration

    /// \brief Remove all points from the array, which are not part of the
    ///     convex hull.
    /// \details The algorithm moves the points on the convex hull to the front
    ///     of the _points array. The other points are overwritten.
    /// \param [in] _threshold Discard vertices which are closer to the
    ///     previous convex hull than this threshold. This also includes
    ///     duplicates of vertices.
    /// \return Number of points in the convex set (these are the first elements
    ///     in _points after call).
    inline uint32 convexSet(Vec3* _points, uint32 _numPoints, float _threshold = 0.0f)
    {
        float tSq = _threshold * _threshold; // Compare squared numbers

        // Find duplicate points (brute force)
        for(uint32 i = 0; i < _numPoints-1; ++i)
        {
            for(uint32 j = i+1; j < _numPoints; ++j)
                if(lensq(_points[j] - _points[i]) <= tSq)
                {
                    // Whoops, j is the same as i -> remove j
                    _points[j--] = _points[--_numPoints];
                }
        }

        if(_numPoints <= 2) return _numPoints;

        // Find some extremal points for initial triangle.
        // First get two out of the six which define the AABox.
        uint32 extrema[6] = {0};
        for(uint32 i = 1; i < _numPoints; ++i)
        {
            if(_points[i].x < _points[extrema[0]].x) extrema[0] = i;
            if(_points[i].x > _points[extrema[1]].x) extrema[1] = i;
            if(_points[i].y < _points[extrema[2]].y) extrema[2] = i;
            if(_points[i].y > _points[extrema[3]].y) extrema[3] = i;
            if(_points[i].z < _points[extrema[4]].z) extrema[4] = i;
            if(_points[i].z > _points[extrema[5]].z) extrema[5] = i;
        }
        uint32 idx[2];
        float dmax = 0.0f;
        for(uint i = 0; i < 5; ++i)
        {
            for(uint j = i+1; j < 6; ++j)
            {
                float d = lensq(_points[extrema[i]] - _points[extrema[j]]);
                if(d > dmax)
                {
                    dmax = d;
                    idx[0] = extrema[i];
                    idx[1] = extrema[j];
                }
            }
        }
        eiAssert(dmax > tSq, "Points should have been deleted in duplicate search.");
        // Move the two points to the front. Order could be an issue if idx[1] == 0
        // or idx[0] == 1. If there is a cross reference swap that element first.
        if(idx[0] == 1) { std::swap(_points[0], _points[1]); std::swap(_points[1], _points[idx[1]]); }
        else { std::swap(_points[1], _points[idx[1]]); std::swap(_points[0], _points[idx[0]]); }
        // Get third point as most distant one to the line of the first two.
        dmax = 0.0f;
        Segment tsegment(_points[0], _points[1]);
        for(uint32 i = 2; i < _numPoints; ++i)
        {
            float d = distance(_points[i], tsegment);
            if(d > dmax)
            {
                dmax = d;
                idx[0] = i;
            }
        }
        if(dmax <= _threshold)
        {
            // All other points are convex combinations of the two extrema
            // (colinear point set).
            return 2;
        }
        std::swap(_points[2], _points[idx[0]]);
        if(_numPoints == 3) return 3;

        // Precondition: we have 3 non-colinear points in the front of the list.
        uint32 nconvex = 3;

        // All coplanar?
        Plane groundPlane(_points[0], _points[1], _points[2]);
        bool coplanar = true;
        for(uint32 i = nconvex; i < _numPoints && coplanar; ++i)
        {
            float d = dot(groundPlane.n, _points[i]) + groundPlane.d;
            if( abs(d) > _threshold ) coplanar = false;
        }
        if(coplanar)
        {
            // Perform an 2D variant of the quickhull algorithm.
            UVec2* edges = new UVec2[_numPoints];
            // Start with a loop of three edges
            _numPoints -= nconvex; // Decrement because offsetting requires smaller bounds
            nconvex += details::quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[0], _points[1]);
            nconvex += details::quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[1], _points[2]);
            nconvex += details::quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[2], _points[0]);

            delete[] edges;
            return nconvex;
        }

        // Setup memory for a temporary mesh
        details::CHFace* buf = new details::CHFace[_numPoints * 4]; // TODO: sufficient space?
        uint qs = 0, qe = 0, qn = _numPoints * 2;
        uint* queue = new uint[qn];
        uint* visibleStack = new uint[_numPoints * 2 + 2];
        UVec3* edges = new UVec3[_numPoints];

        // Put two opposite faces into our mesh and into a queue of
        // unvisited half spaces.
        buf[0] = details::CHFace(0, 1, 2, UVec3(1), _points);
        buf[1] = details::CHFace(0, 2, 1, UVec3(0), _points);
        uint nf = 2; // number of faces
        queue[qe] = 0; qe = (qe + 1) % qn;
        queue[qe] = 1; qe = (qe + 1) % qn;

        while(qs != qe && nconvex < _numPoints)
        {
            // Take a face from the queue. It can happen that deleted faces are listed.
            uint fidx = queue[qs]; qs = (qs + 1) % qn;
            while(buf[fidx].isDeleted && (qs != qe)) { fidx = queue[qs]; qs = (qs + 1) % qn; }
            if(buf[fidx].isDeleted && (qs == qe)) break;

            // Find the maximum distance point on the positive side.
            dmax = 0.0f;
            for(uint32 i = nconvex; i < _numPoints; ++i)
            {
                float d = dot(buf[fidx].p.n, _points[i]) + buf[fidx].p.d;
                if( d > dmax ) { dmax = d; idx[0] = i; }
            }
            if(dmax <= _threshold) continue;
            if(nconvex != idx[0]) std::swap(_points[nconvex], _points[idx[0]]);

            // For the new point find the set of visible faces.
            visibleStack[0] = fidx;
            int sp = 1;
            int ep = 0;
            while(sp > 0)
            {
                uint vidx = visibleStack[--sp];
                while(buf[vidx].isDeleted && sp > 0) vidx = visibleStack[--sp];
                if(buf[vidx].isDeleted && (sp <= 0)) break;

                // The point and the current visible triangle span a tetrahedron.
                // Remove all points inside the tetrahedron.
                Plane sides[4];
                sides[0].n = -buf[vidx].p.n; sides[0].d = -buf[vidx].p.d;
                sides[1] = Plane(_points[buf[vidx].indices.x], _points[buf[vidx].indices.y], _points[nconvex]);
                sides[2] = Plane(_points[buf[vidx].indices.y], _points[buf[vidx].indices.z], _points[nconvex]);
                sides[3] = Plane(_points[buf[vidx].indices.z], _points[buf[vidx].indices.x], _points[nconvex]);
                for(uint32 i = nconvex+1; i < _numPoints; ++i)
                {
                    bool inside = true;
                    for(int j = 0; j < 4 && inside; ++j)
                        inside &= (dot(sides[j].n, _points[i]) + sides[j].d) <= 0.0f;
                    if(inside)
                        std::swap(_points[i--], _points[--_numPoints]);
                }

                // Remove the visible face. It is now interior.
                buf[vidx].isDeleted = true;

                // Check the (up to) three neighbors.
                for(int i = 0; i < 3; ++i)
                {
                    uint32 nidx = buf[vidx].neighbors[i];
                    bool seesNeighbor = (dot(buf[nidx].p.n, _points[nconvex]) + buf[nidx].p.d) > 0.0f;

                    if(seesNeighbor)
                    {
                        visibleStack[sp++] = nidx;
                    } else {
                        // If there is an edge between visible and invisible face it will
                        // span a new face with the current point. Add to horizon list.
                        edges[ep++] = UVec3(buf[vidx].indices[i], buf[vidx].indices[(i+1)%3], nidx);
                    }
                }
            }

            // Sort all edges in the horizon list such that they give a loop.
            for(int i = 0; i < ep-1; ++i)
            {
                for(int j = i+1; j < ep; ++j)
                {
                    if(edges[i].y == edges[j].x)
                    {
                        if(i+1 != j) std::swap(edges[i+1], edges[j]);
                        break;
                    }
                }
            }

            // Create new faces to the horizon loop. To define the adjacencies we
            // need the sorted order.
            for(int e = 0; e < ep; ++e)
            {
                buf[nf + e].indices = UVec3(edges[e].x, edges[e].y, nconvex);
                buf[nf + e].isDeleted = false;
                buf[nf + e].neighbors = UVec3(edges[e].z, nf + (e + 1) % ep, nf + (e + ep - 1) % ep);
                buf[nf + e].p = Plane(_points[edges[e].x], _points[edges[e].y], _points[nconvex]);
                queue[qe] = nf + e; qe = (qe + 1) % qn;
            }
            nf += ep;

            nconvex++;
        }

        delete[] buf;
        delete[] queue;
        delete[] visibleStack;
        delete[] edges;
        return nconvex;
    }

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


    inline Box::Box( const OBox& _box ) noexcept
    {
        // Effectively generate all 8 corners and find min/max coordinates.
        // Relative to the center two diagonal opposite corners only differ
        // in the sign (even after rotation).

        Mat3x3 rot(_box.orientation);
        // Rows of the matrix are the aabox face directions in obox-space.
        // Choose diagonal entries via sign to get the largest coordinate into
        // face direction. Then project this maxCoord onto the direction.
        //Vec3 maxCoord = sgn(invRot(0)) * _box.halfSides;
        //max.x = dot(maxCoord, invRot(0));
        //max.x = dot(_box.halfSides, abs(rot(0))); // Equivalent to the sign stuff
        max = abs(rot) * _box.halfSides;
        min = -max;

        min += _box.center;
        max += _box.center;
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
