#pragma once

#include "vector.hpp"

namespace ei {

    // Predeclarations for all types to enable conversion operators.
    struct Disc2D;
    struct Rect2D;
    struct ORect2D;
    struct Triangle2D;
    struct Ellipse2D;
    struct OEllipse2D;
    struct Capsule2D;
    struct Segment2D;
    struct Ray2D;

    /// \brief A list of all supported 2d types
    enum struct Types2D
    {
        DISC,
        RECT,
        ORECT,
        TRIANGLE,
        ELLIPSE,
        OELLIPSE,
        CAPSULE,
        LINE,
        RAY,

        NUM_TYPES
    };

    struct Rect2D
    {
        Vec2 min;
        Vec2 max;

        /// \brief Create uninitialized rectangle.
        EIAPI Rect2D() noexcept {}

        /// \brief Construct from minimal and maximal coordinates
        EIAPI Rect2D( const Vec2& _min, const Vec2& _max ) noexcept :                // TESTED
            min(_min),
            max(_max)
        {
            eiAssert( _min <= _max,
                "Minimum coordinates must be smaller or equal the maximum." );
        }

        /// \brief Create an optimal box for a set of points
        EIAPI Rect2D( const Vec2* _points, uint32 _numPoints ) noexcept
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

    struct ORect2D
    {
        Vec2 center;            ///< Center of the oriented rect.
        Vec2 size;              ///< Width and height.
        float angle;            ///< Angle in radiant.

        /// \brief Create uninitialized oriented rectangle.
        EIAPI ORect2D() noexcept {}

        /// \brief Create from center, side lengths and rotation
        EIAPI ORect2D( const Vec2& _center, const Vec2& _size, float _angle ) noexcept :  // TESTED
            center(_center),
            size(_size),
            angle(_angle)
        {
            eiAssert( _size >= 0.0f, "Side lengths must be positive." );
        }
    };

    /// \brief A triangle with counter clock wise vertices.
    struct Triangle2D
    {
        Vec2 v0;
        Vec2 v1;
        Vec2 v2;

        /// \brief Create uninitialized triangle.
        EIAPI Triangle2D() noexcept {}

        /// \brief Construct from three vertices.
        EIAPI Triangle2D( const Vec2& _v0, const Vec2& _v1, const Vec2& _v2 ) noexcept :  // TESTED
            v0(_v0),
            v1(_v1),
            v2(_v2)
        {}
    };

    struct Ellipse2D
    {
        Vec2 center;
        Vec2 radii;

        /// \brief Create uninitialized ellipse.
        EIAPI Ellipse2D() noexcept {}

        /// \brief Create from center and radii.
        EIAPI Ellipse2D( const Vec2& _center, const Vec2& _radii ) noexcept :        // TESTED
            center(_center),
            radii(_radii)
        {
            eiAssert( _radii >= 0,
                "Radii of an ellipse should never be negative." );
        }
    };

    struct OEllipse2D
    {
        Vec2 center;
        Vec2 radii;
        float angle;

        /// \brief Create uninitialized oriented ellipse.
        EIAPI OEllipse2D() noexcept {}

        /// \brief Create from center, radii and angle.
        EIAPI OEllipse2D( const Vec2& _center, const Vec2& _radii, float _angle ) noexcept :  // TESTED
            center(_center),
            radii(_radii),
            angle(_angle)
        {
            eiAssert( _radii >= 0,
                "Radii of an ellipse should never be negative." );
        }
    };

    /// \brief A line segment is the finite connection between two points
    struct Segment2D
    {
        Vec2 a;
        Vec2 b;

        /// \brief Create uninitialized line.
        EIAPI Segment2D() noexcept {}

        /// \brief Create from two points (start and end).
        EIAPI Segment2D( const Vec2& _a, const Vec2& _b ) noexcept :                 // TESTED
            a(_a),
            b(_b)
        {}
    };

    struct Ray2D
    {
        Vec2 origin;
        Vec2 direction;         ///< Normalized direction vector.

        /// \brief Create uninitialized ray.
        EIAPI Ray2D() noexcept {}

        /// \brief Create from point (origin) and direction
        /// \param [in] _direction Normalized direction vector.
        EIAPI Ray2D( const Vec2& _origin, const Vec2& _direction ) noexcept :        // TESTED
            origin(_origin),
            direction(_direction)
        {
            eiAssert( approx(len(_direction), 1.0f),
                "The direction vector must be normalized." );
        }
    };

    /// \brief A capsule is a region around a line.
    struct Capsule2D
    {
        Segment2D seg;
        float radius;           ///< Defines the region around the line.

        /// \brief Create uninitialized capsule.
        EIAPI Capsule2D() noexcept {}

        /// \brief Create from two points and a radius for the surrounding
        ///    region.
        EIAPI Capsule2D( const Vec2& _a, const Vec2& _b, float _radius ) noexcept :   // TESTED
            seg(_a, _b),
            radius(_radius)
        {}
    };

    /// \brief Standard disc (solid filled circle)
    /// \details You may change the position and radius arbitrary.
    ///
    ///    The radius should never be negative.
    ///    A radius can be zero.
    struct Disc2D
    {
        Vec2 center;
        float radius;

        /// \brief Create uninitialized circle.
        EIAPI Disc2D() noexcept {}

        /// \brief Create from center and radius
        EIAPI Disc2D( Vec2 _center, float _radius ) noexcept :                       // TESTED
            center(_center),
            radius(_radius)
        {}

        /// \brief Create a circle which encloses two points
        EIAPI Disc2D( Vec2 _p0, Vec2 _p1 ) noexcept                                  // TESTED
        {
            center = (_p0 + _p1) * 0.5f;
            radius = len(_p0 - _p1) * 0.5f;
        }

        /// \brief Create a circle which encloses three points
        EIAPI Disc2D( Vec2 _p0, Vec2 _p1, Vec2 _p2 ) noexcept                        // TESTED
        {
            // The center of the circumscribed circle is at (barycentric coords)
            // v0*sin(2 alpha) + v1*sin(2 beta) + v2*sin(2 gamma) and has the radius
            // abc/4A.
            Vec2 c = _p0 - _p1;	float csq = lensq(c);
            Vec2 a = _p1 - _p2;	float asq = lensq(a);
            Vec2 b = _p2 - _p0;	float bsq = lensq(b);

            // One of the sides could be the longest side - the minimum sphere is
            // defined through only two points.
            // This can also handle the coplanar case.
            if( csq + bsq <= asq ) *this = Disc2D(_p1, _p2);
            else if( asq + bsq <= csq ) *this = Disc2D(_p1, _p0);
            else if( asq + csq <= bsq ) *this = Disc2D(_p2, _p0);
            else {
                float area2Sq = 2;// * lensq(cross(a, c));
                center = 
                    _p0 * (-dot(c,b) * asq / area2Sq)
                    + _p1 * (-dot(c,a) * bsq / area2Sq)
                    + _p2 * (-dot(b,a) * csq / area2Sq);
                radius = sqrt(asq * bsq * csq / (2 * area2Sq));
            }
        }

        /// \brief Create circumcircle of a rect
        EIAPI explicit Disc2D( const Rect2D& _rect ) noexcept                        // TESTED
        {
            eiAssert( _rect.max >= _rect.min,
                "The input rect is degenerated! All components of max should be larger than in min." );
            // For a regular n-gone the center is simply the center of all vertices.
            // The radius is then equal to all vertices.
            center = (_rect.max + _rect.min) * 0.5f;
            radius = len(_rect.max - center);
        }

        /// \brief Create circumcircle of a oriented rect
        EIAPI explicit Disc2D( const ORect2D& _rect ) noexcept                       // TESTED
        {
            eiAssert( _rect.size >= 0,
                "Side lengths of a rectangle should never be negative." );
            // The center is already given and the distance to all corners is equal.
            center = _rect.center;
            radius = len(_rect.size) * 0.5f;
        }

        /// \brief Create circumcircle of a triangle
        EIAPI explicit Disc2D( const Triangle2D& _triangle ) noexcept :              // TESTED
            Disc2D( _triangle.v0, _triangle.v1, _triangle.v2 )
        {}

        /// \brief Create circumcircle of an ellipse
        EIAPI explicit Disc2D( const Ellipse2D& _ellipse ) noexcept                  // TESTED
        {
            eiAssert( _ellipse.radii >= 0,
                "Radii of an ellipse should never be negative." );
            center = _ellipse.center;
            radius = max( _ellipse.radii );
        }

        /// \brief Create circumcircle of an oriented ellipse
        EIAPI explicit Disc2D( const OEllipse2D& _ellipse ) noexcept                 // TESTED
        {
            eiAssert( _ellipse.radii >= 0,
                "Radii of an ellipse should never be negative." );
            // Rotation does not change anything
            center = _ellipse.center;
            radius = max( _ellipse.radii );
        }

        /// \brief Create circumcircle of a capsule
        EIAPI explicit Disc2D( const Capsule2D& _capsule ) noexcept                  // TESTED
        {
            center = (_capsule.seg.a + _capsule.seg.b) * 0.5f;
            radius = len(_capsule.seg.a - _capsule.seg.b) * 0.5f + _capsule.radius;
        }

        /// \brief Create circumcircle of a line
        EIAPI explicit Disc2D( const Segment2D& _line ) noexcept                     // TESTED
        {
            center = (_line.a + _line.b) * 0.5f;
            radius = len(_line.a - _line.b) * 0.5f;
        }

        /// \brief Compare on binary identity
        EIAPI bool operator== ( const Disc2D& _circle ) const noexcept               // TESTED
        {
            return ( center == _circle.center ) && radius == _circle.radius;
        }
    };

    // TODO: approx, area, center

    // ************************************************************************* //
    // AREA METHODS
    // ************************************************************************* //
    EIAPI float area(const Disc2D& _disc) noexcept                            // TESTED
    {
        return _disc.radius * _disc.radius * PI;
    }

    EIAPI float area(const Rect2D& _rect) noexcept                            // TESTED
    {
        eiAssert(_rect.max >= _rect.min, "Rect max-boundary must be larger than its min-boundary!");
        return (_rect.max.x - _rect.min.x) * (_rect.max.y - _rect.min.y);
    }

    EIAPI float area(const ORect2D& _orect) noexcept                          // TESTED
    {
        eiAssert(_orect.size >= 0.0f, "Rect must have positive side lengths!");
        return _orect.size.x * _orect.size.y;
    }

    EIAPI float area(const Triangle2D& _triangle) noexcept                    // TESTED
    {
        // Use determinant rule for the cross product of two sides
        return 0.5f * ((_triangle.v1.x-_triangle.v0.x) * (_triangle.v2.y-_triangle.v0.y)
                     - (_triangle.v2.x-_triangle.v0.x) * (_triangle.v1.y-_triangle.v0.y));
    }

    EIAPI float area(const Ellipse2D& _ellipse) noexcept                      // TESTED
    {
        return PI * _ellipse.radii.x * _ellipse.radii.y;
    }

    EIAPI float area(const OEllipse2D& _oellipse) noexcept                    // TESTED
    {
        return PI * _oellipse.radii.x * _oellipse.radii.y;
    }

    EIAPI float area(const Segment2D&) noexcept                               // TESTED
    {
        return 0.0f;
    }

    EIAPI float area(const Ray2D&) noexcept                                   // TESTED
    {
        return 0.0f;
    }

    EIAPI float area(const Capsule2D& _capsule) noexcept                      // TESTED
    {
        return PI * _capsule.radius * _capsule.radius + 2.0f * _capsule.radius * len(_capsule.seg.b-_capsule.seg.a);
    }

    // ************************************************************************* //
    // CENTROID METHODS
    // ************************************************************************* //
    EIAPI Vec2 center(const Disc2D& _disc) noexcept                           // TESTED
    {
        return _disc.center;
    }

    EIAPI Vec2 center(const Rect2D& _rect) noexcept                           // TESTED
    {
        return 0.5f * (_rect.min + _rect.max);
    }

    EIAPI Vec2 center(const ORect2D& _orect) noexcept                         // TESTED
    {
        return _orect.center;
    }

    EIAPI Vec2 center(const Triangle2D& _triangle) noexcept                   // TESTED
    {
        return (_triangle.v0 + _triangle.v1 + _triangle.v2) / 3.0f;
    }

    EIAPI Vec2 center(const Ellipse2D& _ellipse) noexcept                     // TESTED
    {
        return _ellipse.center;
    }

    EIAPI Vec2 center(const OEllipse2D& _oellipse) noexcept                   // TESTED
    {
        return _oellipse.center;
    }

    EIAPI Vec2 center(const Segment2D& _segment) noexcept                     // TESTED
    {
        return (_segment.a + _segment.b) * 0.5f;
    }

    EIAPI Vec2 center(const Capsule2D& _capsule) noexcept                     // TESTED
    {
        return (_capsule.seg.a + _capsule.seg.b) * 0.5f;
    }
}
