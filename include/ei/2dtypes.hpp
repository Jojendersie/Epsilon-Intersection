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
        Disc2D() noexcept {}

        /// \brief Create from center and radius
        Disc2D( Vec2 _center, float _radius ) noexcept;                        // TESTED

        /// \brief Create a circle which encloses two points
        Disc2D( Vec2 _p0, Vec2 _p1 ) noexcept;                                 // TESTED

        /// \brief Create a circle which encloses three points
        Disc2D( Vec2 _p0, Vec2 _p1, Vec2 _p2 ) noexcept;                       // TESTED

        /// \brief Create circumcircle of a rect
        explicit Disc2D( const Rect2D& _rect ) noexcept;                       // TESTED

        /// \brief Create circumcircle of a oriented rect
        explicit Disc2D( const ORect2D& _rect ) noexcept;                      // TESTED

        /// \brief Create circumcircle of a triangle
        explicit Disc2D( const Triangle2D& _triangle ) noexcept;               // TESTED

        /// \brief Create circumcircle of an ellipse
        explicit Disc2D( const Ellipse2D& _ellipse ) noexcept;                 // TESTED

        /// \brief Create circumcircle of an oriented ellipse
        explicit Disc2D( const OEllipse2D& _ellipse ) noexcept;                // TESTED

        /// \brief Create circumcircle of a capsule
        explicit Disc2D( const Capsule2D& _capsule ) noexcept;                 // TESTED

        /// \brief Create circumcircle of a line
        explicit Disc2D( const Segment2D& _line ) noexcept;                    // TESTED

        /// \brief Compare on binary identity
        bool operator== ( const Disc2D& _circle ) const noexcept;              // TESTED
    };

    struct Rect2D
    {
        Vec2 min;
        Vec2 max;

        /// \brief Create uninitialized rectangle.
        Rect2D() noexcept {}

        /// \brief Construct from minimal and maximal coordinates
        Rect2D( const Vec2& _min, const Vec2& _max ) noexcept;                 // TESTED

        /// \brief Create an optimal box for a set of points
        Rect2D( const Vec2* _points, uint32 _numPoints ) noexcept;
    };

    struct ORect2D
    {
        Vec2 center;            ///< Center of the oriented rect.
        Vec2 size;              ///< Width and height.
        float angle;            ///< Angle in radiant.

        /// \brief Create uninitialized oriented rectangle.
        ORect2D() noexcept {}

        /// \brief Create from center, side lengths and rotation
        ORect2D( const Vec2& _center, const Vec2& _size, float _angle ) noexcept;   // TESTED
    };

    /// \brief A triangle with counter clock wise vertices.
    struct Triangle2D
    {
        Vec2 v0;
        Vec2 v1;
        Vec2 v2;

        /// \brief Create uninitialized triangle.
        Triangle2D() noexcept {}

        /// \brief Construct from three vertices.
        Triangle2D( const Vec2& _v0, const Vec2& _v1, const Vec2& _v2 ) noexcept;   // TESTED
    };

    struct Ellipse2D
    {
        Vec2 center;
        Vec2 radii;

        /// \brief Create uninitialized ellipse.
        Ellipse2D() noexcept {}

        /// \brief Create from center and radii.
        Ellipse2D( const Vec2& _center, const Vec2& _radii ) noexcept;         // TESTED
    };

    struct OEllipse2D
    {
        Vec2 center;
        Vec2 radii;
        float angle;

        /// \brief Create uninitialized oriented ellipse.
        OEllipse2D() noexcept {}

        /// \brief Create from center, radii and angle.
        OEllipse2D( const Vec2& _center, const Vec2& _radii, float _angle ) noexcept;   // TESTED
    };

    /// \brief A line segment is the finite connection between two points
    struct Segment2D
    {
        Vec2 a;
        Vec2 b;

        /// \brief Create uninitialized line.
        Segment2D() noexcept {}

        /// \brief Create from two points (start and end).
        Segment2D( const Vec2& _a, const Vec2& _b ) noexcept;                  // TESTED
    };

    struct Ray2D
    {
        Vec2 origin;
        Vec2 direction;         ///< Normalized direction vector.

        /// \brief Create uninitialized ray.
        Ray2D() noexcept {}

        /// \brief Create from point (origin) and direction
        /// \param [in] _direction Normalized direction vector.
        Ray2D( const Vec2& _origin, const Vec2& _direction ) noexcept;         // TESTED
    };

    /// \brief A capsule is a region around a line.
    struct Capsule2D
    {
        Segment2D seg;
        float radius;           ///< Defines the region around the line.

        /// \brief Create uninitialized capsule.
        Capsule2D() noexcept {}

        /// \brief Create from two points and a radius for the surrounding
        ///    region.
        Capsule2D( const Vec2& _a, const Vec2& _b, float _radius ) noexcept;   // TESTED
    };

    // TODO: approx, area, center
    float area(const Disc2D& _disc) noexcept;                                  // TESTED
    float area(const Rect2D& _rect) noexcept;                                  // TESTED
    float area(const ORect2D& _orect) noexcept;                                // TESTED
    float area(const Triangle2D& _triangle) noexcept;                          // TESTED
    float area(const Ellipse2D& _ellipse) noexcept;                            // TESTED
    float area(const OEllipse2D& _oellipse) noexcept;                          // TESTED
    float area(const Segment2D& _segment) noexcept;                            // TESTED
    float area(const Ray2D& _ray) noexcept;                                    // TESTED
    float area(const Capsule2D& _capsule) noexcept;                            // TESTED

    Vec2 center(const Disc2D& _disc) noexcept;                                 // TESTED
    Vec2 center(const Rect2D& _rect) noexcept;                                 // TESTED
    Vec2 center(const ORect2D& _orect) noexcept;                               // TESTED
    Vec2 center(const Triangle2D& _triangle) noexcept;                         // TESTED
    Vec2 center(const Ellipse2D& _ellipse) noexcept;                           // TESTED
    Vec2 center(const OEllipse2D& _oellipse) noexcept;                         // TESTED
    Vec2 center(const Segment2D& _segment) noexcept;                           // TESTED
    Vec2 center(const Capsule2D& _capsule) noexcept;                           // TESTED

    // Include implementation of short constructors.
#   include "details/2dtypes.inl"
}
