#pragma once

#include "matrix.hpp"

namespace ei {

    // Predeclarations for all types to enable conversion operators.
    struct Circle2D;
    struct Rect2D;
    struct ORect2D;
    struct Triangle2D;
    struct Ellipse2D;
    struct OEllipse2D;
    struct Capsule2D;
    struct Line2D;
    struct Ray2D;

    struct Circle2D
    {
        Vec2 center;
        float radius;

        /// \brief Create from center and radius
        Circle2D( Vec2 _center, float _radius );

        /// \brief Create a circle which encloses two points
        Circle2D( Vec2 _p0, Vec2 _p1 );

        /// \brief Create a circle which encloses three points
        Circle2D( Vec2 _p0, Vec2 _p1, Vec2 _p2 );

        /// \brief Create circumcircle of a triangle
        explicit Circle2D( const Triangle2D& _triangle );

        /// \brief Create circumcircle of a rect
        explicit Circle2D( const Rect2D& _rect );

        /// \brief Create circumcircle of a oriented rect
        explicit Circle2D( const ORect2D& _rect );

        /// \brief Create circumcircle of an ellipse
        explicit Circle2D( const Ellipse2D& _ellipse );

        /// \brief Create circumcircle of an oriented ellipse
        explicit Circle2D( const OEllipse2D& _ellipse );

        /// \brief Create circumcircle of a capsule
        explicit Circle2D( const Capsule2D& _capsule );

        /// \brief Create circumcircle of a line
        explicit Circle2D( const Line2D& _line );
    };

    struct Rect2D
    {
        Vec2 min;
        Vec2 max;
    };

    struct ORect2D
    {
        Vec2 center;            ///< Center of the oriented rect.
        Vec2 size;              ///< Width and height.
        float angle;            ///< Angle in radiant.
    };

    struct Triangle2D
    {
        Vec2 v0;
        Vec2 v1;
        Vec2 v2;
    };

    struct Ellipse2D
    {
        Vec2 center;
        Vec2 radius;
    };

    struct OEllipse2D
    {
        Vec2 center;
        Vec2 radius;
        float angle;
    };

    struct Line2D
    {
        Vec2 p0;
        Vec2 p1;
        float radius;
    };

    struct Ray2D
    {
        Vec2 origin;
        Vec2 direction;         ///< Normalized direction vector.
        float radius;
    };

    /// \brief A capsule is a region around a line.
    struct Capsule2D
    {
        Vec2 p0;
        Vec2 p1;
        float radius;           ///< Defines the region around the line.
    };
};