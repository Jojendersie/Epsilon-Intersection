#pragma once

#include "2dfunctions.hpp"

namespace ei {

    /// \brief Does a point lie inside a circle/on the boundary?
    /// \details Performance index: TODO
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec2& _point, const Disc2D& _disc );                // TESTED
    inline bool intersects( const Disc2D& _disc, const Vec2& _point )  { return intersects( _point, _disc ); }

    /// \brief Do two discs intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the discs intersect in two or one point or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _disc0, const Disc2D& _disc1 );             // TESTED

    /// \brief Do two discs intersect/touch?
    /// \details Performance index: TODO
    /// \param [out] _outInfo The central point of the intersection area and
    ///    undefined if no intersection occurred. If one circle is contained
    ///    in the other one this info is its center.
    /// \return true if the discs intersect in two or one point or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _disc0, const Disc2D& _disc1, Vec2& _outInfo ); // TESTED

    /// \brief Do two rectangles intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the rectangles intersect or one lies
    ///    fully inside the other one.
    bool intersects( const Rect2D& _rect0, const Rect2D& _rect1 );             // TESTED

    /// \brief Does a point lie inside a rectangle/on the boundary?
    /// \details Performance index: TODO
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec2& _point, const Rect2D& _rect );                // TESTED
    inline bool intersects( const Rect2D& _rect, const Vec2& _point )  { return intersects( _point, _rect ); }

    /// \brief Do a rectangle and disc intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the rectangle and disc intersect or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _disc, const Rect2D& _rect );               // TESTED
    inline bool intersects( const Rect2D& _rect, const Disc2D& _disc )  { return intersects(_disc, _rect); }

    /// \brief Do two lines intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the lines intersect.
    bool intersects( const Line2D& _line0, const Line2D& _line1 );             // TESTED

    /// \brief Where do two lines intersect/touch?
    /// \details Performance index: TODO
    /// \param [out] _outInfo The point of the intersection. If the two lines
    ///    overlap the returned point is the center of this overlap.
    /// \return true if the lines intersect.
    bool intersects( const Line2D& _line0, const Line2D& _line1, Vec2& _outInfo ); // TESTED
}