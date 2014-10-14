#pragma once

#include "2dtypes.hpp"

namespace ei {

    /// \brief Does a point lies inside a circle/on the boundary?
    /// \details Performance index: TODO
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec2& _point, const Disc2D& _disc );                // TESTED
    inline bool intersects( const Disc2D& _disc, const Vec2& _point )  { return intersects( _point, _disc ); }

    /// \brief Does two discs intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the discs intersect in two or one point or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _disc0, const Disc2D& _disc1 );             // TESTED

    /// \brief Does two discs intersect/touch?
    /// \details Performance index: TODO
    /// \param [out] _outInfo The central point of the intersection area and
    ///    undefined if no intersection occurred. If one circle is contained
    ///    in the other one this info is its center.
    /// \return true if the discs intersect in two or one point or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _disc0, const Disc2D& _disc1, Vec2& _outInfo ); // TESTED

    /// \brief Does two rectangles intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the rectangles intersect or one lies
    ///    fully inside the other one.
    bool intersects( const Rect2D& _rect0, const Rect2D& _rect1 );             // TESTED

    /// \brief Does a rectangle and disc intersect/touch?
    /// \details Performance index: TODO
    /// \return true if the rectangle and disc intersect or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _disc, const Rect2D& _rect );               // TESTED
    inline bool intersects( const Rect2D& _rect, const Disc2D& _disc )  { return intersects(_disc, _rect); }
}