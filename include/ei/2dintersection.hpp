#pragma once

#include "2dtypes.hpp"

namespace ei {

    /// \brief Does two circles intersect?
    /// \details Performance index: TODO
    /// \return true if the circles intersect in two or one point or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _circle0, const Disc2D& _circle1 );         // TESTED

    /// \brief Does two circles intersect?
    /// \details Performance index: TODO
    /// \param [out] _outInfo The central point of the intersection area and
    ///    undefined if no intersection occurred. If one circle is contained
    ///    in the other one this info is its center.
    /// \return true if the circles intersect in two or one point or one lies
    ///    fully inside the other one.
    bool intersects( const Disc2D& _circle0, const Disc2D& _circle1, Vec2& _outInfo );
}