#pragma once

#include "3dtypes.hpp"

namespace ei {

    /// \brief Does a point lie inside an ellipsoid/on the boundary?
    /// \details Performance index: TODO
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid );        // TESTED
    inline bool intersects( const Ellipsoid& _ellipsoid, const Vec3& _point )  { return intersects( _point, _ellipsoid ); }

}
