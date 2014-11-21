#pragma once

#include "3dfunctions.hpp"

namespace ei {

    /// \brief Does a point lie inside an ellipsoid/on the boundary?
    /// \details Performance index: 8.9
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid );        // TESTED
    inline bool intersects( const Ellipsoid& _ellipsoid, const Vec3& _point )  { return intersects( _point, _ellipsoid ); }

    /// \brief Do a ray and an ellipsoid intersect or touch?
    /// \details Performance index: 12.8
    /// \return true if there is at least one comon point between ray and ellipsoid
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid );           // TESTED
    inline bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray )  { return intersects( _ray, _ellipsoid ); }

    /// \brief Do a ray and a box intersect or touch?
    /// \details Performance index: 17.1
    /// \return true if there is at least one comon point between ray and box
    bool intersects( const Ray& _ray, const Box& _box );                       // TESTED
    inline bool intersects( const Box& _box, const Ray& _ray )  { return intersects( _ray, _box ); }
}
