#pragma once

#include "3dfunctions.hpp"

namespace ei {

    /// \brief Does a point lie inside an ellipsoid/on the boundary?
    /// \details Performance index: 5.68
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid );        // TESTED
    inline bool intersects( const Ellipsoid& _ellipsoid, const Vec3& _point )  { return intersects( _point, _ellipsoid ); }

    /// \brief Do a ray and an ellipsoid intersect or touch?
    /// \details Performance index: 23.1, 23.9 (Numbers for increasing
    ///     usage of optional parameters)
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \return true if there is at least one comon point between ray and
    ///     ellipsoid.
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid );           // TESTED
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid, float& _distance );
    inline bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray )  { return intersects( _ray, _ellipsoid ); }
    inline bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray, float& _distance )  { return intersects( _ray, _ellipsoid, _distance ); }

    /// \brief Do a ray and a box intersect or touch?
    /// \details Performance index: 17.1, 18.0 (Numbers for increasing
    ///     usage of optional parameters)
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \return true if there is at least one comon point between ray and box
    bool intersects( const Ray& _ray, const Box& _box );                       // TESTED
    bool intersects( const Ray& _ray, const Box& _box, float& _distance );
    inline bool intersects( const Box& _box, const Ray& _ray )  { return intersects( _ray, _box ); }
    inline bool intersects( const Box& _box, const Ray& _ray, float& _distance )  { return intersects( _ray, _box, _distance ); }

    /// \brief Do a ray and a triangle intersect or touch?
    /// \details Performance index: 12.3, 13.7, TODO (Numbers for increasing
    ///     usage of optional parameters)
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    /// \param [out,opt] _barycentric The barycentric coordinates of the hit point
    ///     on the triangle.
    /// \return true if there is at least one comon point between ray and triangle.
    ///     This point has a ray parameter >= 0 (no negative direction).
    bool intersects( const Ray& _ray, const Triangle& _triangle );                                          // TESTED
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance );                        // TESTED
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance, Vec3& _barycentric );    // TESTED
    inline bool intersects( const Triangle& _triangle, const Ray& _ray )  { return intersects( _ray, _triangle ); }
    inline bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance )  { return intersects( _ray, _triangle, _distance ); }
    inline bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance, Vec3& _barycentric )  { return intersects( _ray, _triangle, _distance, _barycentric ); }
}
