#pragma once

#include "3dfunctions.hpp"

namespace ei {

    /// \brief Do two spheres intersect, touch or is one inside the other?
    /// \details Performance index: 3.95
    /// \return true if both spheres have at least one point in common.
    bool intersects( const Sphere& _sphere0, const Sphere& _sphere1 );         // TESTED

    /// \brief Does a point lie inside a sphere/on the boundary?
    /// \details Performance index: 3.33
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Sphere& _sphere );              // TESTED
    inline bool intersects( const Sphere& _sphere, const Vec3& _point )  { return intersects( _point, _sphere ); }

    /// \brief Does a sphere lie inside a box/touches the boundary?
    /// \details Performance index: 12.2
    /// \return true if the sphere and the box have at least one point in common.
    bool intersects( const Sphere& _sphere, const Box& _box );                 // TESTED
    inline bool intersects( const Box& _box, const Sphere& _sphere )  { return intersects( _sphere, _box ); }

    /// \brief Do two boxes intersect, touch or is one inside the other?
    /// \details Performance index: 16.7
    /// \return true if both boxes have at least one point in common.
    bool intersects( const Box& _box0, const Box& _box1 );                     // TESTED

    /// \brief Do two capsules intersect, touch or is one inside the other?
    /// \details Performance index: TODO
    /// \return true if both spheres have at least one point in common.
    bool intersects( const Capsule& _capsule0, const Capsule& _capsule1 );     // TESTED

    /// \brief Does a point lie inside an ellipsoid/on the boundary?
    /// \details Performance index: 6.9
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid );        // TESTED
    inline bool intersects( const Ellipsoid& _ellipsoid, const Vec3& _point )  { return intersects( _point, _ellipsoid ); }

    /// \brief Do a ray and an ellipsoid intersect or touch?
    /// \details Performance index: 22.7, 25.1 (Numbers for increasing
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
    /// \details Performance index: 27.6, 27.2 (Numbers for increasing
    ///     usage of optional parameters)
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \param [out,opt] _distanceExit The ray parameter (distance) for the
    ///     second intersection point in positive direction.
    ///
    ///     This is always the exit point. If the ray starts on the boundary
    ///     and shows away _distance and _distanceExit are the same (0).
    /// \return true if there is at least one comon point between ray and box
    bool intersects( const Ray& _ray, const Box& _box );                       // TESTED
    bool intersects( const Ray& _ray, const Box& _box, float& _distance );
    bool intersects( const Ray& _ray, const Box& _box, float& _distance, float& _distanceExit );
    inline bool intersects( const Box& _box, const Ray& _ray )  { return intersects( _ray, _box ); }
    inline bool intersects( const Box& _box, const Ray& _ray, float& _distance )  { return intersects( _ray, _box, _distance ); }
    inline bool intersects( const Box& _box, const Ray& _ray, float& _distance, float& _distanceExit )  { return intersects( _ray, _box, _distance, _distanceExit ); }

    /// \brief Do a ray and a triangle intersect or touch?
    /// \details Performance index: 11.9, 16.1, TODO (Numbers for increasing
    ///     usage of optional parameters)
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    /// \param [out,opt] _barycentric The barycentric coordinates of the hit point
    ///     on the triangle.
    /// \return true if there is at least one common point between ray and triangle.
    ///     This point has a ray parameter >= 0 (no negative direction).
    bool intersects( const Ray& _ray, const Triangle& _triangle );                                          // TESTED
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance );                        // TESTED
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance, Vec3& _barycentric );    // TESTED
    inline bool intersects( const Triangle& _triangle, const Ray& _ray )  { return intersects( _ray, _triangle ); }
    inline bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance )  { return intersects( _ray, _triangle, _distance ); }
    inline bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance, Vec3& _barycentric )  { return intersects( _ray, _triangle, _distance, _barycentric ); }

    /// \brief Does the sphere touches the triangle?
    /// \details Performance index: TODO
    /// \return true if the sphere and the triangle have at least on point in common.
    bool intersects( const Sphere& _sphere, const Triangle& _triangle );       // TESTED
    inline bool intersects( const Triangle& _triangle, const Sphere& _sphere )          { return intersects(_sphere, _triangle); }

    /// \brief Intersection test between sphere and capsule.
    /// \details Performance index: TODO
    /// \return true if the sphere and the capsule have at least on point in common.
    bool intersects( const Sphere& _sphere, const Capsule& _capsule );         // TESTED
    inline bool intersects( const Capsule& _capsule, const Sphere& _sphere )          { return intersects(_sphere, _capsule); }
}
