#pragma once

#include "3dtypes.hpp"

namespace ei {

    /// \brief Descriptor for hit results with boxes.
    /// \details Sides are relative to the local system of the box.
    ///     For OBoxs this is not the absolute coordinate. E.g. if a box is
    ///     rotated by π/2 two coordinates seem to be exchanged.
    enum struct HitSide
    {
        X_NEG = 0x01,         ///< Left.
        X_POS = 0x02,         ///< Right.
        X     = 0x03,         ///< Any side in x-direction
        Y_NEG = 0x04,         ///< Bottom.
        Y_POS = 0x08,         ///< Top.
        Y     = 0x0c,         ///< Any side in y-direction
        Z_NEG = 0x10,         ///< Front.
        Z_POS = 0x20,         ///< Back.
        Z     = 0x30,         ///< Any side in z-direction
    };

    /// \brief Get the euclidean distance between two objects.
    /// \details The distance for point-solid queries can be negative. All
    ///     other geometries return 0 if they intersect.
    float distance(const Vec3& _point0, const Vec3& _point1);                  // TESTED
    float distance(const Vec3& _point, const Segment& _line);                  // TESTED
    float distance(const Vec3& _point, const Triangle& _triangle);             // TESTED
    float distance(const Vec3& _point, const Sphere& _sphere);                 // TESTED
    float distance(const Vec3& _point, const Capsule& _capsule);               // TESTED
    float distance(const Vec3& _point, const Ray& _ray);                       
    float distance(const Vec3& _point, const Box& _box);                       // TESTED
    float distance(const Vec3& _point, const Plane& _plane);                   // TESTED
    /// \returns -x if dot(n,_point) larger then both, x if smaller and 0 if it
    ///     is between both planes
    float distance(const Vec3& _point, const DOP& _dop);                   
    float distance(const Sphere& _sphere, const Segment& _segment);            // TESTED
    float distance(const Sphere& _sphere, const Capsule& _capsule);            // TESTED
    float distance(const Sphere& _sphere, const Box& _box);                    // TESTED
    float distance(const Sphere& _sphere, const Plane& _plane);                // TESTED
    float distance(const Segment& _line0, const Segment& _line1);              // TESTED
    float distance(const Capsule& _capsule0, const Capsule& _capsule1);        // TESTED
    inline float distance(const Segment& _line, const Vec3& _point)            { return distance(_point, _line); }
    inline float distance(const Triangle& _triangle, const Vec3& _point)       { return distance(_point, _triangle); }
    inline float distance(const Sphere& _sphere, const Vec3& _point)           { return distance(_point, _sphere); }
    inline float distance(const Capsule& _capsule, const Vec3& _point)         { return distance(_point, _capsule); }
    inline float distance(const Capsule& _capsule, const Sphere& _sphere)      { return distance(_sphere, _capsule); }
    inline float distance(const Box& _box, const Vec3& _point)                 { return distance(_point, _box); }
    inline float distance(const Box& _box, const Sphere& _sphere)              { return distance(_sphere, _box); }
    inline float distance(const Plane& _plane, const Vec3& _point)             { return distance(_point, _plane); }
    inline float distance(const Plane& _plane, const Sphere& _sphere)          { return distance(_sphere, _plane); }
    inline float distance(const DOP& _dop, const Vec3& _point)                 { return distance(_point, _dop); }
    inline float distance(const Ray& _ray, const Vec3& _point)                 { return distance(_point, _ray); }

    /// \brief Do two spheres intersect, touch or is one inside the other?
    /// \return true if both spheres have at least one point in common.
    bool intersects( const Sphere& _sphere0, const Sphere& _sphere1 );         // TESTED

    /// \brief Does a point lie inside a sphere/on the boundary?
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Sphere& _sphere );              // TESTED
    inline bool intersects( const Sphere& _sphere, const Vec3& _point )  { return intersects( _point, _sphere ); }

    /// \brief Do a sphere and a box intersect?
    /// \return true if the sphere and the box have at least one point in common.
    bool intersects( const Sphere& _sphere, const Box& _box );                 // TESTED
    inline bool intersects( const Box& _box, const Sphere& _sphere )  { return intersects( _sphere, _box ); }

    /// \brief Does a point lie inside a box/on the boundary?
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Box& _box );                    // TESTED
    inline bool intersects( const Box& _box, const Vec3& _point )  { return intersects( _point, _box ); }

    /// \brief Does a point lie inside a box/on the boundary?
    /// \return true if both boxes have at least one point in common.
    bool intersects( const Box& _box0, const Box& _box1 );                     // TESTED

    /// \brief Do two capsules intersect, touch or is one inside the other?
    /// \return true if both spheres have at least one point in common.
    bool intersects( const Capsule& _capsule0, const Capsule& _capsule1 );     // TESTED

    /// \brief Does a point lie inside an ellipsoid/on the boundary?
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const Ellipsoid& _ellipsoid );        // TESTED
    inline bool intersects( const Ellipsoid& _ellipsoid, const Vec3& _point )  { return intersects( _point, _ellipsoid ); }

    /// \brief Does a point lie inside an oriented ellipsoid/on the boundary?
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const OEllipsoid& _oellipsoid );      // TESTED
    inline bool intersects( const OEllipsoid& _oellipsoid, const Vec3& _point )  { return intersects( _point, _oellipsoid ); }

    /// \brief Does a point lie inside a DOP or on the boundary?
    /// \return true if the point is inside or on the boundary.
    bool intersects( const Vec3& _point, const DOP& _dop );                    // TESTED
    inline bool intersects( const DOP& _dop, const Vec3& _point )  { return intersects( _point, _dop ); }

    /// \brief Intersection test between point and oriented box.
    /// \return true if the point and the oriented box have at least on point in common.
    bool intersects( const Vec3& _point, const OBox& _obox );
    inline bool intersects( const OBox& _obox, const Vec3& _point )          { return intersects(_point, _obox); }

    /// \brief Do a ray and a sphere intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \return true if the ray has at least one point in common with thp sphere
    bool intersects( const Ray& _ray, const Sphere& _sphere );
    bool intersects( const Ray& _ray, const Sphere& _sphere, float& _distance );
    inline bool intersects( const Sphere& _sphere, const Ray& _ray ) { return intersects(_ray, _sphere); }
    inline bool intersects( const Sphere& _sphere, const Ray& _ray, float& _distance ) { return intersects(_ray, _sphere, _distance); }

    /// \brief Do a ray and an ellipsoid intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \return true if there is at least one common point between ray and
    ///     ellipsoid.
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid );           // TESTED
    bool intersects( const Ray& _ray, const Ellipsoid& _ellipsoid, float& _distance );
    inline bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray )  { return intersects( _ray, _ellipsoid ); }
    inline bool intersects( const Ellipsoid& _ellipsoid, const Ray& _ray, float& _distance )  { return intersects( _ray, _ellipsoid, _distance ); }

    /// \brief Do a ray and a box intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    ///
    ///     If the ray starts inside 0 is returned.
    /// \param [out,opt] _distanceExit The ray parameter (distance) for the
    ///     second intersection point in positive direction.
    ///
    ///     This is always the exit point. If the ray starts on the boundary
    ///     and shows away _distance and _distanceExit are the same (0).
    /// \param [out,opt] _side Returns which side of the box was hit.
    /// \return true if there is at least one common point between ray and box
    bool intersects( const Ray& _ray, const Box& _box );                       // TESTED
    bool intersects( const Ray& _ray, const Box& _box, float& _distance );
    bool intersects( const Ray& _ray, const Box& _box, float& _distance, HitSide& _side );
    bool intersects( const Ray& _ray, const Box& _box, float& _distance, float& _distanceExit );
    inline bool intersects( const Box& _box, const Ray& _ray )  { return intersects( _ray, _box ); }
    inline bool intersects( const Box& _box, const Ray& _ray, float& _distance )  { return intersects( _ray, _box, _distance ); }
    inline bool intersects( const Box& _box, const Ray& _ray, float& _distance, HitSide& _side )  { return intersects( _ray, _box, _distance, _side ); }
    inline bool intersects( const Box& _box, const Ray& _ray, float& _distance, float& _distanceExit )  { return intersects( _ray, _box, _distance, _distanceExit ); }

    /// \brief Do an oriented box and a ray intersect or touch?
    bool intersects( const Ray& _ray, const OBox& _obox );                     // TESTED
    bool intersects( const Ray& _ray, const OBox& _obox, float& _distance );
    inline bool intersects( const OBox& _obox, const Ray& _ray ) { return intersects( _ray, _obox ); }
    inline bool intersects( const OBox& _obox, const Ray& _ray, float& _distance ) { return intersects( _ray, _obox, _distance ); }

    /// \brief Do a ray and a triangle intersect or touch?
    /// \param [out,opt] _distance The ray parameter (distance) for the first
    ///     intersection point in positive direction.
    /// \param [out,opt] _barycentric The barycentric coordinates of the hit point
    ///     on the triangle.
    /// \return true if there is at least one common point between ray and triangle.
    ///     This point has a ray parameter >= 0 (no negative direction).
    bool intersects( const Ray& _ray, const Triangle& _triangle );                                          // TESTED
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance );                        // TESTED
    bool intersects( const Ray& _ray, const FastTriangle& _triangle, float& _distance );                    // TESTED
    bool intersects( const Ray& _ray, const Triangle& _triangle, float& _distance, Vec3& _barycentric );    // TESTED
    bool intersects( const Ray& _ray, const FastTriangle& _triangle, float& _distance, Vec3& _barycentric );// TESTED
    inline bool intersects( const Triangle& _triangle, const Ray& _ray )  { return intersects( _ray, _triangle ); }
    inline bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance )  { return intersects( _ray, _triangle, _distance ); }
    inline bool intersects( const FastTriangle& _triangle, const Ray& _ray, float& _distance )  { return intersects( _ray, _triangle, _distance ); }
    inline bool intersects( const Triangle& _triangle, const Ray& _ray, float& _distance, Vec3& _barycentric )  { return intersects( _ray, _triangle, _distance, _barycentric ); }
    inline bool intersects( const FastTriangle& _triangle, const Ray& _ray, float& _distance, Vec3& _barycentric )  { return intersects( _ray, _triangle, _distance, _barycentric ); }

    /// \brief Do a sphere and a plane intersect or touch?
    /// \return true if there is at least one point in common.
    bool intersects( const Sphere& _sphere, const Plane& _plane );             // TESTED
    inline bool intersects( const Plane& _plane, const Sphere& _sphere )       { return intersects(_sphere, _plane); }

    /// \brief Does the sphere touches the triangle?
    /// \return true if the sphere and the triangle have at least one point in common.
    bool intersects( const Sphere& _sphere, const Triangle& _triangle );       // TESTED
    inline bool intersects( const Triangle& _triangle, const Sphere& _sphere ) { return intersects(_sphere, _triangle); }

    /// \brief Intersection test between sphere and capsule.
    /// \return true if the sphere and the capsule have at least one point in common.
    bool intersects( const Sphere& _sphere, const Capsule& _capsule );         // TESTED
    inline bool intersects( const Capsule& _capsule, const Sphere& _sphere )   { return intersects(_sphere, _capsule); }

    /// \brief Intersection test between point and capsule.
    /// \return true if the point and the capsule have at least one point in common.
    bool intersects( const Vec3& _point, const Capsule& _capsule );            // TESTED
    inline bool intersects( const Capsule& _capsule, const Vec3& _point )      { return intersects(_point, _capsule); }

    /// \brief Intersection test between point and frustum.
    /// \return true if the point and the frustum have at least one point in common.
    bool intersects( const Vec3& _point, const FastFrustum& _frustum );        // TESTED
    inline bool intersects( const FastFrustum& _frustum, const Vec3& _point )  { return intersects(_point, _frustum); }

    /// \brief Intersection test between sphere and frustum.
    /// \return true if the sphere and the frustum have at least one point in common.
    bool intersects( const Sphere& _sphere, const FastFrustum& _frustum );     // TESTED
    inline bool intersects( const FastFrustum& _frustum, const Sphere& _sphere ) { return intersects(_sphere, _frustum); }

    /// \brief Intersection test between box and frustum.
    /// \return true if the box and the frustum have at least one point in common.
    bool intersects( const Box& _box, const FastFrustum& _frustum );           // TESTED
    inline bool intersects( const FastFrustum& _frustum, const Box& _box )     { return intersects(_box, _frustum); }

    /// \brief Intersection test between point and tetrahedron.
    /// \return true if the point and the tetrahedron have at least one point in common.
    bool intersects( const Vec3& _point, const Tetrahedron& _tetrahedron );    // TESTED
    inline bool intersects( const Tetrahedron& _tetrahedron, const Vec3& _point ) { return intersects(_point, _tetrahedron); }

    /// \brief Intersection test between triangle and box (based on SAT).
    /// \return true if the triangle and the box have at least one point in common.
    bool intersects( const Triangle& _triangle, const Box& _box );             // TESTED
    inline bool intersects( const Box& _box, const Triangle& _triangle ) { return intersects(_triangle, _box); }

    /// \brief Intersection test between triangle and oriented box (based on SAT).
    /// \return true if the triangle and the box have at least one point in common.
    bool intersects( const Triangle& _triangle, const OBox& _obox );           // TESTED
    inline bool intersects( const OBox& _obox, const Triangle& _triangle ) { return intersects(_triangle, _obox); }

    /// \brief Intersection test between plane and box.
    /// \return true if the plane and the box have at least one point in common.
    bool intersects( const Plane& _plane, const Box& _box );
    inline bool intersects( const Box& _box, const Plane& _plane ) { return intersects(_plane, _box); }

    /// \brief Intersection test between plane and oriented box.
    /// \return true if the plane and the oriented box have at least one point in common.
    bool intersects( const Plane& _plane, const OBox& _obox );
    inline bool intersects( const OBox& _obox, const Plane& _plane ) { return intersects(_plane, _obox); }

    /// \brief Intersection test between cone and point.
    bool intersects( const Vec3& _point, const Cone& _cone );
    inline bool intersects( const Cone& _cone, const Vec3& _point ) { return intersects(_point, _cone); }
    bool intersects( const Vec3& _point, const FastCone& _cone );
    inline bool intersects( const FastCone& _cone, const Vec3& _point ) { return intersects(_point, _cone); }

    /// \brief Intersection test between cone and triangle.
    bool intersects( const Triangle& _triangle, const Cone& _cone );
    inline bool intersects( const Cone& _cone, const Triangle& _triangle ) { return intersects(_triangle, _cone); }
}
