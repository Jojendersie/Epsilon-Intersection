#include "ei/3dtypes.hpp"

namespace ei {

    // ********************************************************************* //
    Box::Box( const OBox& _box )
    {
        // Effectively generate all 8 corners and find min/max coordinates.
        // Relative to the center two diagonal opposite corners only differ
        // in the sign (even after rotation).
        Vec3 diag = _box.sides * 0.5f;
        Vec3 trDiag;
        Mat3x3 rot(_box.orientation);
        trDiag = rot * Vec3(diag.x, diag.y, diag.z);
        min = ei::min(trDiag, -trDiag);
        max = ei::max(trDiag, -trDiag);
        trDiag = rot * Vec3(diag.x, diag.y, -diag.z);
        min = ei::min(ei::min(trDiag, -trDiag), min);
        max = ei::max(ei::max(trDiag, -trDiag), max);
        trDiag = rot * Vec3(diag.x, -diag.y, diag.z);
        min = ei::min(ei::min(trDiag, -trDiag), min);
        max = ei::max(ei::max(trDiag, -trDiag), max);
        trDiag = rot * Vec3(diag.x, -diag.y, -diag.z);
        min = ei::min(ei::min(trDiag, -trDiag), min);
        max = ei::max(ei::max(trDiag, -trDiag), max);
        min += _box.center;
        max += _box.center;
    }

    // ********************************************************************* //
    FastFrustum::FastFrustum(const Frustum& _frustum)
    {
        // Initialization of planes is difficult in the list, so the const-cast
        // only defers the initialization a bit.

        // Compute all 8 vertices (first get more help vectors).
        Vec3* v = const_cast<Vec3*>(vertices);
        Vec3 far = _frustum.f * _frustum.direction + _frustum.apex;
        Vec3 near = _frustum.n * _frustum.direction + _frustum.apex;
        // Get third axis and central off point
        Vec3 xAxis = cross(_frustum.up, _frustum.direction);
        Vec3 bottom = _frustum.b*_frustum.up;
        Vec3 top = _frustum.t*_frustum.up;
        Vec3 left =  _frustum.l * xAxis;
        Vec3 right =  _frustum.r * xAxis;
        float fton = _frustum.n / _frustum.f;
        v[0] = near + (left  + bottom) * fton;
        v[1] = near + (left  + top) * fton;
        v[2] = near + (right + bottom) * fton;
        v[3] = near + (right + top) * fton;
        v[4] = far  + left  + bottom;
        v[5] = far  + left  + top;
        v[6] = far  + right + bottom;
        v[7] = far  + right + top;

        // Create planes
        const_cast<DOP&>(nf) = DOP(_frustum.direction, near, far);
        // Use two vectors in the planes to derive the normal.
        const_cast<Plane&>(l) = Plane(normalize(cross(_frustum.up, v[4]-_frustum.apex)), v[4]);
        const_cast<Plane&>(r) = Plane(normalize(cross(v[7]-_frustum.apex, _frustum.up)), v[7]);
        const_cast<Plane&>(b) = Plane(normalize(cross(v[4]-_frustum.apex, xAxis)), v[4]);
        const_cast<Plane&>(t) = Plane(normalize(cross(xAxis, v[7]-_frustum.apex)), v[7]);
    }

}