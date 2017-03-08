#include "ei/3dtypes.hpp"
#include "ei/3dintersection.hpp"

namespace ei {

    // ********************************************************************* //
    Sphere::Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2 )
    {
        // The center of the circumscribed circle is at (barycentric coordinates)
        // v0*sin(2 alpha) + v1*sin(2 beta) + v2*sin(2 gamma) and has the radius
        // abc/4A.
        Vec3 c = _p0 - _p1;    float csq = lensq(c);
        Vec3 a = _p1 - _p2;    float asq = lensq(a);
        Vec3 b = _p2 - _p0;    float bsq = lensq(b);

        // One of the sides could be the longest side - the minimum sphere is
        // defined through only two points.
        // This can also handle the coplanar case.
        if( csq + bsq <= asq ) *this = Sphere(_p1, _p2);
        else if( asq + bsq <= csq ) *this = Sphere(_p1, _p0);
        else if( asq + csq <= bsq ) *this = Sphere(_p2, _p0);
        else {
            float area2Sq = 2*lensq(cross(a, c));
            center =
                     _p0 * (-dot(c,b)*asq/area2Sq)
                   + _p1 * (-dot(c,a)*bsq/area2Sq)
                   + _p2 * (-dot(b,a)*csq/area2Sq);
            radius = sqrt(asq*bsq*csq/(2*area2Sq));
        }
    }

    // ********************************************************************* //
    Sphere::Sphere( const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const Vec3& _p3 )
    {
        // It is possible that not all 4 points lie on the surface of the sphere.
        // Just two of them could already define a sphere enclosing all other.
        // So we need to compute any combination of possible spheres (14), but
        // luckily we know a direct solution for any combination of 3 points.
        // This reduces the work to 4 cases: build a bounding sphere for 3 points
        // and have a look if the fourth point is inside.
        *this = Sphere( _p1, _p2, _p3 );
        float rsq = radius*radius;
        if( lensq(_p3 - center) > rsq ) {
            *this = Sphere( _p0, _p1, _p3 );
            if( lensq(_p2 - center) > rsq ) {
                *this = Sphere( _p0, _p2, _p3 );
                if( lensq(_p1 - center) > rsq ) {
                    *this = Sphere( _p1, _p2, _p3 );
                    if( lensq(_p0 - center) > rsq ) {
                        // All 4 points are on the boundary -> construct sphere
                        // from 4 points.
                        Vec3 a = _p1 - _p0;
                        Vec3 b = _p2 - _p0;
                        Vec3 c = _p3 - _p0;

                        Mat3x3 m = Mat3x3( _p1[0], _p1[1], _p1[2],
                                           _p2[0], _p2[1], _p2[2],
                                           _p3[0], _p3[1], _p3[2] );

                        float denominator = 0.5f / determinant( m );

                        Vec3 o = (lensq(c) * cross(a,b) +
                                  lensq(b) * cross(c,a) +
                                  lensq(a) * cross(b,c)) * denominator;

                        center = _p0 + o;
                        radius = len(o);
                    }
                }
            }
        }
    }

    // ********************************************************************* //
    struct SingleLinkedPointList
    {
        Vec3 p;
        int next;   ///< -1 for the last element
    };

    static Sphere minimalBoundingSphere( SingleLinkedPointList* _points, uint32 _first, uint32 _n, uint32 _boundarySet )
    {
        Sphere mbs;
        eiAssertWeak(_boundarySet > 0, "Expected at least one point.");

        // If the boundary list is full or all points where added stop
        switch(_boundarySet)
        {
        case 1: mbs = Sphere(_points[_first].p, 0.0f);
            break;
        case 2: {
            Vec3 v0 = _points[_first].p; uint32 next = _points[_first].next;
            mbs = Sphere(v0, _points[next].p);
            break;
        }
        case 3: {
            Vec3 v0 = _points[_first].p; uint32 next = _points[_first].next;
            Vec3 v1 = _points[next].p; next = _points[next].next;
            mbs = Sphere(v0, v1, _points[next].p);
            break;
        }
        case 4: {
            Vec3 v0 = _points[_first].p; uint32 next = _points[_first].next;
            Vec3 v1 = _points[next].p; next = _points[next].next;
            Vec3 v2 = _points[next].p; next = _points[next].next;
            return Sphere(v0, v1, v2, _points[next].p);
        }
        }

        uint32 it = _first;
        uint32 last = _first; // At the beginning this is wrong but does not cause damage, afterwards it will be ovewritten by real numbers
        for(uint32 i = 0; i < _boundarySet; ++i) {last = it; it = _points[it].next;}
        for(uint32 i = _boundarySet; i < _n; ++i)
        {
            // Save next pointer to advance from this point even if the list is changed
            uint32 next = _points[it].next;
            eiAssert(it != 0xffffffff, "Iteration should not have stopped.");
            if(lensq(mbs.center - _points[it].p) > mbs.radius * mbs.radius)
            {
                // Move to front
                _points[last].next = _points[it].next;
                _points[it].next = _first;
                _first = it;
                // Rebuild the first i elements
                mbs = minimalBoundingSphere(_points, it, i, _boundarySet+1);
            }
            last = it;
            it = next;
        }

        return mbs;
    }

    Sphere::Sphere( const Vec3* _points, uint32 _numPoints )
    {
        eiAssert( _numPoints > 0, "The point list must have at least one point." );
        // Create a single linked list for the move to front heuristic
        SingleLinkedPointList* list = new SingleLinkedPointList[_numPoints];
        for(uint32 i = 0; i < _numPoints; ++i) {
            list[i].p = _points[i];
            list[i].next = i+1;
        }
        list[_numPoints-1].next = -1;
        *this = minimalBoundingSphere(list, 0, _numPoints, 1);
        delete[] list;
    }

    // ********************************************************************* //
    Box::Box( const OBox& _box )
    {
        // Effectively generate all 8 corners and find min/max coordinates.
        // Relative to the center two diagonal opposite corners only differ
        // in the sign (even after rotation).

        Vec3 diag = (_box.sides * 0.5f);
        Mat3x3 rot(_box.orientation);
        // Rows of the matrix are the aabox face directions in obox-space.
        // Choose diagonal entries via sign to get the largest coordinate into
        // face direction. Then project this maxCoord onto the direction.
        //Vec3 maxCoord = sgn(invRot(0)) * diag;
        //max.x = dot(maxCoord, invRot(0));
        //max.x = dot(diag, abs(rot(0))); // Equivalent to the sign stuff
        max = abs(rot) * diag;
        min = -max;

        min += _box.center;
        max += _box.center;
    }

    // ********************************************************************* //
    Box::Box( const Vec3* _points, uint32 _numPoints )
    {
        eiAssert( _numPoints > 0, "The point list must have at least one point." );
        min = max = *_points++;
        for( uint32 i = 1; i < _numPoints; ++i, ++_points )
        {
            min = ei::min(min, *_points);
            max = ei::max(max, *_points);
        }
    }

    // ********************************************************************* //
    OBox::OBox( const Quaternion& _orientation, const Box& _box ) :
        center((_box.min + _box.max) * 0.5f),
        orientation(_orientation)
    {
        // Use the same farthest plane search like in Box(OBox), but using
        // columns instead rows.
        Mat3x3 rotation(~_orientation);
        Vec3 bmax = _box.max - center;
        sides = abs(rotation) * bmax;

        sides *= 2.0f;
    }

    // ********************************************************************* //
    OBox::OBox( const Mat3x3& _orientation, const Box& _box ) :
        center((_box.min + _box.max) * 0.5f),
        orientation(_orientation)
    {
        // Use the same farthest plane search like in Box(OBox), but using
        // columns instead rows.
        Vec3 bmax = _box.max - center;
        sides = abs(transpose(_orientation)) * bmax;

        sides *= 2.0f;
    }

    // ********************************************************************* //
    OBox::OBox( const Quaternion& _orientation, const Vec3* _points, uint32 _numPoints ) :
        orientation(_orientation)
    {
        eiAssert( _numPoints > 0, "The point list must have at least one point." );

        // Project all points to the cube sides by transforming them into local
        // space, such that the box is axis aligned again.
        Mat3x3 rotation(conjugate(_orientation));
        Vec3 min, max;
        min = max = rotation * _points[0];
        for(uint32 i = 1; i < _numPoints; ++i)
        {
            Vec3 p = rotation * _points[i];
            min = ei::min(min, p);
            max = ei::max(max, p);
        }

        // Center known with respect to local rotation, go back to world space.
        center = transform((min + max) * 0.5f, conjugate(_orientation));
        sides = max - min;
    }

    // ********************************************************************* //
    OBox::OBox( const Vec3* _points, uint32 _numPoints )
    {
        if(_numPoints == 1)
        {
            sides = Vec3(0.0f);
            center = *_points;
            orientation = qidentity();
        } else if(_numPoints == 2) {
            Vec3 connection = _points[1] - _points[0];
            sides = Vec3(len(connection), 0.0f, 0.0f);
            orientation = Quaternion(connection/sides.x, Vec3(1.0f, 0.0f, 0.0f));
            center = _points[0] + 0.5f * connection;
        } else {
            sides = Vec3(1e12f);
            float volume = 1e36f;
            Mat3x3 finalRotation;
            // Try each combination of three vertices to setup an orientation
            for(uint32 i = 0; i < _numPoints-2; ++i)
            {
                for(uint32 j = i+1; j < _numPoints-1; ++j)
                {
                    Vec3 xAxis = normalize(_points[i] - _points[j]);
                    for(uint32 k = j+1; k < _numPoints; ++k)
                    {
                        Mat3x3 rotation;
                        Vec3 yAxis = cross(xAxis, _points[i] - _points[k]);
                        float l = len(yAxis);
                        if( l < 1e-6f ) // Colinear points
                            rotation = Mat3x3(Quaternion(xAxis, Vec3(1.0f, 0.0f, 0.0f)));
                        else {
                            yAxis /= l;
                            rotation = Mat3x3(xAxis, yAxis, cross(xAxis, yAxis));
                        }
                        // Refit a box with the current rotation. Since the rotation of
                        // all points is the most expensive part try to early out.
                        Vec3 min, max;
                        min = max = rotation * _points[0];
                        for(uint32 i = 1; i < _numPoints; ++i)
                        {
                            Vec3 p = rotation * _points[i];
                            if(p.x < min.x) min.x = p.x;
                            if(p.y < min.y) min.y = p.y;
                            if(p.z < min.z) min.z = p.z;
                            if(p.x > max.x) max.x = p.x;
                            if(p.y > max.y) max.y = p.y;
                            if(p.z > max.z) max.z = p.z;
                            if(prod(max-min) >= volume) goto NextRotation;
                        }
                        // The new box is better than the current, otherwise
                        // the goto would have skipped this section before.
                        //center = transform((min + max) * 0.5f, conjugate(testOrientation));
                        center = transpose(rotation) * ((min + max) * 0.5f);
                        sides = max - min;
                        finalRotation = rotation;
                        //orientation = conjugate(Quaternion(rotation));
                        volume = prod(sides);
                        NextRotation:;
                    }
                }
            }
            orientation = conjugate(Quaternion(finalRotation));
        }
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


    // ********************************************************************* //
    // Test if a point _x is inside a triangle (_v0, _v1, _v2) slab with a
    // thickness of _threshold.
    // Since the normal (and therefore the plane) must be computed anyway it
    // is returned as additional value.
    /*static bool isConvexComb(const Vec3& _x, const Vec3& _v0, const Vec3& _v1, const Vec3& _v2, float _threshold, Plane& _outPlane, float& _d)
    {
        Vec3 a = _v1 - _v0;
        Vec3 b = _v2 - _v0;
        Vec3 scaledNormal = cross(b, a);
        float area2 = len(scaledNormal);
        _outPlane.n = scaledNormal / area2;
        _outPlane.d = -dot(_outPlane.n, _v0);

        // If point is far away from the plane we can early out
        _d = dot(_outPlane.n, _x) + _outPlane.d;
        if(abs(_d) > _threshold)
            return false;

        Vec3 o = _v0 - _x;
        Vec3 d = cross( _outPlane.n, o );

        float barycentricCoord1 = dot( d, b ) / area2;
        if(barycentricCoord1 < 0.0f) return false;
        float barycentricCoord2 = dot( d, a ) / area2;
        if(barycentricCoord2 < 0.0f) return false;
        return barycentricCoord1 + barycentricCoord2 <= 1.0f;
    }*/


    // Helper struct for faces on the convex hull.
    struct CHFace
    {
        UVec3 indices;
        UVec3 neighbors; // triangle neighbor indices for the edges (0,1), (1,2), (2,0).
        Plane p;
        bool isDeleted;
        CHFace() = default;
        CHFace(uint32 _i0, uint32 _i1, uint32 _i2, const UVec3& _neighbors, const Vec3* _points/*, uint32 _insidePoint*/) :
            indices(_i0, _i1, _i2),
            neighbors(_neighbors),
            p(_points[_i0], _points[_i1], _points[_i2]),
            isDeleted(false)
        {
        }
    };

    template<typename T>
    static void swap(T& _a, T& _b)
    {
        T tmp = _a;
        _a = _b;
        _b = tmp;
    }

    uint32 quickHull2D(Vec3* _points, uint32& _numPoints, const Vec3& _normal, Vec3 _a, Vec3 _b)
    {
        if(_numPoints == 0) return 0;
        // Find the furthest point on the right side relative to line a-b
        uint32 idx = _numPoints;
        Vec3 e0 = _b - _a;
        Vec3 edgeNormal = cross(e0, _normal);
        float dmax = 0.0f;
        for(uint32 i = 0; i < _numPoints; ++i)
        {
            float d = dot(_points[i] - _a, edgeNormal);
            if(d > dmax) { dmax = d; idx = i; }
            // TODO: list could be sorted after sign. Only positive side must be investigated by recursive calls.
        }

        if(idx != _numPoints)
        {
            if(idx != 0) swap(_points[0], _points[idx]);
            // Discard all points within the triangle a,b,p[i].
            Vec3 e1 = _points[0] - _a;
            float d00 = dot(e0, e0);
            float d01 = dot(e0, e1);
            float d11 = dot(e1, e1);
            float denom = d00 * d11 - d01 * d01;
            for(uint32 i = 1; i < _numPoints; ++i)
            {
                float d20 = dot(_points[i] - _a, e0);
                float d21 = dot(_points[i] - _a, e1);
                float bary0 = (d11 * d20 - d01 * d21) / denom;
                float bary1 = (d00 * d21 - d01 * d20) / denom;
                if((bary0 >= 0.0f) && (bary1 >= 0.0f) && (bary0 + bary1 <= 1.0f))
                    swap(_points[i--], _points[--_numPoints]);
            }

            // Call recursive for the two new edges
            uint32 nconvex = 1;
            eiAssert(_numPoints > 0, "At least the new extreme point must exist!");
            _numPoints--; // Decrement because offsetting requires smaller bounds
            nconvex += quickHull2D(_points + nconvex, _numPoints, _normal, _a, _points[idx]);
            nconvex += quickHull2D(_points + nconvex, _numPoints, _normal, _points[idx], _b);
            return nconvex;
        }
        return 0;
    }

    uint32 convexSet(Vec3* _points, uint32 _numPoints, float _threshold)
    {
        float tSq = _threshold * _threshold; // Compare squared numbers

        // Find duplicate points (brute force)
        for(uint32 i = 0; i < _numPoints-1; ++i)
        {
            for(uint32 j = i+1; j < _numPoints; ++j)
                if(lensq(_points[j] - _points[i]) <= tSq)
                {
                    // Whoops, j is the same as i -> remove j
                    _points[j--] = _points[--_numPoints];
                }
        }

        if(_numPoints <= 2) return _numPoints;

        // Find some extremal points for initial triangle.
        // First get two out of the six which define the AABox.
        uint32 extrema[6] = {0};
        for(uint32 i = 1; i < _numPoints; ++i)
        {
            if(_points[i].x < _points[extrema[0]].x) extrema[0] = i;
            if(_points[i].x > _points[extrema[1]].x) extrema[1] = i;
            if(_points[i].y < _points[extrema[2]].y) extrema[2] = i;
            if(_points[i].y > _points[extrema[3]].y) extrema[3] = i;
            if(_points[i].z < _points[extrema[4]].z) extrema[4] = i;
            if(_points[i].z > _points[extrema[5]].z) extrema[5] = i;
        }
        uint32 idx[2];
        float dmax = 0.0f;
        for(uint i = 0; i < 5; ++i)
        {
            for(uint j = i+1; j < 6; ++j)
            {
                float d = lensq(_points[extrema[i]] - _points[extrema[j]]);
                if(d > dmax)
                {
                    dmax = d;
                    idx[0] = extrema[i];
                    idx[1] = extrema[j];
                }
            }
        }
        eiAssert(dmax > tSq, "Points should have been deleted in duplicate search.");
        // Move the two points to the front. Order could be an issue if idx[1] == 0
        // or idx[0] == 1. If there is a cross reference swap that element first.
        if(idx[0] == 1) { swap(_points[0], _points[1]); swap(_points[1], _points[idx[1]]); }
        else { swap(_points[1], _points[idx[1]]); swap(_points[0], _points[idx[0]]); }
        // Get third point as most distant one to the line of the first two.
        dmax = 0.0f;
        Segment tsegment(_points[0], _points[1]);
        for(uint32 i = 2; i < _numPoints; ++i)
        {
            float d = distance(_points[i], tsegment);
            if(d > dmax)
            {
                dmax = d;
                idx[0] = i;
            }
        }
        if(dmax <= _threshold)
        {
            // All other points are convex combinations of the two extrema
            // (colinear point set).
            return 2;
        }
        swap(_points[2], _points[idx[0]]);
        if(_numPoints == 3) return 3;

        // Precondition: we have 3 non-colinear points in the front of the list.
        uint32 nconvex = 3;

        // All coplanar?
        Plane groundPlane(_points[0], _points[1], _points[2]);
        bool coplanar = true;
        for(uint32 i = nconvex; i < _numPoints && coplanar; ++i)
        {
            float d = dot(groundPlane.n, _points[i]) + groundPlane.d;
            if( abs(d) > _threshold ) coplanar = false;
        }
        if(coplanar)
        {
            // Perform an 2D variant of the quickhull algorithm.
            UVec2* edges = new UVec2[_numPoints];
            // Start with a loop of three edges
            _numPoints -= nconvex; // Decrement because offsetting requires smaller bounds
            nconvex += quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[0], _points[1]);
            nconvex += quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[1], _points[2]);
            nconvex += quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[2], _points[0]);

            delete[] edges;
            return nconvex;
        }

        // Setup memory for a temporary mesh
        CHFace* buf = new CHFace[_numPoints * 4]; // TODO: sufficient space?
        uint qs = 0, qe = 0, qn = _numPoints * 2;
        uint* queue = new uint[qn];
        uint* visibleStack = new uint[_numPoints * 2 + 2];
        UVec3* edges = new UVec3[_numPoints];

        // Put two opposite faces into our mesh and into a queue of
        // unvisited half spaces.
        buf[0] = CHFace(0, 1, 2, UVec3(1), _points);
        buf[1] = CHFace(0, 2, 1, UVec3(0), _points);
        uint nf = 2; // number of faces
        queue[qe] = 0; qe = (qe + 1) % qn;
        queue[qe] = 1; qe = (qe + 1) % qn;

        while(qs != qe && nconvex < _numPoints)
        {
            // Take a face from the queue. It can happen that deleted faces are listed.
            uint fidx = queue[qs]; qs = (qs + 1) % qn;
            while(buf[fidx].isDeleted && (qs != qe)) { fidx = queue[qs]; qs = (qs + 1) % qn; }
            if(buf[fidx].isDeleted && (qs == qe)) break;

            // Find the maximum distance point on the positive side.
            dmax = 0.0f;
            for(uint32 i = nconvex; i < _numPoints; ++i)
            {
                float d = dot(buf[fidx].p.n, _points[i]) + buf[fidx].p.d;
                if( d > dmax ) { dmax = d; idx[0] = i; }
            }
            if(dmax <= _threshold) continue;
            if(nconvex != idx[0]) swap(_points[nconvex], _points[idx[0]]);

            // For the new point find the set of visible faces.
            visibleStack[0] = fidx;
            int sp = 1;
            int ep = 0;
            while(sp > 0)
            {
                uint vidx = visibleStack[--sp];
                while(buf[vidx].isDeleted && sp > 0) vidx = visibleStack[--sp];
                if(buf[vidx].isDeleted && (sp <= 0)) break;

                // The point and the current visible triangle span a tetrahedron.
                // Remove all points inside the tetrahedron.
                Plane sides[4];
                sides[0].n = -buf[vidx].p.n; sides[0].d = -buf[vidx].p.d;
                sides[1] = Plane(_points[buf[vidx].indices.x], _points[buf[vidx].indices.y], _points[nconvex]);
                sides[2] = Plane(_points[buf[vidx].indices.y], _points[buf[vidx].indices.z], _points[nconvex]);
                sides[3] = Plane(_points[buf[vidx].indices.z], _points[buf[vidx].indices.x], _points[nconvex]);
                for(uint32 i = nconvex+1; i < _numPoints; ++i)
                {
                    bool inside = true;
                    for(int j = 0; j < 4 && inside; ++j)
                        inside &= (dot(sides[j].n, _points[i]) + sides[j].d) <= 0.0f;
                    if(inside)
                        swap(_points[i--], _points[--_numPoints]);
                }

                // Remove the visible face. It is now interior.
                buf[vidx].isDeleted = true;

                // Check the (up to) three neighbors.
                for(int i = 0; i < 3; ++i)
                {
                    uint32 nidx = buf[vidx].neighbors[i];
                    bool seesNeighbor = (dot(buf[nidx].p.n, _points[nconvex]) + buf[nidx].p.d) > 0.0f;

                    if(seesNeighbor)
                    {
                        visibleStack[sp++] = nidx;
                    } else {
                        // If there is an edge between visible and invisible face it will
                        // span a new face with the current point. Add to horizon list.
                        edges[ep++] = UVec3(buf[vidx].indices[i], buf[vidx].indices[(i+1)%3], nidx);
                    }
                }
            }

            // Sort all edges in the horizon list such that they give a loop.
            for(int i = 0; i < ep-1; ++i)
            {
                for(int j = i+1; j < ep; ++j)
                {
                    if(edges[i].y == edges[j].x)
                    {
                        if(i+1 != j) swap(edges[i+1], edges[j]);
                        break;
                    }
                }
            }

            // Create new faces to the horizon loop. To define the adjacencies we
            // need the sorted order.
            for(int e = 0; e < ep; ++e)
            {
                buf[nf + e].indices = UVec3(edges[e].x, edges[e].y, nconvex);
                buf[nf + e].isDeleted = false;
                buf[nf + e].neighbors = UVec3(edges[e].z, nf + (e + 1) % ep, nf + (e + ep - 1) % ep);
                buf[nf + e].p = Plane(_points[edges[e].x], _points[edges[e].y], _points[nconvex]);
                queue[qe] = nf + e; qe = (qe + 1) % qn;
            }
            nf += ep;

            nconvex++;
        }

        delete[] buf;
        delete[] queue;
        delete[] visibleStack;
        delete[] edges;
        return nconvex;
    }
}
