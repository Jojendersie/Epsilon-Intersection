#pragma once

#include "3dintersection.hpp"

namespace ei {
    namespace details {
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
        EIAPI void swap(T& _a, T& _b) {
            T t = _a;
            _a = _b;
            _b = t;
        }

        EIAPI uint32 quickHull2D(Vec3* _points, uint32& _numPoints, const Vec3& _normal, Vec3 _a, Vec3 _b)
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
    } // namespace details

    /// \brief Remove all points from the array, which are not part of the
    ///     convex hull.
    /// \details The algorithm moves the points on the convex hull to the front
    ///     of the _points array. The other points are overwritten.
    /// \param [in] _threshold Discard vertices which are closer to the
    ///     previous convex hull than this threshold. This also includes
    ///     duplicates of vertices.
    /// \return Number of points in the convex set (these are the first elements
    ///     in _points after call).
    EIAPI uint32 convexSet(Vec3* _points, uint32 _numPoints, float _threshold = 0.0f)
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
        if(idx[0] == 1) { details::swap(_points[0], _points[1]); details::swap(_points[1], _points[idx[1]]); }
        else { details::swap(_points[1], _points[idx[1]]); details::swap(_points[0], _points[idx[0]]); }
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
        details::swap(_points[2], _points[idx[0]]);
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
            nconvex += details::quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[0], _points[1]);
            nconvex += details::quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[1], _points[2]);
            nconvex += details::quickHull2D(_points + nconvex, _numPoints, groundPlane.n, _points[2], _points[0]);

            delete[] edges;
            return nconvex;
        }

        // Setup memory for a temporary mesh
        details::CHFace* buf = new details::CHFace[_numPoints * 4]; // TODO: sufficient space?
        uint qs = 0, qe = 0, qn = _numPoints * 2;
        uint* queue = new uint[qn];
        uint* visibleStack = new uint[_numPoints * 2 + 2];
        UVec3* edges = new UVec3[_numPoints];

        // Put two opposite faces into our mesh and into a queue of
        // unvisited half spaces.
        buf[0] = details::CHFace(0, 1, 2, UVec3(1), _points);
        buf[1] = details::CHFace(0, 2, 1, UVec3(0), _points);
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
            if(nconvex != idx[0]) details::swap(_points[nconvex], _points[idx[0]]);

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
                        details::swap(_points[i--], _points[--_numPoints]);
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
                        if(i+1 != j) details::swap(edges[i+1], edges[j]);
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

} // namespace ei