#pragma once

#include "2dtypes.hpp"

namespace ei
{
// TODO: approx, area, center
    float area(const Disc2D& _disc);                                           // TESTED
    float area(const Rect2D& _rect);                                           // TESTED
    float area(const ORect2D& _orect);                                         // TESTED
    float area(const Triangle2D& _triangle);                                   // TESTED
    float area(const Ellipse2D& _ellipse);                                     // TESTED
    float area(const OEllipse2D& _oellipse);                                   // TESTED
    float area(const Segment2D& _segment);                                     // TESTED
    float area(const Ray2D& _ray);                                             // TESTED
    float area(const Capsule2D& _capsule);                                     // TESTED

    Vec2 center(const Disc2D& _disc);                                          // TESTED
    Vec2 center(const Rect2D& _rect);                                          // TESTED
    Vec2 center(const ORect2D& _orect);                                        // TESTED
    Vec2 center(const Triangle2D& _triangle);                                  // TESTED
    Vec2 center(const Ellipse2D& _ellipse);                                    // TESTED
    Vec2 center(const OEllipse2D& _oellipse);                                  // TESTED
    Vec2 center(const Segment2D& _segment);                                    // TESTED
    Vec2 center(const Capsule2D& _capsule);                                    // TESTED

    // Include inline implementations
#   include "details/2dfunctions.inl"
}