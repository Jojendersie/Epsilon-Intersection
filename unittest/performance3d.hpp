#pragma once

#include <iostream>
#include <vector>


// ************************************************************************* //
// Functions to create reasonable random geometry in the [-1,1] cube
template<class T> void random(T& _out)
{
    throw "Type not implemented!";
}

template<> inline void random<ei::Vec3>(ei::Vec3& _out)
{
    _out = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
}

template<> inline void random<ei::Quaternion>(ei::Quaternion& _out)
{
    float cosTheta = rnd() * 2.0f - 1.0f;
    float phi = rnd() * 2.0f * ei::PI;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    float angle = rnd() * 2.0f * ei::PI;
    _out = ei::Quaternion(ei::Vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi)), angle);
}

template<> inline void random<ei::Sphere>(ei::Sphere& _out)
{
    _out = ei::Sphere(ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f), rnd() * 0.2f + 0.05f);
}

template<> inline void random<ei::Ellipsoid>(ei::Ellipsoid& _out)
{
    _out.center = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.radii = ei::Vec3(rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f);
}

template<> inline void random<ei::OEllipsoid>(ei::OEllipsoid& _out)
{
    _out.center = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.radii = ei::Vec3(rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f);
    random(_out.orientation);
}

template<> inline void random<ei::Ray>(ei::Ray& _out)
{
    _out.origin = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    float cosTheta = rnd() * 2.0f - 1.0f;
    float phi = rnd() * 2.0f * ei::PI;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    _out.direction = ei::Vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
}

template<> inline void random<ei::Box>(ei::Box& _out)
{
    ei::Vec3 w(rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f);
    _out.min = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.max = _out.min + w;
    _out.min = _out.min - w;
}

template<> inline void random<ei::OBox>(ei::OBox& _out)
{
    _out.sides = ei::Vec3(rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f);
    _out.center = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    random(_out.orientation);
}

template<> inline void random<ei::Triangle>(ei::Triangle& _out)
{
    ei::Vec3 pos(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.v0 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v1 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v2 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
}

template<> inline void random<ei::Tetrahedron>(ei::Tetrahedron& _out)
{
    ei::Vec3 pos(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.v0 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v1 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v2 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v3 = pos + ei::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
}


template<> inline void random<ei::Capsule>(ei::Capsule& _out)
{
    _out.seg.a = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.seg.b = ei::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.radius = ei::sq(rnd()) * 0.5f;
}

template<> inline void random<ei::Plane>(ei::Plane& _out)
{
    float cosTheta = rnd() * 2.0f - 1.0f;
    float phi = rnd() * 2.0f * ei::PI;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    _out.n = ei::Vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
    _out.d = rnd() - 0.5f;
}

template<> inline void random<ei::DOP>(ei::DOP& _out)
{
    float cosTheta = rnd() * 2.0f - 1.0f;
    float phi = rnd() * 2.0f * ei::PI;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    _out.n = ei::Vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
    _out.d0 = rnd() - 0.5f;
    _out.d1 = _out.d0 + rnd() - 0.5f;
}

// Functions to assign names to types
template<class T> const char* name() { return typeid(T).name(); }
template<> inline const char* name<ei::Vec3>() { return "Point"; }
template<> inline const char* name<ei::Ray>() { return "Ray"; }
template<> inline const char* name<ei::Sphere>() { return "Sphere"; }
template<> inline const char* name<ei::Ellipsoid>() { return "Ellipsoid"; }
template<> inline const char* name<ei::OEllipsoid>() { return "OEllipsoid"; }
template<> inline const char* name<ei::Box>() { return "Box"; }
template<> inline const char* name<ei::OBox>() { return "Oriented Box"; }
template<> inline const char* name<ei::Triangle>() { return "Triangle"; }
template<> inline const char* name<ei::Tetrahedron>() { return "Tetrahedron"; }
template<> inline const char* name<ei::Capsule>() { return "Capsule"; }
template<> inline const char* name<ei::Plane>() { return "Plane"; }
template<> inline const char* name<ei::DOP>() { return "DOP"; }

// ************************************************************************* //
//                          PERFORMANCE TESTING                              //
// ************************************************************************* //
#ifdef _DEBUG
    const int PERF_ITERATIONS = 100;
    const int TEST_PER_ITERATION = 100;
#else
    const int PERF_ITERATIONS = 5000;
    const int TEST_PER_ITERATION = 1000;
#endif

/// \brief Generic performance testing method for 2 parameters
template<class P0, class P1, class R> void performance(R (*_func)(const P0&, const P1&), const char* _funcName)
{
    float perfIndex = 0.0f;
    uint64 totalTicks = 0;
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        std::vector<P0> geo0(TEST_PER_ITERATION);
        std::vector<P1> geo1(TEST_PER_ITERATION);
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
        {
            random(geo0[i]);
            random(geo1[i]);
        }

        // Start
        uint64 a = ticks();
        volatile R res;
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            res = _func(geo0[i], geo1[i]);
        uint64 b = ticks();
        // Measure dot products to compare
        Vec3* source0 = reinterpret_cast<Vec3*>(geo0.data());
        Vec3* source1 = reinterpret_cast<Vec3*>(geo1.data());
        volatile float xres;
        for(int j = 0; j < 8; ++j)
            for(int i = 0; i < TEST_PER_ITERATION; ++i)
                 xres = sum(cross(source0[i], source1[i]));
        uint64 c = ticks();
        perfIndex += (b-a) * 8.0f / (c-b);
        totalTicks += b-a;
    }
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << "): "
        << perfIndex/PERF_ITERATIONS << " absolute ticks: " << totalTicks / float(PERF_ITERATIONS * TEST_PER_ITERATION) << std::endl;
}

// ************************************************************************* //
/// \brief Generic performance testing method for 3 parameters
template<class P0, class P1, class P2, class R> void performance(R (*_func)(const P0&, const P1&, P2&), const char* _funcName)
{
    float perfIndex = 0.0f;
    uint64 totalTicks = 0;
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        std::vector<P0> geo0(TEST_PER_ITERATION);
        std::vector<P1> geo1(TEST_PER_ITERATION);
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
        {
            random(geo0[i]);
            random(geo1[i]);
        }

        // Start
        uint64 a = ticks();
        volatile R res;
        P2 rt;
        volatile P2 eat;
        for(int i = 0; i < TEST_PER_ITERATION; ++i) {
            res = _func(geo0[i], geo1[i], rt);
            eat = rt;
        }
        uint64 b = ticks();
        // Measure dot products to compare
        Vec3* source0 = reinterpret_cast<Vec3*>(geo0.data());
        Vec3* source1 = reinterpret_cast<Vec3*>(geo1.data());
        volatile float xres;
        for(int j = 0; j < 8; ++j)
            for(int i = 0; i < TEST_PER_ITERATION; ++i)
                 xres = sum(cross(source0[i], source1[i]));
        uint64 c = ticks();
        perfIndex += (b-a) * 8.0f / (c-b);
        totalTicks += b-a;
    }
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << ", " << name<P2>() << "): "
        << perfIndex/PERF_ITERATIONS << " absolute ticks: " << totalTicks / float(PERF_ITERATIONS * TEST_PER_ITERATION) << std::endl;
}

