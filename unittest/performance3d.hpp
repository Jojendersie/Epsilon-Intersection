#pragma once

#include <iostream>
#include <vector>


// ************************************************************************* //
// Functions to create reasonable random geometry in the [-1,1] cube
template<class T> void random(T& _out)
{
}

template<> inline void random<ε::Vec3>(ε::Vec3& _out)
{
    _out = ε::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
}

template<> inline void random<ε::Sphere>(ε::Sphere& _out)
{
    _out = ε::Sphere(ε::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f), rnd() * 0.2f + 0.05f);
}

template<> inline void random<ε::Ellipsoid>(ε::Ellipsoid& _out)
{
    _out.center = ε::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.radii = ε::Vec3(rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f);
}

template<> inline void random<ε::Ray>(ε::Ray& _out)
{
    _out.origin = ε::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    float cosTheta = rnd() * 2.0f - 1.0f;
    float phi = rnd() * 2.0f * ε::PI;
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    _out.direction = ε::Vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
}

template<> inline void random<ε::Box>(ε::Box& _out)
{
    ε::Vec3 w(rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f);
    _out.min = ε::Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.max = _out.min + w;
    _out.min = _out.min - w;
}

template<> inline void random<ε::Triangle>(ε::Triangle& _out)
{
    ε::Vec3 pos(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.v0 = pos + ε::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v1 = pos + ε::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
    _out.v2 = pos + ε::Vec3(rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f, rnd() * 0.5f - 0.25f);
}

// Functions to assign names to types
template<class T> const char* name() { return typeid(T).name(); }
template<> inline const char* name<ε::Vec3>() { return "Point"; }
template<> inline const char* name<ε::Ray>() { return "Ray"; }
template<> inline const char* name<ε::Sphere>() { return "Sphere"; }
template<> inline const char* name<ε::Ellipsoid>() { return "Ellipsoid"; }
template<> inline const char* name<ε::Box>() { return "Box"; }
template<> inline const char* name<ε::Triangle>() { return "Triangle"; }

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
        for(int j = 0; j < 10; ++j)
            for(int i = 0; i < TEST_PER_ITERATION; ++i)
                 xres = dot(source0[i], source1[i]);
        uint64 c = ticks();
        perfIndex += (b-a) * 10.0f / (c-b);
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
        for(int j = 0; j < 10; ++j)
            for(int i = 0; i < TEST_PER_ITERATION; ++i)
                 xres = dot(source0[i], source1[i]);
        uint64 c = ticks();
        perfIndex += (b-a) * 10.0f / (c-b);
        totalTicks += b-a;
    }
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << ", " << name<P2>() << "): "
        << perfIndex/PERF_ITERATIONS << " absolute ticks: " << totalTicks / float(PERF_ITERATIONS * TEST_PER_ITERATION) << std::endl;
}

