#include "ei/3dintersection.hpp"
#include "unittest.hpp"

#include <iostream>
#include <vector>

using namespace ei;
using namespace std;

// ************************************************************************* //
// Function to create reasonable random geometry in the [-1,1] cube
template<class T> void random(T& _out)
{
}

template<> void random<Vec3>(Vec3& _out)
{
    _out = Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
}

template<> void random<Ellipsoid>(Ellipsoid& _out)
{
    _out.center = Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.radii = Vec3(rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f, rnd() * 0.2f + 0.05f);
}

template<> void random<Ray>(Ray& _out)
{
    _out.origin = Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    float cosTheta = rnd() * 2.0f - 1.0f;
    float phi = rnd() * 2.0f * PI;
    float sinTheta = sqrt(1 - sq(cosTheta));
    _out.direction = Vec3(sinTheta * cos(phi), cosTheta, sinTheta * sin(phi));
}

template<> void random<Box>(Box& _out)
{
    Vec3 w(rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f);
    _out.min = Vec3(rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f, rnd() * 2.0f - 1.0f);
    _out.max = _out.min + w;
    _out.min = _out.min - w;
}

template<class T> const char* name() { return typeid(T).name(); }
template<> const char* name<Vec3>() { return "Point"; }
template<> const char* name<Ray>() { return "Ray"; }
template<> const char* name<Ellipsoid>() { return "Ellipsoid"; }
template<> const char* name<Box>() { return "Box"; }
template<> const char* name<Triangle>() { return "Triangle"; }

#ifdef _DEBUG
    const int PERF_ITERATIONS = 100;
    const int TEST_PER_ITERATION = 100;
#else
    const int PERF_ITERATIONS = 5000;
    const int TEST_PER_ITERATION = 1000;
#endif

// ************************************************************************* //
// Benchmark for basic intersection types
template<class T0, class T1> void performance()
{
    float perfIndex = 0.0f;
    uint64 totalTicks = 0;
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        std::vector<T0> geo0(TEST_PER_ITERATION);
        std::vector<T1> geo1(TEST_PER_ITERATION);
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
        {
            random(geo0[i]);
            random(geo1[i]);
        }

        // Start
        uint64 a = ticks();
        volatile bool res;
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            res = intersects(geo0[i], geo1[i]);
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
    std::cerr << "Performance " << name<T0>() << " <-> " << name<T1>() << ": "
        << perfIndex/PERF_ITERATIONS << " absolute ticks: " << totalTicks / float(PERF_ITERATIONS * TEST_PER_ITERATION) << std::endl;
}

// ************************************************************************* //
// Benchmark for basic intersection types with single float return value
template<class T0, class T1> void performanceRet1f()
{
    float perfIndex = 0.0f;
    uint64 totalTicks = 0;
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        std::vector<T0> geo0(TEST_PER_ITERATION);
        std::vector<T1> geo1(TEST_PER_ITERATION);
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
        {
            random(geo0[i]);
            random(geo1[i]);
        }

        // Start
        uint64 a = ticks();
        volatile bool res;
        float rt;
        volatile float eat;
        for(int i = 0; i < TEST_PER_ITERATION; ++i) {
            res = intersects(geo0[i], geo1[i], rt);
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
    std::cerr << "Performance " << name<T0>() << " <-> " << name<T1>() << " with ret1f: "
        << perfIndex/PERF_ITERATIONS << " absolute ticks: " << totalTicks / float(PERF_ITERATIONS * TEST_PER_ITERATION) << std::endl;
}


bool test_3dintersections()
{
    bool result = true;

    // Test ellipsoid <-> point intersection
    {
        Vec3 v0( 0.0f, 1.0f, 2.0f );
        Vec3 v1( 1e-8f, 1.0f, 2.0f );
        Ellipsoid ell0( Vec3(0.0f, 1.0f, 0.0f), Vec3(0.1f, 0.0f, 2.0f) );
        TEST( intersects( v0, ell0 ), "Point in degenerated ellipsoid failed!" );
        TEST( !intersects( v1, ell0 ), "Point outside degenerated ellipsoid failed!" );

        performance<Vec3,Ellipsoid>();
    }

    // Test ellipsoid <-> ray intersection
    {
        Ray ray0( Vec3(0.0f, 1.0f, 2.0f), normalize(Vec3(0.2f, -0.9f, 0.0f)) ); /// Starts inside ell1 and misses ell0
        Ray ray1( Vec3(100000.0f, 50.0f, -256.f), normalize(Vec3(-100000.0f, -50.0f, 256.0f)) ); /// Starting far away and intersecting with the origin
        Ray ray2( Vec3(-5.0f, 10.0f, 23.f), normalize(Vec3(5.0f, -9.0f, -22.0f)) ); /// Targets center of ell1
        Ray ray3( Vec3(25.0f, -40.0f, 30.f), normalize(Vec3(5.0f, 1.0f, -2.0f)) );  /// Should not hit
        Ellipsoid ell0( Vec3(0.0f, 0.0f, 0.0f), Vec3(0.5f, 0.0f, 2.0f) );
        Ellipsoid ell1( Vec3(0.0f, 1.0f, 1.0f), Vec3(1.0f, 3.0f, 2.0f) );
        TEST( !intersects( ray0, ell0 ), "ray0 should miss ell0!" );
        TEST( intersects( ray0, ell1 ), "ray0 should hit ell1!" );
        TEST( intersects( ray0.origin, ell1 ), "ray0.origin should be in ell1!" );
        TEST( intersects( Vec3(0.0f, 0.0f, 0.0f), ell1 ), "0 should be in ell1!" );
        // The next one is numerical to unstable
        TEST( intersects( ray1, ell0 ), "ray1 should hit ell0!" );
        TEST( intersects( ray1, ell1 ), "ray1 should hit ell1!" );
        TEST( !intersects( ray2, ell0 ), "ray2 should miss ell0!" );
        TEST( intersects( ray2, ell1 ), "ray2 should hit ell1!" );
        TEST( !intersects( ray3, ell0 ), "ray3 should miss ell0!" );
        TEST( !intersects( ray3, ell1 ), "ray3 should miss ell1!" );

        float d;
        TEST( intersects( ray0, ell1, d ), "ray0 should hit ell1!" );
        //TEST( intersects( ray0.origin + ray0.direction * (d-1e-6f), ell1 ), "Hit point outside ell1!" );
        //TEST( !intersects( ray0.origin + ray0.direction * (d+1e-6f), ell1 ), "Point should be outside ell1!" );
        TEST( intersects( ray2, ell1, d ), "ray2 should hit ell1!" );
        TEST( intersects( ray2.origin + ray2.direction * (d+1e-6f), ell1 ), "Hit point outside ell1!" );
        TEST( !intersects( ray2.origin + ray2.direction * (d-1e-6f), ell1 ), "Point should be outside ell1!" );

        performance<Ray,Ellipsoid>();
        performanceRet1f<Ray,Ellipsoid>();
    }

    // Test box <-> ray intersection
    {
        Ray ray0( Vec3(-1.0f, 0.0f, 0.0f), normalize(Vec3(0.5f, 0.5f, 0.0f)) );
        Ray ray1( Vec3(3.0f, 5.0f, 8.0f), normalize(Vec3(-2.5f, -4.5f, -7.5f)) );
        Ray ray2( Vec3(1.25f, 1.25f, 1.25f), normalize(Vec3(1.5f, 0.5f, 0.5f)) );
        Ray ray3( Vec3(1.25f, 1.25f, 1.25f), normalize(Vec3(0.0f, 1.0f, 0.0f)) );
        Box box0( Vec3(0.0f, 0.0f, 0.0f), Vec3(0.5f, 1.0f, 2.0f) );
        Box box1( Vec3(1.0f, 1.0f, 1.0f), Vec3(1.5f, 1.5f, 1.5f) );
        TEST( intersects( ray0, box0 ), "ray0 should hit box0!" );
        TEST( !intersects( ray0, box1 ), "ray0 should miss box1!" );
        TEST( intersects( ray1, box0 ), "ray1 should hit box0!" );
        TEST( !intersects( ray1, box1 ), "ray1 should miss box1!" );
        TEST( !intersects( ray2, box0 ), "ray2 should miss box0!" );
        TEST( intersects( ray2, box1 ), "ray2 should hit box1!" );
        TEST( intersects( ray3, box1 ), "ray3 should hit box1!" );

        float d;
        //TEST( intersects( ray3, box1, d ) && d == 0.25f, "ray3 should hit box1 with a distance of 0.25!" );
        TEST( intersects( ray0, box0, d ) && d == sqrt(2.0f), "ray0 should hit box0 with a distance of sqrt(2)!" );

        performance<Ray,Box>();
        performanceRet1f<Ray,Box>();
    }

    // Test triangle <-> ray intersection
    {
        Ray ray0( Vec3(-1.0f, 0.0f, 0.0f), normalize(Vec3(0.5f, 0.5f, 0.0f)) );
        Ray ray1( Vec3(90.0f, 100.0f, -110.0f), normalize(Vec3(-88.75f, -99.5f, 111.16666f)) );
        Triangle tri0( Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 2.0f, 1.0f), Vec3(0.0f, 2.0f, -1.0f) );
        Triangle tri1( Vec3(1.0f, 0.0f, 1.0f), Vec3(1.5f, 0.5f, 1.0f), Vec3(1.25f, 1.0f, 1.5f) );
        TEST( intersects( ray0, tri0 ), "ray0 should hit tri0!" );
        TEST( !intersects( ray0, tri1 ), "ray0 should miss tri1!" );
        TEST( !intersects( ray1, tri0 ), "ray1 should miss tri0!" );
        TEST( intersects( ray1, tri1 ), "ray1 should hit tri1!" );

        float d;
        TEST( intersects( ray0, tri0, d ) && d == sqrt(2.0f), "ray0 should hit tri0 with a distance of sqrt(2)!" );
        TEST( intersects( ray1, tri1, d ) && d == 173.593903f, "ray1 should hit tri1 with a distance of 173.593903!" );

        performance<Ray,Triangle>();
        performanceRet1f<Ray,Triangle>();
    }

    return result;
}