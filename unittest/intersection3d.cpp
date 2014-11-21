#include "ei/3dintersection.hpp"
#include "unittest.hpp"

#include <iostream>

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

// ************************************************************************* //
// Benchmark for basic intersection types
template<class T0, class T1> void performance()
{
    float perfIndex = 0.0f;
#ifdef _DEBUG
    const int PERF_ITERATIONS = 100;
    const int TEST_PER_ITERATION = 100;
#else
    const int PERF_ITERATIONS = 10000;
    const int TEST_PER_ITERATION = 200;
#endif
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        T0 geo0[TEST_PER_ITERATION];
        T1 geo1[TEST_PER_ITERATION];
        Vec<bool, TEST_PER_ITERATION> res;
        Vec<float, TEST_PER_ITERATION> resDot;
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
        {
            random(geo0[i]);
            random(geo1[i]);
        }

        // Start
        uint64 a = ticks();
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            res[i] = intersects(geo0[i], geo1[i]);
        uint64 b = ticks();
        // Measure dot products to compare
        Vec3* source0 = reinterpret_cast<Vec3*>(geo0);
        Vec3* source1 = reinterpret_cast<Vec3*>(geo1);
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            resDot[i] = dot(source0[i], source1[i]);
        uint64 c = ticks();
        perfIndex += float(b-a)/(c-b);
        eatMyDummy((float)sum(res));
        eatMyDummy(sum(resDot));
    }
    std::cerr << "Performance " << name<T0>() << " <-> " << name<T1>() << ": " << perfIndex/PERF_ITERATIONS << std::endl;
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
        //TEST( intersects( ray1, ell0 ), "ray1 should hit ell0!" );
        TEST( intersects( ray1, ell1 ), "ray1 should hit ell1!" );
        TEST( !intersects( ray2, ell0 ), "ray2 should miss ell0!" );
        TEST( intersects( ray2, ell1 ), "ray2 should hit ell1!" );
        TEST( !intersects( ray3, ell0 ), "ray3 should miss ell0!" );
        TEST( !intersects( ray3, ell1 ), "ray3 should miss ell1!" );

        performance<Ray,Ellipsoid>();
    }

    return result;
}