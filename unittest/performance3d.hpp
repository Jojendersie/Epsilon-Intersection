#pragma once

#include <iostream>
#include <vector>
#include <typeinfo>
#include <algorithm>
#include <numeric>


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
    _out.halfSides = ei::Vec3(rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f, rnd() * 0.3f + 0.05f) * 0.5f;
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

template<> inline void random<ei::FastTriangle>(ei::FastTriangle& _out)
{
    ei::Triangle rndTriangle;
    random<ei::Triangle>(rndTriangle);
    _out = ei::FastTriangle(rndTriangle);
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

template<> inline void random<ei::Cone>(ei::Cone& _out)
{
    random<ei::Ray>(_out.centralRay);
    _out.tanTheta = rnd();
    _out.height = ei::sq(rnd()) * 3.0f;
}

template<> inline void random<ei::FastCone>(ei::FastCone& _out)
{
    ei::Cone rndCone;
    random<ei::Cone>(rndCone);
    _out = ei::FastCone(rndCone);
}

// Functions to assign names to types
template<class T> const char* name() { return typeid(T).name(); }
template<> inline const char* name<ei::Vec3>()          { return "       Point"; }
template<> inline const char* name<ei::Ray>()           { return "         Ray"; }
template<> inline const char* name<ei::Cone>()          { return "        Cone"; }
template<> inline const char* name<ei::FastCone>()      { return "   Fast Cone"; }
template<> inline const char* name<ei::Sphere>()        { return "      Sphere"; }
template<> inline const char* name<ei::Ellipsoid>()     { return "   Ellipsoid"; }
template<> inline const char* name<ei::OEllipsoid>()    { return "  OEllipsoid"; }
template<> inline const char* name<ei::Box>()           { return "         Box"; }
template<> inline const char* name<ei::OBox>()          { return "Oriented Box"; }
template<> inline const char* name<ei::Triangle>()      { return "    Triangle"; }
template<> inline const char* name<ei::FastTriangle>()  { return "FastTriangle"; }
template<> inline const char* name<ei::Tetrahedron>()   { return " Tetrahedron"; }
template<> inline const char* name<ei::Capsule>()       { return "     Capsule"; }
template<> inline const char* name<ei::Plane>()         { return "       Plane"; }
template<> inline const char* name<ei::DOP>()           { return "         DOP"; }

// ************************************************************************* //
//                          PERFORMANCE TESTING                              //
// ************************************************************************* //
#ifdef _DEBUG
    const int PERF_ITERATIONS = 100;
    const int TEST_PER_ITERATION = 10;
#else
    const int PERF_ITERATIONS = 5000;
    const int TEST_PER_ITERATION = 50;
#endif

/// \brief Preparation of performance measurements to get stable timings later.
/// The idea is to determine relative performance (much more stable) and multiply
/// that with the time for the reference test. This reference time is determined
/// by the following method.
inline float getTimePer1MElements()
{
    static float s_timePer1MElements;
    if(s_timePer1MElements == 0.0f) // Not determined yet?
    {
        uint64 totalTicks = 0;
        uint64 totalTicksSq = 0;
        const int iterations = PERF_ITERATIONS * 20;
        const int ntest = TEST_PER_ITERATION;
        std::vector<uint64> iterationTime(iterations);
        for(int t = 0; t < iterations; ++t)
        {
            std::vector<ei::Vec3> v0(ntest);
            std::vector<ei::Vec3> v1(ntest);
            for(int i = 0; i < ntest; ++i)
            {
                random(v0[i]);
                random(v1[i]);
            }
            volatile float xres;
            uint64 a = ticks();
            for(int i = 0; i < ntest; ++i)
                for(int j = 0; j < ntest; ++j)
                    xres = lensq(cross(v0[i], v1[j]));
            uint64 b = ticks();
            totalTicks += b-a;
            totalTicksSq += (b-a) * (b-a);
            iterationTime[t] = b-a;
        }
        std::sort(iterationTime.begin(), iterationTime.end());
        // Tons of tests to get the most stable solution
        // PROBLEM: the current measurement has a stddev between 25% and 50%.
        std::cerr << "Reference timing median: " << iterationTime[iterations/2] / float(ntest * ntest) << std::endl;
        /*float ticksPerTest = totalTicks / float(iterations);
        float tickVariance = totalTicksSq / float(iterations) - ticksPerTest * ticksPerTest;
        ticksPerTest /= ntest * ntest;
        float tickStddev = sqrt(tickVariance) / (ntest * ntest);
        std::cerr << "Timing for reference test: " << ticksPerTest << " +- " << tickStddev << std::endl;
        std::cerr << "                           " << s_timePer1MElements << " ms per 1M elements" << std::endl;
        double mean = std::accumulate(iterationTime.begin() + iterations/4, iterationTime.begin() + iterations*3/4, 0.0);
        mean /= iterations/2;
        double s = 0.0;
        for(int i = iterations/4; i < iterations*3/4; ++i)
            s += ei::sq(iterationTime[i] - mean);
        s /= iterations/2 - 1;
        std::cerr << "Timing for reference test 2: " << mean / (ntest * ntest) << " +- " << sqrt(s) / (ntest * ntest) << std::endl;*/

        // The winner is: the median
        s_timePer1MElements = deltaTicksToMilliSeconds(ei::uint64(iterationTime[iterations/2] * 1000000.0 / (ntest * ntest)));
    }
    return s_timePer1MElements;
}

/// \brief Generic performance testing method for 2 parameters
template<class P0, class P1, class R> void performance(R (*_func)(const P0&, const P1&), const char* _funcName)
{
    /*double perfIndex = 0.0f;
    double perfIndexSq = 0.0f;
    uint64 totalTicks = 0;
    uint64 totalTicksSq = 0;*/
    std::vector<uint64> iterationTimes(PERF_ITERATIONS);
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        // Since the Fast... types have no default ctor I use the dummyMem+cast here.
        char dummyMem[512];
        std::vector<P0> geo0(TEST_PER_ITERATION, *reinterpret_cast<P0*>(dummyMem));
        std::vector<P1> geo1(TEST_PER_ITERATION, *reinterpret_cast<P1*>(dummyMem));
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
        {
            random(geo0[i]);
            random(geo1[i]);
        }

        // Start
        uint64 a = ticks();
        volatile R res;
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            for(int j = 0; j < TEST_PER_ITERATION; ++j)
                res = _func(geo0[i], geo1[j]);
        uint64 b = ticks();
        // Measure cross products to compare
        /*ei::Vec3* source0 = reinterpret_cast<ei::Vec3*>(geo0.data());
        ei::Vec3* source1 = reinterpret_cast<ei::Vec3*>(geo1.data());
        volatile float xres;
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            for(int j = 0; j < TEST_PER_ITERATION; ++j)
                xres = lensq(cross(source0[i], source1[j]));
        uint64 c = ticks();
        perfIndex += (b-a) / double(c-b);
        perfIndexSq += ei::sq((b-a) / double(c-b));
        totalTicks += b-a;
        totalTicksSq += (b-a) * (b-a);*/
        iterationTimes[t] = b-a;
        //std::cerr << (b-a) << std::endl;
    }
    /*float ticksPerTest = totalTicks / float(PERF_ITERATIONS);
    float tickVariance = totalTicksSq / float(PERF_ITERATIONS) - ticksPerTest * ticksPerTest;
    ticksPerTest /= TEST_PER_ITERATION * TEST_PER_ITERATION;
    float tickStddev = sqrt(tickVariance) / (TEST_PER_ITERATION * TEST_PER_ITERATION);
    float perfPerTest = float(perfIndex / PERF_ITERATIONS);
    float perfVariance = float(perfIndexSq / PERF_ITERATIONS) - perfPerTest * perfPerTest;
    float perfStddev = sqrt(perfVariance);// / (TEST_PER_ITERATION * TEST_PER_ITERATION);
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << "):\t"
        //<< getTimePer1MElements() * perfPerTest << " ms/M\t"
        << deltaTicksToMilliSeconds(ei::uint64(ticksPerTest * 1000000.0))
        << perfPerTest << " +- " << perfStddev << " absolute ticks: " << ticksPerTest
        << " +- " << tickStddev << std::endl;*/

    std::sort(iterationTimes.begin(), iterationTimes.end());
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << "):\t"
        << deltaTicksToMilliSeconds(ei::uint64(iterationTimes[PERF_ITERATIONS/2] * 1000000.0 / (TEST_PER_ITERATION * TEST_PER_ITERATION)))
        << " ms/M\n";
}

// ************************************************************************* //
/// \brief Generic performance testing method for 3 parameters
template<class P0, class P1, class P2, class R> void performance(R (*_func)(const P0&, const P1&, P2&), const char* _funcName)
{
    /*double perfIndex = 0.0f;
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
        // Measure cross products to compare
        ei::Vec3* source0 = reinterpret_cast<ei::Vec3*>(geo0.data());
        ei::Vec3* source1 = reinterpret_cast<ei::Vec3*>(geo1.data());
        volatile float xres;
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
             xres = lensq(cross(source0[i], source1[i]));
        uint64 c = ticks();
        perfIndex += (b-a) / double(c-b);
        totalTicks += b-a;
    }
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << ", " << name<P2>() << "): "
        << float(perfIndex/PERF_ITERATIONS) << " absolute ticks: " << totalTicks / float(PERF_ITERATIONS * TEST_PER_ITERATION) << std::endl;*/

    std::vector<uint64> iterationTimes(PERF_ITERATIONS);
    for(int t = 0; t < PERF_ITERATIONS; ++t)
    {
        // Fill data set for x intersections
        // Since the Fast... types have no default ctor I use the dummyMem+cast here.
        char dummyMem[512];
        std::vector<P0> geo0(TEST_PER_ITERATION, *reinterpret_cast<P0*>(dummyMem));
        std::vector<P1> geo1(TEST_PER_ITERATION, *reinterpret_cast<P1*>(dummyMem));
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
        for(int i = 0; i < TEST_PER_ITERATION; ++i)
            for(int j = 0; j < TEST_PER_ITERATION; ++j)
            {
                res = _func(geo0[i], geo1[j], rt);
                eat = rt;
            }
        uint64 b = ticks();
        iterationTimes[t] = b-a;
    }

    std::sort(iterationTimes.begin(), iterationTimes.end());
    std::cerr << "Performance " << _funcName << '(' << name<P0>() << ", " << name<P1>() << ", " << name<P2>() << "):\t"
        << deltaTicksToMilliSeconds(ei::uint64(iterationTimes[PERF_ITERATIONS/2] * 1000000.0 / (TEST_PER_ITERATION * TEST_PER_ITERATION)))
        << " ms/M\n";
}

