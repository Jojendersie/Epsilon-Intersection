#define NOMINMAX
#include <Windows.h>

#include "ei/elementarytypes.hpp"

/// Also test name clashes in this file: independent of the
/// EI_GLOBAL_ELEMENTARIES option it should compile.
using namespace ei;

ei::uint64 ticks()
{
    uint64 tickCount;
    QueryPerformanceCounter((LARGE_INTEGER*)&tickCount);
    return tickCount;
}

float deltaTicksToMicroSeconds(ei::uint64 _deltaTicks)
{
    uint64 frequency;
    QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);
    return float(_deltaTicks * 1000000.0 / frequency);
}

float deltaTicksToMilliSeconds(ei::uint64 _deltaTicks)
{
    uint64 frequency;
    QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);
    return float(_deltaTicks * 1000.0 / frequency);
}

float deltaTicksToSeconds(ei::uint64 _deltaTicks)
{
    uint64 frequency;
    QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);
    return float(_deltaTicks / double(frequency));
}

#pragma warning(push)
#pragma warning(disable:4189)
void eatMyDummy(float _dummy)
{
    volatile float devour = _dummy;
}

void eatMyDummy(bool _dummy)
{
    volatile bool devour = _dummy;
}
#pragma warning(pop)