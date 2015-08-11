#define NOMINMAX
#include <Windows.h>

#include "ei/elementarytypes.hpp"

/// Also test name clashes in this file: independent of the
/// USE_ELEMENTARIES_WITHOUT_NAMESPACE option it should compile.
using namespace ei;

ei::uint64 ticks()
{
    uint64 tickCount;
    QueryPerformanceCounter((LARGE_INTEGER*)&tickCount);
    return tickCount;
}

void eatMyDummy(float _dummy)
{
    volatile float devour = _dummy;
}