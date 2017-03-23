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

#pragma warning(push)
#pragma warning(disable:4189)
void eatMyDummy(float _dummy)
{
    volatile float devour = _dummy;
}
#pragma warning(pop)