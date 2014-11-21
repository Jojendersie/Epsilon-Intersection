#include "ei/elementarytypes.hpp"

#include <Windows.h>

ε::uint64 ticks()
{
    ε::uint64 tickCount;
    QueryPerformanceCounter((LARGE_INTEGER*)&tickCount);
    return tickCount;
}

void eatMyDummy(float _dummy)
{
    volatile float devour = _dummy;
}