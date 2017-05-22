#pragma once

#define TEST(expression, message)      \
    if( !(expression) )                \
    {                                  \
        std::cerr << message << '\n';  \
        result = false;                \
    }

ei::uint64 ticks();
float deltaTicksToMicroSeconds(ei::uint64 _deltaTicks);
float deltaTicksToMilliSeconds(ei::uint64 _deltaTicks);
float deltaTicksToSeconds(ei::uint64 _deltaTicks);
/// \brief A function to avoid too heavy code optimization (code cannot be deleted)
void eatMyDummy(float _dummy);
void eatMyDummy(bool _dummy);


/// \brief Returns a number in [0,1]
inline float rnd()
{
    static ei::uint64 s_rndState = 47;
    s_rndState ^= s_rndState >> 12;
    s_rndState ^= s_rndState << 25;
    s_rndState ^= s_rndState >> 27;
    s_rndState *= 2685821657736338717ULL;
    return static_cast<float>(s_rndState % 8388593) / 8388592.0f;
}
