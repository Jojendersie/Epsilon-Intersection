#include "ei/elementarytypes.hpp"

namespace ei {

    // ********************************************************************* //
    float successor(float _number)
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == INF) return _number;
        if(_number == -INF) return details::ReinterpretFloat(0xff7fffffu).f;

        details::ReinterpretFloat bits = _number;
        unsigned sign = bits.i & 0x80000000u;
        unsigned mantissa = bits.i & 0x7fffffffu;
        if( sign && mantissa )
            return details::ReinterpretFloat( sign | (mantissa - 1) ).f;
        else return details::ReinterpretFloat( mantissa + 1 ).f;
    }

    double successor(double _number)
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == INF_D) return _number;
        if(_number == -INF_D) return details::ReinterpretDouble(0xffeffffffffffffful).f;

        details::ReinterpretDouble bits = _number;
        uint64 sign = bits.i & 0x8000000000000000ul;
        uint64 mantissa = bits.i & 0x7ffffffffffffffful;
        if( sign && mantissa )
            return details::ReinterpretDouble( sign | (mantissa - 1) ).f;
        else return details::ReinterpretDouble( mantissa + 1 ).f;
    }

    // ********************************************************************* //
    float predecessor(float _number)
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == -INF) return _number;
        if(_number == INF) return details::ReinterpretFloat(0x7f7fffffu).f;

        details::ReinterpretFloat bits = _number;
        unsigned sign = bits.i & 0x80000000u;
        unsigned mantissa = bits.i & 0x7fffffffu;
        if( sign || !mantissa )
            return details::ReinterpretFloat( 0x80000000u | (mantissa + 1) ).f;
        else return details::ReinterpretFloat( mantissa - 1 ).f;
    }

    double predecessor(double _number)
    {
        // Handling NaN
        if(_number != _number) return _number;
        // Infinity
        if(_number == -INF) return _number;
        if(_number == INF) return details::ReinterpretDouble(0x7feffffffffffffful).f;

        details::ReinterpretDouble bits = _number;
        uint64 sign = bits.i & 0x8000000000000000ul;
        uint64 mantissa = bits.i & 0x7ffffffffffffffful;
        if( sign || !mantissa )
            return details::ReinterpretDouble( 0x8000000000000000ul | (mantissa + 1) ).f;
        else return details::ReinterpretDouble( mantissa - 1 ).f;
    }

}