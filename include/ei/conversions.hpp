#pragma once

#include "vector.hpp"

/// \brief Utility functions for conversions of colors and vector representations.

namespace ei {

    // Conversion of sRGB to linear RGB (single component)
    inline float sRgbToRgb(float _c)
    {
        if(_c <= 0.0031308f)
            return _c * 12.92f;
        return 1.055f * pow(_c, 1.0f/2.4f) - 0.055f;;
    }

    // Conversion of sRGB to linear RGB (two components)
    inline Vec2 sRgbToRgb(const Vec2 & _sRg)
    {
        return Vec2{sRgbToRgb(_sRg.r), sRgbToRgb(_sRg.g)};
    }

    // Conversion of sRGB to linear RGB
    inline Vec3 sRgbToRgb(const Vec3 & _sRgb)
    {
        return Vec3{sRgbToRgb(_sRgb.r), sRgbToRgb(_sRgb.g), sRgbToRgb(_sRgb.b)};
    }

    // Conversion of linear RGB to sRGB (single component)
    inline float rgbToSRgb(float _c)
    {
        if(_c <= 0.04045f)
            return _c / 12.92f;
        return pow((_c + 0.055f) / 1.055f, 2.4f);
    }

    // Conversion of linear RGB to sRGB (two components)
    inline Vec2 rgbToSRgb(const Vec2 & _rg)
    {
        return Vec2{rgbToSRgb(_rg.r), rgbToSRgb(_rg.g)};
    }

    // Conversion of linear RGB to sRGB
    inline Vec3 rgbToSRgb(const Vec3 & _rgb)
    {
        return Vec3{rgbToSRgb(_rgb.r), rgbToSRgb(_rgb.g), rgbToSRgb(_rgb.b)};
    }


    // Conversion of CIE XYZ to linear RGB
    inline Vec3 xyzToRgb(const Vec3 & _xyz)
    {
        static constexpr Mat3x3 XYZ_TO_RGB {
             3.2406f, -1.5372f, -0.4986f,
            -0.9689f,  1.8758f,  0.0415f,
             0.0557f, -0.2040f,  1.0570f
        };
        return XYZ_TO_RGB * _xyz;
    }

    // Conversion of linear RGB to CIE XYZ
    inline Vec3 rgbToXyz(const Vec3 & _rgb)
    {
        static constexpr Mat3x3 RGB_TO_XYZ {
            0.4124f, 0.3576f, 0.1805f,
            0.2126f, 0.7152f, 0.0722f,
            0.0193f, 0.1192f, 0.9505f
        };
        return RGB_TO_XYZ * _rgb;
    }


    // Conversion of YCgCo color space to linear RGB
    inline Vec3 yCgCoToRgb(const Vec3 & _yCgCo)
    {
        float tmp = _yCgCo.x - _yCgCo.y;
        return Vec3(tmp + _yCgCo.z, _yCgCo.x + _yCgCo.y, tmp - _yCgCo.z);
    }

    // Conversion of linear RGB to YCgCo color space
    inline Vec3 rgbToYCgCo(const Vec3 & _rgb)
    {
        return Vec3{(_rgb.r + _rgb.b) * 0.25f + _rgb.g * 0.5f,
                    (_rgb.r + _rgb.b) * -0.25f + _rgb.g * 0.5f,
                    (_rgb.r - _rgb.b) * 0.5f};
    }




    // Discretize a [0,1]^3 vector into a single integer with 11.11.10 bits
    // for the components.
    // Note: This is not the R11G11B10F (float) format from textures!
    inline uint32 packR11G11B10(Vec3 _v)
    {
        eiAssertWeak(all(greatereq(_v, 0.0f)) && all(lesseq(_v, 1.0f)), "Unclamped color cannot be converted into R11G11B10 format!");
        _v *= Vec3{2047.0, 2047.0, 1023.0};
        return (uint32(_v.r) << 21) | (uint32(_v.g) << 10) | uint32(_v.b);
    }

    // Unpack a 11.11.10 bit descretized vector into a full Vec3
    // Note: This is not the R11G11B10F (float) format from textures!
    inline Vec3 unpackR11G11B10(uint32 _code)
    {
        return Vec3 {
            float(_code >> 21) / 2047.0,
            float((_code >> 10) & 0x7ff) / 2047.0,
            float(_code & 0x3ff) / 1023.0
        };
    }


    // Pack an HDR color value into RGB9E5 shared exponent format
    // See https://www.khronos.org/registry/OpenGL/extensions/EXT/EXT_texture_shared_exponent.txt
    inline uint32 packRGB9E5(const Vec3 & _v)
    {
        eiAssertWeak(all(greatereq(_v, 0.0f)), "Vector must be positive to be packed into RGB9E5");
        float maxComp = ei::max(_v);
        if(maxComp <= 6e-8f) return 0;
        // Code from spec:
        //int e = ei::max(0, ei::floor(log2f(maxComp)) + 1); // Without bias
        //float norm = exp2f(-(e - 9));	// - #MantissaBits
        //if(maxComp * norm >= 511.5f) {	// <=> floor(maxComp * norm + 0.5f) == 512
        //    ++e;
        //    norm *= 0.5f;
        //}
        // Equivalent, faster code?:
        int e;
        float norm = frexpf(maxComp, &e) * 511.5f / maxComp;
        e += 15;
        if(e < 0) { norm *= exp2f(e); e = 0; }	// Denormalized values
        return floor(_v.r * norm + 0.5f)
            | (floor(_v.g * norm + 0.5f) << 9)
            | (floor(_v.b * norm + 0.5f) << 18)
            | (ei::min(e, 31) << 27); // Clamp the number if value range of RGB9E5 is exceeded
    }

    // Unpack a RGB9E5 value into a full Vec3
    inline Vec3 unpackRGB9E5(uint32 _code)
    {
        float e = exp2f((_code>>27) - 15 - 9); // - ExponentBias - #MantissaBits
        return Vec3 {
            (_code & 0x1ff) * e,
            ((_code>>9) & 0x1ff) * e,
            ((_code>>18) & 0x1ff) * e
        };
    }

    // Pack an HDR color value into RGB8E8 shared exponent format.
    // This is used for example by the *.hdr file format
    inline uint32 packRGB8E8(const Vec3 & _v)
    {
        eiAssertWeak(all(greatereq(_v, 0.0f)), "Vector must be positive to be packed into RGB8E8");
        float maxComp = ei::max(_v);
        if(maxComp <= 2e-41f) return 0;
        int e;
        float norm = frexpf(maxComp, &e) * 255.5f / maxComp;
        e += 128;
        if(e < 0) { norm *= exp2f(e); e = 0; }	// Denormalized values
        return floor(_v.r * norm + 0.5f)
            | (floor(_v.g * norm + 0.5f) << 8)
            | (floor(_v.b * norm + 0.5f) << 16)
            | (ei::min(e, 255) << 24); // Clamp the number if value range of RGB8E8 is exceeded
    }

    // Unpack a RGB8E8 value into a full Vec3
    inline Vec3 unpackRGB8E8(uint32 _code)
    {
        float e = exp2f((_code>>24) - 128 - 8); // - ExponentBias - #MantissaBits
        return Vec3 {
            (_code & 0xff) * e,
            ((_code>>8) & 0xff) * e,
            ((_code>>16) & 0xff) * e
        };
    }

    // Pack a direction vector into a 32 bit integer using octahedral mapping
    inline uint32 packOctahedral32(const Vec3 & _d)
    {
        eiAssertWeak(approx(len(_d), 1.0f), "Can only pack direction vectors using octahedral mapping");
        float l1norm = abs(_d.x) + abs(_d.y) + abs(_d.z);
        float u,v;
        if(_d.z >= 0) {
            u = _d.x / l1norm;
            v = _d.y / l1norm;
        } else { // warp lower hemisphere
            u = (1 - _d.y / l1norm) * (_d.x >= 0 ? 1 : -1);
            v = (1 - _d.x / l1norm) * (_d.y >= 0 ? 1 : -1);
        }
        return floor((u / 2.0f + 0.5f) * 65535.49f + 0.5f)
            | (floor((v / 2.0f + 0.5f) * 65535.49f + 0.5f) << 16);
    }

    // Unpack a direction vector from octahedral mapping
    inline Vec3 unpackOctahedral32(uint32 _code)
    {
        float u = (_code & 0xff) / 65535.0f;
        float v = (_code >> 16) / 65535.0f;
        u = u * 2 - 1;
        v = v * 2 - 1;
        float x, y, z = 1 - abs(u) - abs(v);
        if(z >= 0) {
            x = u;
            y = v;
        } else {
            x = (1 - abs(v)) * (u >= 0 ? 1 : -1);
            y = (1 - abs(u)) * (v >= 0 ? 1 : -1);
        }
        return normalize(Vec3{x,y,z})
    }

} // namespace ei