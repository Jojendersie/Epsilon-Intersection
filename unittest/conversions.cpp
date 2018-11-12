#include "ei/conversions.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_conversions()
{
    bool result = true;

    { // RGB <-> sRGB
        Vec3 rgb0{0.5f, 0.5f, 0.5f};
        Vec3 sRgb0 = rgbToSRgb(rgb0);
        Vec3 sRgb1{0.5f, 0.5f, 0.5f};
        Vec3 rgb1 = sRgbToRgb(sRgb1);

        TEST(approx(rgb0, sRgbToRgb(sRgb0)), "Conversion RGB->sRGB->RGB not reciprocal.");
        TEST(approx(sRgb1, rgbToSRgb(rgb1)), "Conversion sRGB->RGB->sRGB not reciprocal.");
    }

    { // XYZ <-> RGB
        Vec3 rgb0{0.5f, 0.5f, 0.5f};
        Vec3 xyz0 = rgbToXyz(rgb0);
        Vec3 xyz1{0.5f, 0.5f, 0.5f};
        Vec3 rgb1 = xyzToRgb(xyz1);

        TEST(approx(rgb0, xyzToRgb(xyz0)), "Conversion RGB->XYZ->RGB not reciprocal.");
        TEST(approx(xyz1, rgbToXyz(rgb1)), "Conversion XYZ->RGB->XYZ not reciprocal.");
    }

    { // Vec3 <-> R11G11B10
        Vec3 rgb0{0.5f, 0.5f, 0.5f};
        uint32 pack0 = packR11G11B10(rgb0);
        uint32 pack1 = 0xffffffff;
        Vec3 rgb1 = unpackR11G11B10(pack1);

        TEST(approx(rgb0, unpackR11G11B10(pack0), 1e-3f), "Conversion Vec3->R11G11B10->Vec3 not reciprocal.");
        TEST(pack1 == packR11G11B10(rgb1), "Conversion R11G11B10->Vec3->R11G11B10 not reciprocal.");
    }

    { // Vec3 <-> RGB9E5
        Vec3 rgb0{0.5f, 0.5f, 0.5f};
        uint32 pack0 = packRGB9E5(rgb0);
        uint32 pack1 = 0xffffffff;
        Vec3 rgb1 = unpackRGB9E5(pack1);

        TEST(approx(rgb0, unpackRGB9E5(pack0)), "Conversion Vec3->RGB9E5->Vec3 not reciprocal.");
        TEST(pack1 == packRGB9E5(rgb1), "Conversion RGB9E5->Vec3->RGB9E5 not reciprocal.");
    }

    return result;
}