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
        Vec3 sRgb0{0.735356927f, 0.735356927f, 0.735356927f};
        Vec3 sRgb1{0.5f, 0.5f, 0.5f};
        Vec3 rgb1{0.214041144f, 0.214041144f, 0.214041144f};

        TEST(approx(rgb0, sRgbToRgb(sRgb0)), "Conversion sRGB->RGB not correct.");
        TEST(approx(rgb0, sRgbToRgb(rgbToSRgb(rgb0))), "Conversion RGB->sRGB->RGB not reciprocal.");
        TEST(approx(sRgb1, rgbToSRgb(rgb1)), "Conversion RGB->sRGB not correct.");
        TEST(approx(sRgb1, rgbToSRgb(sRgbToRgb(sRgb1))), "Conversion sRGB->RGB->sRGB not reciprocal.");
    }

    { // XYZ <-> RGB
        Vec3 rgb0{0.5f, 0.5f, 0.5f};
        Vec3 xyz0 = rgbToXyz(rgb0);
        Vec3 xyz1{0.5f, 0.5f, 0.5f};
        Vec3 rgb1 = xyzToRgb(xyz1);

        TEST(approx(rgb0, xyzToRgb(xyz0)), "Conversion RGB->XYZ->RGB not reciprocal.");
        TEST(approx(xyz1, rgbToXyz(rgb1)), "Conversion XYZ->RGB->XYZ not reciprocal.");
    }


    { // HSV
        Vec3 brownRGB{0.36f, 0.18f, 0.09f};
        Vec3 brownHSV{0.0555555555f, 0.75f, 0.36f};
        Vec3 magentaRGB{1.0f, 0.0f, 1.0f};
        Vec3 magentaHSV{0.833333333f, 1.0f, 1.0f};
        Vec3 hdrgreenRGB{0.0f, 12.0f, 0.0f};
        Vec3 hdrgreenHSV{0.333333333f, 1.0f, 12.0f};
        TEST(approx(rgbToHsv(brownRGB), brownHSV), "Conversion RGB->HSV (brown) failed.");
        TEST(approx(hsvToRgb(brownHSV), brownRGB), "Conversion HSV->RGB (brown) failed.");
        TEST(approx(rgbToHsv(magentaRGB), magentaHSV), "Conversion RGB->HSV (magenta) failed.");
        TEST(approx(hsvToRgb(magentaHSV), magentaRGB), "Conversion HSV->RGB (magenta) failed.");
        TEST(approx(rgbToHsv(hdrgreenRGB), hdrgreenHSV), "Conversion RGB->HSV (green) failed.");
        TEST(approx(hsvToRgb(hdrgreenHSV), hdrgreenRGB), "Conversion HSV->RGB (green) failed.");
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

    { // L*a*b*
        const Vec3 xyz0 = rgbToXyz(Vec3{ 1.0f, 1.0f, 1.0f });
        const Vec3 xyz1 = rgbToXyz(Vec3{ 0.0f, 0.0f, 0.0f });
        const Vec3 xyz2 = rgbToXyz(Vec3{ 0.0f, 0.5f, 0.0f });
        const Vec3 xyz3 = rgbToXyz(Vec3{ 3.0f, 3.0f, 0.0f });
        const Vec3 lab0 = xyzToLab(xyz0);
        const Vec3 lab1 = xyzToLab(xyz1);
        const Vec3 lab2 = xyzToLab(xyz2);
        const Vec3 lab3 = xyzToLab(xyz3);
        TEST(approx(lab0, Vec3{1.0f, 0.0f, 0.0f}), "L*a*b* white wrong.");
        TEST(approx(lab1, Vec3{0.0f, 0.0f, 0.0f}), "L*a*b* black wrong.");
        TEST(approx(labToXyz(lab2), xyz2), "XYZ -> L*a*b* -> XYZ for xyz2 wrong.");
        TEST(approx(labToXyz(lab3), xyz3), "XYZ -> L*a*b* -> XYZ for xyz3 wrong.");
    }

    { // Octahedral
        Vec3 v0 = normalize(Vec3{ 1.0f, 0.0f, 0.0f });
        Vec3 v1 = normalize(Vec3{ 0.0f, 1.0f, 0.0f });
        Vec3 v2 = normalize(Vec3{ 0.0f, 0.0f, 1.0f });
        Vec3 v3 = normalize(Vec3{ -1.0f, 0.0f, 0.0f });
        Vec3 v4 = normalize(Vec3{ 0.0f, -1.0f, 0.0f });
        Vec3 v5 = normalize(Vec3{ 0.0f, 0.0f, -1.0f });
        Vec3 v6 = normalize(Vec3{ 1.0f, 1.0f, 1.0f });
        Vec3 v7 = normalize(Vec3{ -1.0f, -1.0f, -1.0f });
        Vec3 v8 = normalize(Vec3{ 0.5f, 0.7f, -0.4f });
        Vec3 v9 = Vec3 {0.0f, 0.0f, 0.0f};
        uint32 oct0 = packOctahedral32(v0);
        uint32 oct1 = packOctahedral32(v1);
        uint32 oct2 = packOctahedral32(v2);
        uint32 oct3 = packOctahedral32(v3);
        uint32 oct4 = packOctahedral32(v4);
        uint32 oct5 = packOctahedral32(v5);
        uint32 oct6 = packOctahedral32(v6);
        uint32 oct7 = packOctahedral32(v7);
        uint32 oct8 = packOctahedral32(v8);
        uint32 oct9 = packOctahedral32(v9);
        TEST(approx(v0, unpackOctahedral32(oct0)), "Octahedral packing of v0 invalid.");
        TEST(approx(v1, unpackOctahedral32(oct1)), "Octahedral packing of v1 invalid.");
        TEST(approx(v2, unpackOctahedral32(oct2)), "Octahedral packing of v2 invalid.");
        TEST(approx(v3, unpackOctahedral32(oct3)), "Octahedral packing of v3 invalid.");
        TEST(approx(v4, unpackOctahedral32(oct4)), "Octahedral packing of v4 invalid.");
        TEST(approx(v5, unpackOctahedral32(oct5)), "Octahedral packing of v5 invalid.");
        TEST(approx(v6, unpackOctahedral32(oct6), 1e-4f), "Octahedral packing of v6 invalid.");
        TEST(approx(v7, unpackOctahedral32(oct7), 1e-4f), "Octahedral packing of v7 invalid.");
        TEST(approx(v8, unpackOctahedral32(oct8), 1e-4f), "Octahedral packing of v8 invalid.");
        TEST(approx(v9, unpackOctahedral32(oct9)), "Octahedral packing of v9 invalid.");
    }

    { // Tangent space
        OrthoSpace o0{Quaternion{normalize(Vec3(1.0f, -0.5f, 0.1f)), 0.234f}};
        OrthoSpace o1{Quaternion{Vec3(1.0f, 0.0f, 0.0f), PI}};
        OrthoSpace o2{qidentity()};
        OrthoSpace o3{diag(Vec3(-1.0f, 1.0f, 1.0f))};
        uint64 op0 = packOrthoSpace64(o0);
        uint64 op1 = packOrthoSpace64(o1);
        uint64 op2 = packOrthoSpace64(o2);
        uint64 op3 = packOrthoSpace64(o3);
        TEST(approx(o0, unpackOrthoSpace64(op0), 0.003f), "OrthoSpace o0 packing to 64 bit failed.");
        TEST(approx(o1, unpackOrthoSpace64(op1)), "OrthoSpace o1 packing to 64 bit failed.");
        TEST(approx(o2, unpackOrthoSpace64(op2)), "OrthoSpace o2 packing to 64 bit failed.");
        TEST(approx(o3, unpackOrthoSpace64(op3)), "OrthoSpace o3 packing to 64 bit failed.");
    }

    return result;
}