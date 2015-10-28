#include "ei/stdextensions.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace ei;
using namespace std;

bool test_stdextensions()
{
    IVec2 v0(6, 1);
    Vec3 v1(0.1f, 0.5f, -10.0f);
    Vec4 v2(1.0f);
    Mat3x3 m0 = rotation(v1);
    
    // Compilation of hash functions
    {
        hash<IVec2> hv2;
        cout << "Hash of IVec2(" << ") is " << hv2(v0);
        hash<Vec3> hv3;
        cout << "Hash of Vec3(" << ") is " << hv3(v1);
        hash<Vec4> hv4;
        cout << "Hash of Vec4(" << ") is " << hv4(v2);
        hash<Mat3x3> hv33;
        cout << "Hash of Mat3x3(" << ") is " << hv33(m0);
    }

    return true;
}