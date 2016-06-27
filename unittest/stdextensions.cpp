#include "ei/stdextensions.hpp"
#include "unittest.hpp"

#include <iostream>
#include <functional>

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
        cout << "Hash of IVec2(..." << ") is " << hv2(v0) << '\n';
        hash<Vec3> hv3;
        cout << "Hash of Vec3(..." << ") is " << hv3(v1) << '\n';
        hash<Vec4> hv4;
        cout << "Hash of Vec4(..." << ") is " << hv4(v2) << '\n';
        //hash<Matrix<int,1,10>> hv11;
        hash<Mat3x3> hv33;
        cout << "Hash of Mat3x3(..." << ") is " << hv33(m0) << '\n';
    }

    // Comparison objects
    {
        equal_to<IVec2> ev2;
        cout << "equal_to of IVec2(..." << ") is " << ev2(v0, v0) << '\n';
        equal_to<Mat3x3> ev33;
        cout << "equal_to of Mat3x3(..." << ") is " << ev33(m0, m0) << '\n';
    }

    // Pretty printers
    {
        IMat2x2 m1(1, 210, 7626, 73);
        cout << "showing pretty printers\n";
        cout << v0 << '\n';
        cout << v1 << '\n';
        cout << v2 << "\n\n";
        cout << m0 << "\n\n";
        cout << m1 << "\n\n";
    }

    return true;
}