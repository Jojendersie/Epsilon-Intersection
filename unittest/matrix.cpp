#include "gam/matrix.hpp"
#include "unittest.hpp"

#include <iostream>

using namespace gam;
using namespace std;

bool test_matrix()
{
	bool result = true;

	// Test type size
	TEST( sizeof(Matrix<byte, 1, 2>) == 2,
		"sizeof(Matrix<byte, 1, 2>) is " << sizeof(Matrix<byte, 1, 2>) << " instead of 2\n"
	);
	TEST( sizeof(Matrix<double, 2, 1>) == 16,
		"sizeof(Matrix<double, 2, 1>) is " << sizeof(Matrix<double, 2, 1>) << " instead of 16\n"
	);
	TEST( sizeof(Matrix<float, 1, 3>) == 12,
		"sizeof(Matrix<float, 1, 3>) is " << sizeof(Matrix<float, 1, 3>) << " instead of 12\n"
	);
	TEST( sizeof(Matrix<byte, 3, 1>) == 3,
		"sizeof(Matrix<byte, 3, 1>) is " << sizeof(Matrix<byte, 3, 1>) << " instead of 3\n"
	);
	TEST( sizeof(Matrix<float, 4, 1>) == 16,
		"sizeof(Matrix<float, 4, 1>) is " << sizeof(Matrix<float, 4, 1>) << " instead of 16\n"
	);
	TEST( sizeof(Matrix<int16, 1, 4>) == 8,
		"sizeof(Matrix<int16, 1, 4>) is " << sizeof(Matrix<int16, 1, 4>) << " instead of 8\n"
	);
	TEST( sizeof(Matrix<uint16, 2, 2>) == 8,
		"sizeof(Matrix<uint16, 2, 2>) is " << sizeof(Matrix<uint16, 2, 2>) << " instead of 8\n"
	);
	TEST( sizeof(Matrix<byte, 3, 3>) == 9,
		"sizeof(Matrix<byte, 3, 3>) is " << sizeof(Matrix<byte, 3, 3>) << " instead of 9\n"
	);
	TEST( sizeof(Matrix<float, 4, 4>) == 64,
		"sizeof(Matrix<float, 4, 4>) is " << sizeof(Matrix<float, 4, 4>) << " instead of 64\n"
	);
	TEST( sizeof(Matrix<float, 10, 10>) == 400,
		"sizeof(Matrix<float, 10, 10>) is " << sizeof(Matrix<float, 10, 10>) << " instead of 400\n"
	);

	// When included the following line must cause a static assertion.
	// Matrix<int, 1, 1> errorMatrix;

	return result;
}