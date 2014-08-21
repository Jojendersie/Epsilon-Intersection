#include "gam/elementarytypes.hpp"
#include "unittest.hpp"

#include <iostream>
#include <limits>

using namespace gam;
using namespace std;

bool test_elementaries()
{
    bool result = true;

    TEST( static_cast<uint>(-1) == static_cast<unsigned int>(-1),
		"The type 'uint' is wrong defined. It should be equal to 'unsigned int'\n"
	);

	TEST( sizeof(uint8) == 1 && sizeof(int8) == 1 && sizeof(byte) == 1,
		"The 8bit integer types 'uint8', 'int8' and 'byte' are not one byte large\n"
	);

	TEST( sizeof(uint16) == 2 && sizeof(int16) == 2,
		"The 16-bit integer types 'uint16' and 'int16' are not two byte large\n"
	);

	TEST( sizeof(uint32) == 4 && sizeof(int32) == 4,
		"The 32-bit integer types 'uint32' and 'int32' are not four byte large\n"
	);	

	TEST( sizeof(uint64) == 8 && sizeof(int64) == 8,
		"The 64-bit integer types 'uint64' and 'int64' are not eight byte large\n"
	);

	TEST( std::numeric_limits<uint8>::is_integer,
		"'uint8' is not an integer\n"
	);
	TEST( std::numeric_limits<int8>::is_integer,
		"'int8' is not an integer\n"
	);
	TEST( std::numeric_limits<byte>::is_integer,
		"'byte' is not an integer\n"
	);
	TEST( std::numeric_limits<uint16>::is_integer,
		"'uint16' is not an integer\n"
	);
	TEST( std::numeric_limits<int16>::is_integer,
		"'int16' is not an integer\n"
	);
	TEST( std::numeric_limits<uint32>::is_integer,
		"'uint32' is not an integer\n"
	);
	TEST( std::numeric_limits<int32>::is_integer,
		"'int32' is not an integer\n"
	);
	TEST( std::numeric_limits<uint64>::is_integer,
		"'uint64' is not an integer\n"
	);
	TEST( std::numeric_limits<int64>::is_integer,
		"'int64' is not an integer\n"
	);

	return result;
}