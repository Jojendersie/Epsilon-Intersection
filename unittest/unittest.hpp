#pragma once

#define TEST(expression, message) \
	if( !(expression) )           \
	{                             \
		cerr << message << '\n';  \
		result = false;           \
	}
