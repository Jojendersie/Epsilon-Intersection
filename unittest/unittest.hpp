#pragma once

#define TEST(expression, message) \
	if( !(expression) )           \
	{                             \
		cerr << message;          \
		result = false;           \
	}
