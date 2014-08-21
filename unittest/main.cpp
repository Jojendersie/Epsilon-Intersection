#include <iostream>
using namespace std;

bool test_elementaries();

int main()
{
    if( test_elementaries() )
		cerr << "Successfully completed: Elementary types." << std::endl;

	return 0;
}
