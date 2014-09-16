#include <iostream>
using namespace std;

bool test_elementaries();
bool test_matrix();

int main()
{
    if( test_elementaries() )
        cerr << "Successfully completed: Elementary types." << std::endl;

    if( test_matrix() )
        cerr << "Successfully completed: Matrix type." << std::endl;

    return 0;
}
