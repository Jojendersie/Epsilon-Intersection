#include <iostream>
using namespace std;

bool test_elementaries();
bool test_matrix();
bool test_2dtypes();

int main()
{
    if( test_elementaries() )
        cerr << "Successfully completed: Elementary types." << std::endl;

    if( test_matrix() )
        cerr << "Successfully completed: Matrix type." << std::endl;

    if( test_2dtypes() )
        cerr << "Successfully completed: 2D types test." << std::endl;

    return 0;
}
