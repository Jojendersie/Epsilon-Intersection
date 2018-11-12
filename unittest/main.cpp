#include <iostream>
using namespace std;

bool test_elementaries();
bool test_matrix();
bool test_quaternion();
bool test_2dtypes();
bool test_2dintersections();
bool test_3dtypes();
bool test_3dintersections();
bool test_stdextensions();
bool test_primes();
bool test_conversions();

int main()
{
    if(  test_conversions() )
        cerr << "Successfully completed: Conversions." << std::endl;

    if( test_elementaries() )
        cerr << "Successfully completed: Elementary types." << std::endl;

    if( test_stdextensions() )
        cerr << "Successfully completed: std:extensions." << std::endl;

    if( test_matrix() )
        cerr << "Successfully completed: Matrix type." << std::endl;

    if( test_quaternion() )
        cerr << "Successfully completed: Quaternion type." << std::endl;

    if( test_2dtypes() )
        cerr << "Successfully completed: 2D types test." << std::endl;

    if( test_2dintersections() )
        cerr << "Successfully completed: 2D intersection test." << std::endl;

    if( test_3dtypes() )
        cerr << "Successfully completed: 3D types test." << std::endl;

    if( test_3dintersections() )
        cerr << "Successfully completed: 3D intersection test." << std::endl;

    if( test_primes() )
        cerr << "Successfully completed: Primes." << std::endl;

    return 0;
}
