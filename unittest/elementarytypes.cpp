#include "elementarytypes.hpp"

bool test_elementaries()
{
    bool result = true;
    if( static_cast<γ::uint>(-1) != static_cast<unsigned int>(-1) )
    {
        cerr << "The type 'uint' is wrong defined. It should be equal to 'unsigned int'\n";
		result = false;
    }
    
    return result;
}