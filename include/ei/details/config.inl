#pragma once

// Just show the namespace once to make sure it exists when aliasing.
namespace ei {}
#ifdef USE_UTF8_NAMESPACES
    namespace ε = ei;
#endif
// An easter egg for german speaking people ;-).
namespace egg = ei;


// Define the assertions depending on the option level.
#if defined(DEBUG) || defined(_DEBUG) || !defined(NDEBUG)
#   if defined(_WIN32)
#       define DEBUG_BREAK __debugbreak()
#   else
#       define DEBUG_BREAK ::kill(0, SIGTRAP)
#   endif
	
	/// \brief Default assert macro for all our needs. Use this instead of <cassert>
	/// \details This form of assertion breaks at the very point of code
	///    where you placed it. Further you are able to step over the assertion
	///    in any condition.
	///
	///    The do while loop is the only way to guarantee robustness of the macro.
	///    It forces you to set a semicolon and can be used within if, for, ...
	///    without curled brackets:
	///        if(...)
	///            assertlvl1(..., "Even if compiled out the if statement does not
	///                             infer with the next line.");
	///        somestatement;
#   define assertlvl1(condition, errorMessage) \
        do {                                   \
            if((condition) == false) {         \
	            DEBUG_BREAK;                   \
	        }                                  \
	    } while((void)0,0)

	/// \brief Optional assert macro for highly frequented code.
	/// \details \see assertlvl1()
#   if ASSERTION_LEVEL > 1
#       define assertlvl2(condition, errorMessage) \
            do {                                   \
                if((condition) == false) {         \
                    DEBUG_BREAK;                   \
                }                                  \
            } while((void)0,0)
#   endif

#else

	/// \brief Empty replacement if switched off.
	/// \details \see assertlvl1()
#   define assertlvl1(condition, errorMessage) do { } while((void)0,0)
#   if ASSERTION_LEVEL > 1
#       define assertlvl2(condition, errorMessage) do { } while((void)0,0)
#   endif

#   undef DEBUG_BREAK
#endif
