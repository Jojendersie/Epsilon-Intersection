#pragma once

/// \brief Include or uncomment lines to customize namespaces and similar.
///


/// \brief All types from elementarytypes.hpp are without namespace.
/// \details The default is on.
//#define USE_ELEMENTARIES_WITHOUT_NAMESPACE

/// \brief Using this option changes the namespace ei:: to ε::.
/// \details Do not use this if the symbol epsilon above is not shown correctly,
///    or enable utf8 in all your files.
///
///    The default is off.
#define USE_UTF8_NAMESPACES


/// \brief Enable a bunch of debug assertions in the vector code.
/// \details As usual the assertions are not contained in release mode. This
///    macro lets you further disable the assertions in debug mode to speed
///    up the application.
///
///    Level 0: Deactivate all assertions. Not recommended.
///    Level 1: Use assertions where errors are likely. Default option.
///    Level 2: Maximum. Permanent tests everywhere.
#define ASSERTION_LEVEL		2

#include "details/config.inl"
