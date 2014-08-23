#pragma once

/// \brief Include or uncomment lines to customize namespaces and similar.
///


/// \brief All types from elementarytypes.hpp are without namespace.
/// \details The default is on.
//#define USE_ELEMENTARIES_WITHOUT_NAMESPACE

/// \brief Using this option changes the namespace gam:: to γ::.
/// \details Do not use this if the symbol gamma above is not shown correctly,
///    or enable utf8 in all your files.
///
///    The default is off.
#define USE_UTF8_NAMESPACES

#include "details/config.inl"
