Epsilon Intersection Library - ε
===============================================================================

1. What you can do!
-------------------------------------------------------------------------------

2. How to use the library?
-------------------------------------------------------------------------------

There are 2 different ways to include γ into your project:

  1. Add all files to your current project and compile as usual.
  2. Compile a library and use that
  
All interfaces are declared in "include/ei/<xyz>.hpp" files. You never need to
include or look into files from the "include/ei/details" or "src" directory.

The hierarchy of header-files is as follows.
config --- elementarytypes --- matrix -|- 2dtypes --- 2dintersection
                                       -- 3dtypes --- 3dintersection
Hence each one includes all its parents. You always need to include only the
latest one.