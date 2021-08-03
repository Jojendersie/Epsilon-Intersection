Epsilon Intersection Library - ε
===============================================================================

ε is a super flexible and simple to use **C++ library** designed for games. It is computing **intersections** of various three dimensional bodies and containes different utilities you could need for culling, ray casting or similar functions.

It provides powerful template based matrix/vector and quaternion classes which come with a ton of useful functions.
Due to its lean design you can easily integrate the parts you need into your project, just by adding the source!


How to use the library?
-------------------------------------------------------------------------------

ε is a header only library. Add all files to your project and include the needed portion. Then add the include path to ``epsilon/include`` to your compile options/IDE.

The best way to include ε into a project is to use git-submodules should your parent project be a git repository too.

All interfaces are declared in ``include/ei/<xyz>.hpp`` files. You never need to
include or look into files from the ``include/ei/details`` or ``src`` directory.

The hierarchy of header-files is as follows.
```
config -- elementarytypes -- vector -|- 2dtypes -- 2dintersection
                                     |- 3dtypes -- 3dintersection
```
Each header includes all its depenencies. So, dependent on what functionality you need you need to include only one of them:

  * *config.hpp*: contains some defines to change basis properies (e.g. utf8 support for things like ε::π).
  * *elementarytypes.hpp*: uint, int8, ... (coming soon: fixed point), lerp, min/max with variable argument count
  * *vector.hpp*: Vec3, IVec2, Mat3x3, ... inclusive all expected functions (dot, sum, min/max, ..., many more)
  * *2dtypes.hpp*: shapes for 2D inclusive area, center, ... functions. Currently 2D support is relative minimalistic (feature request are welcome, but I did not use these functions until now)
  * *2dintersection.hpp*: adds the distance() and intersection() methods to 2dtypes
  * *3dtypes.hpp*: shapes for 3D inclusive surface, volume, center, ... functions. Also adds numerous "conversion" methods e.g. bounding box for a set of points.
  * *3dintersection.hpp*: adds the distance() and intersection() methods to 3dtypes  

The configuration system of epsilon works as follows:

All files expect to include the config.hpp (or depend on it). The repository only contains a defaultconfig.hpp. This should be copied and renamed accordingly. Thus, you can use epsilon as submodule in other repositories. Since the new (renamed) config file is not part of this repository you can change it as you want.
You only need to make sure the file is found from your compiler if parsing the epsilon files.


Pretty Printing in GDB
-------------------------------------------------------------------------------

In the directory ``doc`` you can find three files which enable pretty printing in GDB. To load the pretty printer ``epsilon_pretty_gdb_load.py`` must be imported from GDB. This is achieved by either placing ``.gdbinit`` into the execution directory of your application or into the user's directory and editing the path inside this file. The path must point to the ``epsilon/doc`` directory which contains the python files.

Note, if your user directory already contains a ``.gdbinit`` file, open it and edit it with the lines from ``epsilon/doc/.gdbinit``. The important point is to add the proper search path and then call ``import epsilon_pretty_gdb_load.py``

Range of Functions
-------------------------------------------------------------------------------

  * Elementary types: sized ints uint64... (well meanwhile you can get them from C++ std library, but the int16_t still has this ugly _t)
  * Template row and column vectors of different sizes and elementary types (don't be confused: the library defines vectors as matrices with one dimension set to 1)
  * Matrices
	  * transformations: rotation, translation, projection
	  * determinant and other matrix utilities
	  * LU decomposition, 2D/3D spectral decomposition
  * Lots of 2D and 3D shapes
	  * area(), volume(), surface() and center()
	  * conversion/construction of bounding shapes for other shapes
	  * distance() functions
	  * intersection functions
	  * transform methods


Available Intersection Methods
-------------------------------------------------------------------------------

The following table gives an overview over the implemented intersection methods.
Also, it shows the median performance in milliseconds per 1 million elements.
Even though, 5000 iterations, with 2500 test each, are performed the median (and every other measure I tried) has a noticeable
variance from more than 15% between application runs due to CPU boosts.
Cells with a * are implemented but not measured.

Benchmarkconfigs: VS2015, /O2, **Win32** on an i7-4950S

|                 | Box  | Cap. | Cone | Disc | DOP  | Ell. | Fru. | Line | OBox | OEl. | Pla. | Poi. | Ray  | Sph. | The. | Tri. |
|-----------------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
|**Box**          | 10.1 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Capsule**      |      | 35.0 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Cone**         |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Disc**         |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**DOP**          |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Ellipsoid**    |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Frustum**      | *    |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Line**         |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**OBox**         |      |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**OEllipsoid**   |      |      |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- |
|**Plane**        | 6.27 |      |      |      |      |      |      |      | 9.73 |      |      | ---- | ---- | ---- | ---- | ---- |
|**Point**        | 8.19 | 15.2 | 10.0 | ---- | 8.06 | 6.35 | *    | ---- | 19.1 | 12.9 | ---- | ---- | ---- | ---- | ---- | ---- |
|**Ray**          | 16.4 |      |      |      |      | 17.3 | *    | ---- | 36.0 |      |      | ---- | ---- | ---- | ---- | ---- |
|**Sphere**       | 10.7 | 15.8 |      |      |      |      |      |      |      |      | 5.71 | 3.33 | 5.71 | 6.14 | ---- | ---- |
|**Thetrahedron** |      |      |      |      |      |      |      |      |      |      |      | 24.9 |      |      |      | ---- |
|**Triangle**     | 17.3 |      | 53.5 |      |      |      |      |      | 33.4 |      |      | ---- | 22.4 | 38.2 |      |      |

Benchmarkconfigs: VS2015, /O2, **x64** on an i7-4950S

|                 | Box  | Cap. | Cone | Disc | DOP  | Ell. | Fru. | Line | OBox | OEl. | Pla. | Poi. | Ray  | Sph. | The. | Tri. |
|-----------------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
|**Box**          | 9.60 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Capsule**      |      | 44.0 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Cone**         |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Disc**         |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**DOP**          |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Ellipsoid**    |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Frustum**      | *    |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Line**         |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**OBox**         |      |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**OEllipsoid**   |      |      |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- |
|**Plane**        | 4.99 |      |      |      |      |      |      |      | 12.4 |      |      | ---- | ---- | ---- | ---- | ---- |
|**Point**        | 6.91 | 13.3 | 7.81 | ---- | 6.27 | 5.89 | *    | ---- | 26.6 | 13.1 | ---- | ---- | ---- | ---- | ---- | ---- |
|**Ray**          | 14.5 |      |      |      |      | 16.3 | *    | ---- | 45.3 |      |      | ---- | ---- | ---- | ---- | ---- |
|**Sphere**       | 8.32 | 13.8 |      |      |      |      |      |      |      |      | 1.79 | 1.92 | 4.22 | 3.71 | ---- | ---- |
|**Thetrahedron** |      |      |      |      |      |      |      |      |      |      |      | 18.9 |      |      |      | ---- |
|**Triangle**     | 15.9 |      | 42.8 |      |      |      |      |      | 41.7 |      |      | ---- | 18.6 | 35.1 |      |      |
