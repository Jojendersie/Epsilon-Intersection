Epsilon Intersection Library - ε
===============================================================================

ε is a super flexible and simple to use **C++ library** designed for games. It is computing **intersections** of various three dimensional bodies and containes different utilities you could need for culling, ray casting or similar functions.

It provides powerful template based matrix/vector and quaternion classes which come with a ton of useful functions.
Due to its lean design you can easily integrate the parts you need into your project, just by adding the source!


How to use the library?
-------------------------------------------------------------------------------

There are 2 different ways to include ε into your project:

  1. Add all files to your current project and compile as usual.
  2. Compile a library and use that

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
The numbers are performance indices (numbers relative to lensq(cross(a,b)) timing).
They can be used to compare a method's relative execution speed.
Cells with a * are implemented but not measured.

Benchmarkconfigs: VS2015, /O2, Win32 on an i7-4950S

|                 | Box  | Cap. | Disc | DOP  | Ell. | Fru. | Line | OBox | OEl. | Pla. | Poi. | Ray  | Sph. | The. | Tri. |
|-----------------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
|**Box**          | 3.77 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Capsule**      |      | 13.5 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Disc**         |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**DOP**          |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Ellipsoid**    |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Frustum**      | *    |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**Line**         |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**OBox**         |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|**OEllipsoid**   |      |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- | ---- |
|**Plane**        |      |      |      |      |      |      |      |      |      |      | ---- | ---- | ---- | ---- | ---- |
|**Point**        | 4.41 | 6.60 | ---- | 3.72 | 2.66 | *    | ---- | 8.42 | 5.41 | ---- | ---- | ---- | ---- | ---- | ---- |
|**Ray**          | 8.25 |      |      |      | 8.43 | *    | ---- | 15.2 |      |      | ---- | ---- | ---- | ---- | ---- |
|**Sphere**       | 5.04 | 6.86 |      |      |      |      |      |      |      | 1.28 | 1.32 | 2.48 | 1.42 | ---- | ---- |
|**Thetrahedron** |      |      |      |      |      |      |      |      |      |      | 9.22 |      |      |      | ---- |
|**Triangle**     | 9.53 |      |      |      |      |      |      | 14.8 |      |      | ---- | 9.66 | 16.9 |      |      |
