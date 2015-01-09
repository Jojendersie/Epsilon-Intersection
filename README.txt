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
config -- elementarytypes -- matrix -|- 2dtypes -- 2dfunctions -- 2dintersection
                                     |- 3dtypes -- 3dfunctions -- 3dintersection
Hence each one includes all its parents. You always need to include only the
latest one.


3. Available Methods
-------------------------------------------------------------------------------

The following table gives an overview over the implemented intersection methods.
The numbers are the performance indices (relative numbers) and can be used to
compare the method execution speed.
Benchmarkconfig: i7-4950S, VS2013, /O2, Win32

2D

3D           | Box  | Cap. | Disc | DOP  | Ell. | Fru. | Line | OBox | OEl. | Pla. | Poi. | Ray  | Sph. | The. | Tri. |
-------------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
Box          | 16.7 |      |      |      |      |      |      |      |      |      |      | 27.6 | 12.2 |      |      |
Capsule      | ---- |      |      |      |      |      |      |      |      |      |      |      |      |      |      |
Disc         | ---- | ---- |      |      |      |      |      |      |      |      |      |      |      |      |      |
DOP          | ---- | ---- | ---- |      |      |      |      |      |      |      |      |      |      |      |      |
Ellipsoid    | ---- | ---- | ---- | ---- |      |      |      |      |      |      | 6.90 | 22.7 |      |      |      |
Frustum      | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |      |      |      |
Line         | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |      |      |
OBox         | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |      |
OEllipsoid   | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |
Plane        | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |
Point        | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      | 3.33 |      |      |
Ray          | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      | 11.9 |
Sphere       | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 3.95 |      |      |
Thetrahedron | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |
Triangle     | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |
-------------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
