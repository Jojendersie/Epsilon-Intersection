﻿Epsilon Intersection Library - ε
===============================================================================

ε is a super flexible and simple to use **C++ library** for computing **intersections** of various three dimensional bodies. 

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
config -- elementarytypes -- vector -|- 2dtypes -- 2dfunctions -- 2dintersection
                                     |- 3dtypes -- 3dfunctions -- 3dintersection
```
Hence, each file includes all its parents. You always need to include only the
top one.



Available Methods
-------------------------------------------------------------------------------

The following table gives an overview over the implemented intersection methods.
The numbers are performance indices (relative numbers) and can be used to
compare a method's realtive execution speed.
Benchmarkconfig: i7-4950S, VS2013, /O2, Win32

                 | Box  | Cap. | Disc | DOP  | Ell. | Fru. | Line | OBox | OEl. | Pla. | Poi. | Ray  | Sph. | The. | Tri. |
-----------------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
**Box**          | 16.7 |      |      |      |      |      |      |      |      |      |      | 27.6 | 12.2 |      |      |
**Capsule**      | ---- |      |      |      |      |      |      |      |      |      |      |      |      |      |      |
**Disc**         | ---- | ---- |      |      |      |      |      |      |      |      |      |      |      |      |      |
**DOP**          | ---- | ---- | ---- |      |      |      |      |      |      |      |      |      |      |      |      |
**Ellipsoid**    | ---- | ---- | ---- | ---- |      |      |      |      |      |      | 6.90 | 22.7 |      |      |      |
**Frustum**      | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |      |      |      |
**Line**         | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |      |      |
**OBox**         | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |      |
**OEllipsoid**   | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |      |
**Plane**        | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      |      |      |      |
**Point**        | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      | 3.33 |      |      |
**Ray**          | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |      | 11.9 |
**Sphere**       | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 3.95 |      |      |
**Thetrahedron** | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |      |
**Triangle**     | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |      |