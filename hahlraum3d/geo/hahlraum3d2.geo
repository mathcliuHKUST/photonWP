// Gmsh project created on Sat Nov 19 00:07:11 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, 0.5, 0.5*Pi};
//+
Sphere(2) = {0, 0, 0, 0.2, 0, Pi/2, 0.5*Pi};
//+
Cylinder(3) = {0, 0, 0, 0, 0, 1, 0.52, 0.5*Pi};
//+
Cylinder(4) = {0, 0, 0.5, 1, 1, 0, 0.1, 2*Pi};
//+
Surface Loop(5) = {4, 3, 1, 2, 5};
//+
Surface Loop(6) = {1, 2, 5, 4, 3};
//+
Surface Loop(7) = {15, 16, 17};
//+
Surface Loop(8) = {11, 10, 13, 12, 14};
//+
Volume(5) = {5, 6, 7, 8};
//+
Surface Loop(9) = {12, 13, 10, 11, 14};
//+
Surface Loop(10) = {13, 12, 10, 11, 14};
