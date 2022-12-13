// Gmsh project created on Sat Nov 19 00:40:52 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, 1, 0.5*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 1, 1.03, 0.5*Pi};
//+
Surface Loop(3) = {7, 6, 9, 8, 10};
//+
Surface Loop(4) = {2, 1, 4, 3, 5};
//+
Volume(3) = {3, 4};
//+
Surface Loop(5) = {7, 6, 9, 8, 10};
//+
Surface Loop(6) = {10, 6, 7, 9, 8};
//+
Surface Loop(7) = {10, 6, 7, 9, 8};
//+
Surface Loop(8) = {8, 9, 6, 7, 10};
//+
Surface Loop(9) = {3, 4, 1, 2, 5};
//+
Surface Loop(10) = {5, 1, 2, 4, 3};
//+
Surface Loop(11) = {5, 1, 2, 4, 3};
//+
Surface Loop(12) = {5, 1, 2, 4, 3};
//+
Surface Loop(13) = {5, 1, 2, 4, 3};
