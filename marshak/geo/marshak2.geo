// Gmsh project created on Mon Nov 14 01:30:32 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.5, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0.5, 0, 0, 1.0};
//+
Point(4) = {0.51, 0, 0, 1.0};
//+
Point(5) = {0.51, 0.01, 0, 1.0};
//+
Point(6) = {0.5, 0.01, 0, 1.0};
//+
Point(7) = {0, 0.01, 0, 1.0};
//+
Point(8) = {-0.5, 0.01, 0, 1.0};
//+
Line(1) = {8, 7};
//+
Line(2) = {7, 6};
//+
Line(3) = {6, 5};
//+
Line(4) = {5, 4};
//+
Line(5) = {4, 3};
//+
Line(6) = {3, 2};
//+
Line(7) = {2, 1};
//+
Line(8) = {1, 8};
//+
Curve Loop(1) = {8, 1, 2, 3, 4, 5, 6, 7};
//+
Line(9) = {2, 7};
//+
Line(10) = {6, 3};
//+
Curve Loop(2) = {1, -9, 7, 8};
//+
Plane Surface(1) = {2};
//+
Curve Loop(3) = {6, 9, 2, 10};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {3, 4, 5, -10};
//+
Plane Surface(3) = {4};
//+
Physical Curve(11) = {7, 6, 5, 4, 3, 2, 1, 8};
//+
Physical Curve(11) -= {7, 2, 6, 4, 1, 8, 5, 3};
//+
Physical Curve("boundary", 21) = {8, 1, 7, 6, 5, 4, 3, 2};
//+
Physical Surface("inner1", 31) = {2};
//+
Physical Surface("ghost1", 32) = {1};
//+
Physical Surface("ghost2", 33) = {3};
