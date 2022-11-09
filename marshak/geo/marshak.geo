// Gmsh project created on Wed Nov 09 23:02:03 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.01, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0.2, 0, 0, 1.0};
//+
Point(4) = {0.21, 0, 0, 1.0};
//+
Point(5) = {0.002, 0, 0, 1.0};
//+
Point(6) = {0.202, 0, 0, 1.0};
//+
Point(7) = {0.202, 0.002, 0, 1.0};
//+
Point(8) = {0.2, 0.002, 0, 1.0};
//+
Point(9) = {0, 0.002, 0, 1.0};
//+
Point(10) = {-0.002, 0.002, 0, 1.0};
//+
Point(11) = {-0.002, 0, 0, 1.0};
//+
Recursive Delete {
  Point{1}; 
}
//+
Recursive Delete {
  Point{5}; 
}
//+
Recursive Delete {
  Point{4}; 
}
//+
Line(1) = {10, 11};
//+
Line(2) = {11, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 6};
//+
Line(5) = {7, 6};
//+
Recursive Delete {
  Curve{5}; 
}
//+
Recursive Delete {
  Curve{4}; 
}
//+
Point(12) = {0.202, 0, 0, 1.0};
//+
Point(13) = {0.202, 0.002, 0, 1.0};
//+
Line(4) = {10, 9};
//+
Line(5) = {10, 9};
//+
Line(6) = {9, 8};
//+
Line(7) = {8, 13};
//+
Line(8) = {13, 12};
//+
Line(9) = {12, 3};
//+
Line(10) = {9, 2};
//+
Line(11) = {8, 3};
//+
Curve Loop(1) = {1, 2, -10, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -11, -6, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 8, 9, -11};
//+
Plane Surface(3) = {3};
//+
Physical Curve("boundary", 21) = {1, 2, 3, 9, 8, 7, 6, 4};
//+
Physical Surface("inner1", 31) = {2};
//+
Physical Surface("ghost1", 32) = {1};
//+
Physical Surface("ghost2", 33) = {3};
