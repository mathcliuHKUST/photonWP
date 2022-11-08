// Gmsh project created on Wed Nov 09 00:52:12 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1.5, 0, 0, 1.0};
//+
Point(4) = {5.5, 0, 0, 1.0};
//+
Point(5) = {9.5, 0, 0, 1.0};
//+
Point(6) = {13.5, 0, 0, 1.0};
//+
Point(7) = {14, 0, 0, 1.0};
//+
Point(8) = {15, 0, 0, 1.0};
//+
Point(9) = {-4, 0, 0, 1.0};
//+
Point(10) = {1, 4.5, 0, 1.0};
//+
Point(11) = {1.5, 4.5, 0, 1.0};
//+
Point(12) = {5.5, 4.5, 0, 1.0};
//+
Point(13) = {9.5, 4.5, 0, 1.0};
//+
Point(14) = {13.5, 4.5, 0, 1.0};
//+
Point(15) = {14, 4.5, 0, 1.0};
//+
Point(16) = {15, 4.5, 0, 1.0};
//+
Point(17) = {-4, 6, 0, 1.0};
//+
Point(18) = {1, 6, 0, 1.0};
//+
Point(19) = {13.5, 6, 0, 1.0};
//+
Point(20) = {14, 6, 0, 1.0};
//+
Point(21) = {1, 6.5, 0, 1.0};
//+
Point(22) = {14, 6.5, 0, 1.0};
//+
Point(23) = {-4, 6.5, 0, 1.0};
//+
Point(24) = {15, 6.5, 0, 1.0};
//+
Recursive Delete {
  Point{17}; 
}
//+
Recursive Delete {
  Point{14}; 
}
//+
Recursive Delete {
  Point{15}; 
}
//+
Recursive Delete {
  Point{20}; 
}
//+
Recursive Delete {
  Point{16}; 
}
//+
Point(25) = {0, 6.5, 0, 1.0};
//+
Line(1) = {23, 25};
//+
Line(2) = {25, 21};
//+
Line(3) = {21, 22};
//+
Line(4) = {22, 24};
//+
Line(5) = {24, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 6};
//+
Line(8) = {6, 5};
//+
Line(9) = {5, 4};
//+
Line(10) = {4, 3};
//+
Line(11) = {3, 2};
//+
Line(12) = {2, 1};
//+
Line(13) = {1, 9};
//+
Line(14) = {9, 23};
//+
Line(15) = {1, 25};
//+
Line(16) = {2, 10};
//+
Line(17) = {10, 11};
//+
Line(18) = {11, 3};
//+
Line(19) = {18, 21};
//+
Line(20) = {21, 22};
//+
Line(21) = {22, 7};
//+
Line(22) = {18, 19};
//+
Line(23) = {19, 6};
//+
Line(24) = {4, 12};
//+
Line(25) = {12, 13};
//+
Line(26) = {13, 5};
//+
Curve Loop(1) = {14, 1, -15, 13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {15, 2, -19, 22, 23, 8, -26, -25, -24, 10, -18, -17, -16, 12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {21, -6, -5, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {22, 23, -7, -21, -3, -19};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {16, 17, 18, 11};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {25, 26, 9, 24};
//+
Plane Surface(6) = {6};
//+
Physical Curve("boundary", 21) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
//+
Physical Surface("inner1", 311) = {2};
//+
Physical Surface("inner2", 312) = {4, 5, 6};
//+
Physical Surface("ghost1", 321) = {1};
//+
Physical Surface("ghost2", 322) = {3};
