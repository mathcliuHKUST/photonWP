// Gmsh project created on Tue Nov 08 15:52:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 7, 0, 1.0};
//+
Recursive Delete {
  Point{2}; 
}
//+
Point(2) = {0, 2, 0, 1.0};
//+
Point(3) = {7, 2, 0, 1.0};
//+
Point(4) = {7, 0, 0, 1.0};
//+
Point(5) = {0, 0.5, 0, 1.0};
//+
Point(6) = {2, 0.5, 0, 1.0};
//+
Recursive Delete {
  Point{6}; 
}
//+
Point(6) = {2.5, 0.5, 0, 1.0};
//+
Point(7) = {2.5, 1.5, 0, 1.0};
//+
Point(8) = {3, 0, 0, 1.0};
//+
Point(9) = {3, 1, 0, 1.0};
//+
Point(10) = {4, 1, 0, 1.0};
//+
Point(11) = {4, 0, 0, 1.0};
//+
Point(12) = {4.5, 0.5, 0, 1.0};
//+
Point(13) = {7, 0.5, 0, 1.0};
//+
Point(14) = {4.5, 1.5, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 13};
//+
Line(3) = {13, 12};
//+
Line(4) = {12, 14};
//+
Line(5) = {14, 7};
//+
Line(6) = {7, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {5, 2};
//+
Line(9) = {8, 9};
//+
Line(10) = {9, 10};
//+
Line(11) = {10, 11};
//+
Line(12) = {11, 8};
//+
Line(13) = {5, 1};
//+
Line(14) = {1, 8};
//+
Line(15) = {11, 4};
//+
Line(16) = {4, 13};
//+
Curve Loop(1) = {8, 1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 10, 11, 12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 13, 14, 9, 10, 11, 15, 16, 3, 4, 5, 6};
//+
Surface(3) = {3};
//+
Point(15) = {-2, 0, 0, 1.0};
//+
Point(16) = {-2, 2, 0, 1.0};
//+
Point(17) = {9, 2, 0, 1.0};
//+
Point(18) = {9, 0, 0, 1.0};
//+
Line(17) = {16, 2};
//+
Line(18) = {16, 16};
//+
Line(18) = {16, 15};
//+
Line(19) = {15, 1};
//+
Line(20) = {3, 17};
//+
Line(21) = {17, 18};
//+
Line(22) = {18, 4};
//+
Curve Loop(5) = {18, 19, -13, 8, -17};
//+
Plane Surface(4) = {5};
//+
Recursive Delete {
  Curve{20}; 
}
//+
Recursive Delete {
  Curve{21}; 
}
//+
Recursive Delete {
  Curve{22}; 
}
//+
Point(17) = {7.5, 0, 0, 1.0};
//+
Point(18) = {7.5, 2, 0, 1.0};
//+
Line(20) = {18, 18};
//+
Line(20) = {3, 18};
//+
Line(21) = {18, 17};
//+
Line(22) = {17, 4};
//+
Curve Loop(6) = {2, -16, -22, -21, -20};
//+
Surface(5) = {6};
//+
Physical Curve("boundary", 21) = {18, 17, 1, 20, 21, 22, 15, 12, 14, 19};
//+
Physical Surface("inner1", 311) = {3};
//+
Physical Surface("inner2", 312) = {1, 2};
//+
Physical Surface("ghost1", 321) = {4};
//+
Physical Surface("ghost2", 322) = {5};
