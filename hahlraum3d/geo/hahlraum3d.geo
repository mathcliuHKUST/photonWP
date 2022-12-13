// Gmsh project created on Fri Nov 25 23:11:47 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 0.52, 0.5, 0.5*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 0.52, 0.52, 0.5*Pi};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
//+
Cylinder(3) = {0, 0, 0.05, 0, 0, 0.1, 0.52, 0.5*Pi};
//+
Cylinder(4) = {0, 0, 0.25, 0, 0, 0.1, 0.52, 0.5*Pi};
//+
BooleanFragments{ Volume{2}; Delete; }{ Volume{3}; Volume{4}; Delete; }
//+
Recursive Delete {
  Volume{6}; 
}
//+
Recursive Delete {
  Volume{7}; 
}
//+
Cylinder(6) = {0, 0, 0, 0, 0, 0.52, 0.52, 0.5*Pi};
//+
BooleanFragments{ Volume{6}; Delete; }{ Volume{5}; Volume{4}; Volume{3}; Volume{2}; Volume{1}; Delete; }
//+
Cylinder(7) = {0, 0, 0.5, 0, 0, 0.02, 0.52, 0.5*Pi};
//+
BooleanFragments{ Volume{1}; Volume{6}; Delete; }{ Volume{7}; Delete; }
//+
Sphere(10) = {0, 0, 0, 0.25, 0, Pi/2, 0.5*Pi};
//+
BooleanFragments{ Volume{9}; Delete; }{ Volume{10}; Delete; }
//+
Cylinder(12) = {0, 0, 0, 0, 0, 0.52, 0.52, 2*Pi/32};
//+
BooleanFragments{ Volume{12}; Delete; }{ Volume{10}; Volume{5}; Volume{4}; Volume{3}; Volume{11}; Volume{2}; Volume{7}; Volume{6}; Volume{8}; Delete; }
//+
Recursive Delete {
  Volume{10}; Volume{11}; Volume{12}; Volume{13}; Volume{15}; Volume{16}; Volume{17}; Volume{18}; Volume{14}; 
}
//+
Physical Surface("boundary", 21) = {36, 32, 29, 17, 8, 3, 14, 22, 40, 2, 12, 5, 13, 10, 19, 31, 39, 35, 21, 42, 1, 7, 16, 28, 33, 37, 38, 26, 41};
//+
Physical Volume("ghost1", 301) = {1};
//+
Physical Volume("ghost2", 302) = {8, 7, 4, 6, 2};
//+
Physical Volume("inner1", 311) = {5};
//+
Physical Volume("inner2", 312) = {9, 3};
