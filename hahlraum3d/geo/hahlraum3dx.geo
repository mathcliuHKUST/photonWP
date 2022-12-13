// Gmsh project created on Wed Dec 07 11:48:39 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, 0.5, 0.5*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 0.96, 0.5, 0.5*Pi};
//+
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Cylinder(4) = {0, 0, 0, 0, 0, 1, 0.48, 0.5*Pi};
//+
BooleanFragments{ Volume{3}; Volume{2}; Delete; }{ Volume{4}; Delete; }
//+
Cylinder(5) = {0, 0, 0, 0, 0, 1, 0.46, 0.5*Pi};
//+
BooleanFragments{ Volume{2}; Volume{1}; Volume{4}; Volume{3}; Delete; }{ Volume{5}; Delete; }
//+
Cylinder(8) = {0, 0, 0.96, 0, 0, 0.04, 0.2, 0.5*Pi};
//+
BooleanFragments{ Volume{4}; Delete; }{ Volume{8}; Delete; }
//+
Sphere(10) = {0, 0, -0, 0.2, -0, Pi/2, 0.5*Pi};
//+
Coherence;
//+
Cylinder(12) = {0, 0, 0, 0, 0, 1, 0.5, Pi/16};
//+
BooleanFragments{ Volume{10}; Volume{8}; Volume{11}; Volume{3}; Volume{7}; Volume{9}; Volume{5}; Volume{1}; Delete; }{ Volume{12}; Delete; }
//+
Recursive Delete {
  Volume{2}; Volume{4}; Volume{5}; Volume{7}; Volume{9}; Volume{11}; Volume{13}; Volume{15}; 
}
//+
Physical Volume("ghost0", 90) = {10};
//+
Physical Volume(" ghost0", 90) -= {10};
//+
Physical Volume("ghost0", 301) = {10};
//+
Physical Volume("inner1", 311) = {6};
//+
Physical Volume("inner2", 312) = {1};
//+
Physical Volume("inner3", 313) = {8, 16, 14, 12};
//+
Physical Volume("inner1", 311) += {3};
//+
Physical Surface("boundary", 20) = {21, 2, 47, 10, 54, 59, 28, 40, 32, 62, 12, 49, 63, 56, 42, 34, 22, 4, 9, 48, 64, 55, 3, 25, 33, 43};
