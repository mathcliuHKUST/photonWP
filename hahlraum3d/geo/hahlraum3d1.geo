// Gmsh project created on Tue Dec 06 23:34:30 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, -1, 0, 0, 2, 0.5, 2*Pi};
//+
Cylinder(2) = {0, 0, -1, 0, 0, 0.2, 0.5, 2*Pi};
//+
Cylinder(3) = {0, 0, 0.8, 0, 0, 0.2, 0.5, 2*Pi};
//+
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
BooleanFragments{ Volume{4}; Delete; }{ Volume{3}; Delete; }
//+
Recursive Delete {
  Volume{2}; 
}
