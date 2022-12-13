// Gmsh project created on Wed Dec 07 09:49:57 2022
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, -1, 0, 0, 2, 0.5, 2*Pi};
//+
Cylinder(2) = {0, 0, -1, 0, 0, 0.04, 0.5, 2*Pi};
//+
Cylinder(3) = {0, 0, 0.96, 0, 0, 0.04, 0.5, 2*Pi};


//+
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
BooleanFragments{ Volume{4}; Delete; }{ Volume{3}; Delete; }
//+
Recursive Delete {
  Volume{6}; 
}
//+
Recursive Delete {
  Volume{9}; 
}

//+
Cylinder(5) = {0, 0, -1, 0, 0, 2, 0.48, 2*Pi};

//+
BooleanFragments{ Volume{2}; Volume{4}; Volume{3}; Delete; }{ Volume{5}; Delete; }
//+
Recursive Delete {
  Volume{7}; 
}
//+
Cylinder(7) = {0, 0, -1, 0, 0, 2, 0.46, 2*Pi};
//+
BooleanFragments{ Volume{4}; Volume{2}; Volume{6}; Delete; }{ Volume{7}; Delete; }
//+
Recursive Delete {
  Volume{12}; 
}
//+
Cylinder(12) = {0, 0, -1, 0, 0, 2, 0.2, 2*Pi};
//+
BooleanFragments{ Volume{9}; Delete; }{ Volume{12}; Delete; }
//+
BooleanFragments{ Volume{11}; Delete; }{ Volume{14}; Delete; }
//+
Recursive Delete {
  Volume{16}; 
}
//+
Coherence;
//+
Sphere(16) = {0, 0, 0, 0.2, -Pi/2, Pi/2, 2*Pi};

//+
BooleanFragments{ Volume{7}; Delete; }{ Volume{16}; Delete; }
//+
Recursive Delete {
  Volume{18}; 
}
