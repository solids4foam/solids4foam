// Gmsh .geo file to create a mesh of a cube with a spherical cavity

SetFactory("OpenCASCADE");  // Use OpenCASCADE kernel

// Define parameters
L = 0.5;          // Length of the cube
R = 0.2;          // Radius of the spherical cavity

// Mesh spacing parameters
minDeltaX = 0.01;   // Minimum mesh size near the cavity
maxDeltaX = 0.05;   // Maximum mesh size away from the cavity
minDist = 0.0;       // Distance at which minDeltaX is applied
maxDist = 0.3;       // Distance at which maxDeltaX is applied

// Define the cube
Box(1) = {0, 0, 0, L, L, L};

// Define the spherical cavity
Sphere(2) = {0, 0, 0, R};

// Perform a boolean subtraction to create the cavity in the cube
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Patches
Physical Surface("posX", 19) = {7};
Physical Surface("posY", 20) = {4};
Physical Surface("posZ", 21) = {3};
Physical Surface("negX", 17) = {1};
Physical Surface("negY", 18) = {2};
Physical Surface("negZ", 22) = {5};
Physical Surface("cavity", 23) = {6};

// Define physical volume for the 3D region
Physical Volume("volume") = {1};

// Define a distance field around the cavity
Field[1] = Distance;
Field[1].FacesList = {6};  // Surface ID of the cavity

// Define a threshold field for mesh size transition
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = minDeltaX;   // Minimum mesh size near the cavity
Field[2].LcMax = maxDeltaX;   // Maximum mesh size away from the cavity
Field[2].DistMin = minDist; // Distance at which LcMin is applied
Field[2].DistMax = maxDist;  // Distance at which LcMax is applied

// Use the threshold field to control mesh size
Background Field = 2;

// Mesh the geometry
Mesh 3;
