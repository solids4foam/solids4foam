// Gmsh .geo file to create a mesh of a cube with a spherical cavity

SetFactory("OpenCASCADE");  // Use OpenCASCADE kernel

// Define parameters
L = 0.5;          // Length of the cube
R = 0.2;          // Radius of the spherical cavity

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

// Mesh the geometry
Mesh 3;
