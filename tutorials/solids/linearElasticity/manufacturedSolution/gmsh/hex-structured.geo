// Gmsh .geo file to create a mesh of a cube LxLxL

SetFactory("OpenCASCADE");

// Mesh spacing parameters
Include "meshSpacing.geo";

// Cube side length
L = 0.2;

N = L/dx + 1;

// Define the cube using the built-in Box command.
Box(1) = {0, 0, 0, L, L, L};

// Apply the Transfinite algorithm to the cube's edges.
// This forces a structured grid with N points along each edge.
Transfinite Line { Abs(Boundary{ Surface{ Abs(Boundary{ Volume{1}; }) }; }) } = N;

// Apply the Transfinite algorithm to the cube's surfaces.
Transfinite Surface { Abs(Boundary{ Volume{1}; }) };

// Apply the Transfinite algorithm to the volume. This tells Gmsh to
// create a structured 3D mesh by extruding the surface meshes.
Transfinite Volume {1};

// Recombine the triangular surface elements into quadrilaterals.
// This is a necessary step before Gmsh can create hexahedra in the volume.
Recombine Surface { Abs(Boundary{ Volume{1}; }) };

// Generate the 3D mesh
Mesh 3;
