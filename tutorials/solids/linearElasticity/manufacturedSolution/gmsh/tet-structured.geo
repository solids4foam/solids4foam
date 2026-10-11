// Gmsh .geo file to create a mesh of a cube LxLxL

// Mesh spacing parameters
Include "meshSpacing.geo";

// Cube side length
L = 0.2;

// Number of cells per edge for structured
nCellPerEdge = L/dx;

// Define the points
Point(1) = {0, 0, 0, 1};  // Point at the origin
Point(2) = {L, 0, 0, 1};  // Along the x-axis
Point(3) = {L, L, 0, 1};  // Along the xy-plane
Point(4) = {0, L, 0, 1};  // Along the y-axis
Point(5) = {0, 0, L, 1};  // Along the z-axis
Point(6) = {L, 0, L, 1};  // Along the xz-plane
Point(7) = {L, L, L, 1};  // Opposite corner of the cube
Point(8) = {0, L, L, 1};  // Along the yz-plane

// Define the lines (edges) connecting the points
Line(1) = {1, 2};  // Edge along x-axis
Line(2) = {2, 3};  // Edge along xy-plane
Line(3) = {3, 4};  // Edge along y-axis
Line(4) = {4, 1};  // Edge closing the bottom face

Line(5) = {1, 5};  // Edge along z-axis
Line(6) = {2, 6};  // Edge along z-axis
Line(7) = {3, 7};  // Edge along z-axis
Line(8) = {4, 8};  // Edge along z-axis

Line(9) = {5, 6};  // Edge along xz-plane
Line(10) = {6, 7}; // Edge along xz-plane
Line(11) = {7, 8}; // Edge along yz-plane
Line(12) = {8, 5}; // Edge closing the top face

// Define the surfaces (faces) of the cube
Line Loop(1) = {1, 2, 3, 4};     // Bottom face
Plane Surface(1) = {1};

Line Loop(2) = {5, 9, -6, -1};   // Front face
Plane Surface(2) = {2};

Line Loop(3) = {6, 10, -7, -2};  // Right face
Plane Surface(3) = {3};

Line Loop(4) = {7, 11, -8, -3};  // Back face
Plane Surface(4) = {4};

Line Loop(5) = {8, 12, -5, -4};  // Left face
Plane Surface(5) = {5};

Line Loop(6) = {9, 10, 11, 12};  // Top face
Plane Surface(6) = {6};

// Define the volume (body) of the cube
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Apply transfinite meshing to edges
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = nCellPerEdge Using Progression 1;

// Apply transfinite meshing to surfaces
Transfinite Surface {1, 2, 3, 4, 5, 6};
//Recombine Surface {1, 2, 3, 4, 5, 6};  // Ensure quadrilateral surface meshing

// Apply transfinite meshing to the volume
Transfinite Volume {1};

// Define a single physical boundary patch that includes all 6 surfaces
//Physical Surface("AllBoundaries") = {1, 2, 3, 4, 5, 6};

// Generate the structured mesh
Mesh 3;
