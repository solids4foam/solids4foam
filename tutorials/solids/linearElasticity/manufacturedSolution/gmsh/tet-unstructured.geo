// Gmsh .geo file to create a mesh of a cube LxLxL

// Mesh spacing parameters
Include "meshSpacing.geo";

// Cube side length
L = 0.2;

// Define the points
Point(1) = {0, 0, 0, dx};  // Point at the origin
Point(2) = {L, 0, 0, dx};  // Along the x-axis
Point(3) = {L, L, 0, dx};  // Along the xy-plane
Point(4) = {0, L, 0, dx};  // Along the y-axis
Point(5) = {0, 0, L, dx};  // Along the z-axis
Point(6) = {L, 0, L, dx};  // Along the xz-plane
Point(7) = {L, L, L, dx};  // Opposite corner of the cube
Point(8) = {0, L, L, dx};  // Along the yz-plane

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

// Generate the mesh
Mesh 3;
