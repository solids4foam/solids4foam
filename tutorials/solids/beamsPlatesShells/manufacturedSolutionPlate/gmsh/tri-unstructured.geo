// Gmsh .geo file: one layer of triangular prisms over the L x W plate, for
// the finite area mesh on the front surface

// Mesh spacing
Include "meshSpacing.geo";

// Plate side lengths and mesh thickness
L = 10;
W = 10;
h = 0.1;

Point(1) = {0, 0, 0, dx};
Point(2) = {L, 0, 0, dx};
Point(3) = {L, W, 0, dx};
Point(4) = {0, W, 0, dx};

Line(1) = {1, 2};  // bottom, y = 0
Line(2) = {2, 3};  // right, x = L
Line(3) = {3, 4};  // top, y = W
Line(4) = {4, 1};  // left, x = 0

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Frontal-Delaunay triangles of near-uniform size
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = dx;
Mesh.CharacteristicLengthMax = dx;

// One layer of prisms through the thickness: Recombine keeps the extruded
// elements as prisms with quadrilateral side faces. The extrusion returns
// the front surface, the volume and the four side surfaces in the order of
// the lines of the loop
out[] = Extrude {0, 0, h} { Surface{1}; Layers{1}; Recombine; };

Physical Surface("back") = {1};
Physical Surface("front") = {out[0]};
Physical Surface("bottom") = {out[2]};
Physical Surface("right") = {out[3]};
Physical Surface("top") = {out[4]};
Physical Surface("left") = {out[5]};
Physical Volume("plate") = {out[1]};

Mesh 3;
