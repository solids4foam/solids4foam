// Unstructured triangular-prism mesh of Cook's membrane, created with
// ./Allrun triangular. The verification driver replaces the spacing dx.

// Average mesh spacing
dx = 0.003125;

// Vertices of the trapezoid and out-of-plane depth
L = 0.048;
W = 0.044;
T = 0.060;
d = 1;

Point(1) = {0, 0, 0, dx};
Point(2) = {L, W, 0, dx};
Point(3) = {L, T, 0, dx};
Point(4) = {0, W, 0, dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Unstructured Delaunay triangulation
Mesh.Algorithm = 5;

newEntities[] =
    Extrude {0, 0, d}
    {
        Surface{1};
        Layers{1};
        Recombine;
    };

Physical Volume("internal") = {newEntities[1]};
Physical Surface("frontAndBack") = {1, newEntities[0]};
Physical Surface("topAndBottom") = {newEntities[2], newEntities[4]};
Physical Surface("right") = {newEntities[3]};
Physical Surface("left") = {newEntities[5]};
