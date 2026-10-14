// Average mesh spacing, replaced by the verification driver
dx = 0.1;

// Plate size, hole radius, and out-of-plane depth
L = 2;
r = 0.5;
d = 0.1;

Point(1) = {r, 0, 0, dx};
Point(2) = {L, 0, 0, dx};
Point(3) = {L, L, 0, dx};
Point(4) = {0, L, 0, dx};
Point(5) = {0, r, 0, dx};
Point(6) = {0, 0, 0, dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Circle(5) = {5, 6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

newEntities[] =
    Extrude {0, 0, d}
    {
        Surface{1};
        Layers{1};
        Recombine;
    };

Physical Volume("internal") = {newEntities[1]};
Physical Surface("front") = {newEntities[0]};
Physical Surface("back") = {1};
Physical Surface("down") = {newEntities[2]};
Physical Surface("right") = {newEntities[3]};
Physical Surface("up") = {newEntities[4]};
Physical Surface("left") = {newEntities[5]};
Physical Surface("hole") = {newEntities[6]};
