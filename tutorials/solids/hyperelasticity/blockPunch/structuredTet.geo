// Tetrahedral mesh for the 3D hyperelastic punch benchmark in Xu et al.
//
// Geometry: L1 = L2 = 2, H = 1.
// Loading: pressure on one quarter of the top surface.
//
// The default mesh follows the paper's refinement convention:
// N divisions through the height and 2*N divisions along each in-plane side.

SetFactory("Built-in");

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm3D = 1;

// Geometry
L1 = 2.0;
L2 = 2.0;
H = 1.0;

// Refinement parameter from the paper: regular mesh N x 2N x 2N
N = 8;

x0 = 0.0;
x1 = 0.5*L2;
x2 = L2;

y0 = 0.0;
y1 = 0.5*L1;
y2 = L1;

z0 = 0.0;
z1 = H;

lc = H/N;

// Mesh points per half in-plane segment and through the height
nHalf = N + 1;
nHeight = N + 1;

Point(1) = {x0, y0, z0, lc};
Point(2) = {x1, y0, z0, lc};
Point(3) = {x2, y0, z0, lc};
Point(4) = {x0, y1, z0, lc};
Point(5) = {x1, y1, z0, lc};
Point(6) = {x2, y1, z0, lc};
Point(7) = {x0, y2, z0, lc};
Point(8) = {x1, y2, z0, lc};
Point(9) = {x2, y2, z0, lc};

Point(11) = {x0, y0, z1, lc};
Point(12) = {x1, y0, z1, lc};
Point(13) = {x2, y0, z1, lc};
Point(14) = {x0, y1, z1, lc};
Point(15) = {x1, y1, z1, lc};
Point(16) = {x2, y1, z1, lc};
Point(17) = {x0, y2, z1, lc};
Point(18) = {x1, y2, z1, lc};
Point(19) = {x2, y2, z1, lc};

// Bottom x-lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {7, 8};
Line(6) = {8, 9};

// Bottom y-lines
Line(7) = {1, 4};
Line(8) = {4, 7};
Line(9) = {2, 5};
Line(10) = {5, 8};
Line(11) = {3, 6};
Line(12) = {6, 9};

// Top x-lines
Line(13) = {11, 12};
Line(14) = {12, 13};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {17, 18};
Line(18) = {18, 19};

// Top y-lines
Line(19) = {11, 14};
Line(20) = {14, 17};
Line(21) = {12, 15};
Line(22) = {15, 18};
Line(23) = {13, 16};
Line(24) = {16, 19};

// Vertical lines
Line(25) = {1, 11};
Line(26) = {2, 12};
Line(27) = {3, 13};
Line(28) = {4, 14};
Line(29) = {5, 15};
Line(30) = {6, 16};
Line(31) = {7, 17};
Line(32) = {8, 18};
Line(33) = {9, 19};

Transfinite Line {1:24} = nHalf;
Transfinite Line {25:33} = nHeight;

// Bottom surfaces
Line Loop(101) = {1, 9, -3, -7};
Plane Surface(101) = {101};
Line Loop(102) = {2, 11, -4, -9};
Plane Surface(102) = {102};
Line Loop(103) = {3, 10, -5, -8};
Plane Surface(103) = {103};
Line Loop(104) = {4, 12, -6, -10};
Plane Surface(104) = {104};

// Top surfaces; surface 201 is the loaded quarter
Line Loop(201) = {13, 21, -15, -19};
Plane Surface(201) = {201};
Line Loop(202) = {14, 23, -16, -21};
Plane Surface(202) = {202};
Line Loop(203) = {15, 22, -17, -20};
Plane Surface(203) = {203};
Line Loop(204) = {16, 24, -18, -22};
Plane Surface(204) = {204};

// Boundary side surfaces
Line Loop(301) = {1, 26, -13, -25};
Plane Surface(301) = {301};
Line Loop(302) = {2, 27, -14, -26};
Plane Surface(302) = {302};
Line Loop(303) = {5, 32, -17, -31};
Plane Surface(303) = {303};
Line Loop(304) = {6, 33, -18, -32};
Plane Surface(304) = {304};
Line Loop(305) = {7, 28, -19, -25};
Plane Surface(305) = {305};
Line Loop(306) = {8, 31, -20, -28};
Plane Surface(306) = {306};
Line Loop(307) = {11, 30, -23, -27};
Plane Surface(307) = {307};
Line Loop(308) = {12, 33, -24, -30};
Plane Surface(308) = {308};

// Internal splitting surfaces
Line Loop(401) = {9, 29, -21, -26};
Plane Surface(401) = {401};
Line Loop(402) = {10, 32, -22, -29};
Plane Surface(402) = {402};
Line Loop(403) = {3, 29, -15, -28};
Plane Surface(403) = {403};
Line Loop(404) = {4, 30, -16, -29};
Plane Surface(404) = {404};

Surface Loop(501) = {101, 201, 301, 305, 401, 403};
Volume(501) = {501};
Surface Loop(502) = {102, 202, 302, 307, 401, 404};
Volume(502) = {502};
Surface Loop(503) = {103, 203, 303, 306, 402, 403};
Volume(503) = {503};
Surface Loop(504) = {104, 204, 304, 308, 402, 404};
Volume(504) = {504};

Transfinite Surface {101:104, 201:204, 301:308, 401:404};
Transfinite Volume {501:504};

// Boundary patches for OpenFOAM/gmshToFoam.
Physical Surface("punchLoading") = {201};
Physical Surface("topFree") = {202, 203, 204};
Physical Surface("xMin") = {305, 306};
Physical Surface("xMax") = {307, 308};
Physical Surface("yMin") = {301, 302};
Physical Surface("yMax") = {303, 304};
Physical Surface("bottom") = {101, 102, 103, 104};

Physical Volume("solid") = {501, 502, 503, 504};

Mesh 3;
