// Tetrahedral mesh of the silicone flap of the Hessenthaler et al. (2017)
// FSI experiment: an 11 x 2 x 65 mm brick clamped at z = 0.
//
// lc is the target tetrahedron size in metres; override it with
//     gmsh -3 -setnumber lc 0.0005 flap.geo

SetFactory("OpenCASCADE");

Mesh.MshFileVersion = 2.2;
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;

DefineConstant[ lc = {0.001, Name "Parameters/lc"} ];

Box(1) = {-0.0055, -0.001, 0, 0.011, 0.002, 0.065};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// OpenCASCADE box faces: 1 x-min, 2 x-max, 3 y-min, 4 y-max, 5 z-min, 6 z-max
Physical Surface("fixed") = {5};
Physical Surface("interface") = {1, 2, 3, 4, 6};
Physical Volume("solid") = {1};
