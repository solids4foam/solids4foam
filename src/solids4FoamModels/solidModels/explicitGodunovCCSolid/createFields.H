Info << "Creating mesh parameters ..." << endl;

// Operations class
operations op(mesh());

// Point mesh
pointMesh pMesh(mesh());

// Material face area
const surfaceScalarField& magSf = mesh().magSf();

// Material face area normal vector
const surfaceVectorField& Sf = mesh().Sf();

// Minimum edge length
const dimensionedScalar& h = op.minimumEdgeLength();

// Boundary patches
const polyBoundaryMesh& bm = mesh().boundaryMesh();
const label& symmetricPatchID = bm.findPatchID("symmetric");
const label& symmetricXpatchID = bm.findPatchID("symmetricX");
const label& symmetricYpatchID = bm.findPatchID("symmetricY");
const label& symmetricZpatchID = bm.findPatchID("symmetricZ");


Info << "Creating mesh coordinate fields ..." << endl;

// Material cell center coordinates
const volVectorField& C = mesh().C();

// Spatial cell center coordinates
volVectorField x
(
    IOobject("x", mesh()),
    C
);

// Spatial nodal coordinates
pointVectorField xN
(
    IOobject("xN", mesh()),
    pMesh,
    dimensionedVector("xN", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
xN.primitiveFieldRef() = mesh().points();

// Material nodal coordinates
pointVectorField XN = xN;

// Spatial face center coordinates
surfaceVectorField xF = mesh().Cf();


Info << "Creating mesh normal fields ..." << endl;

// Material normals
const surfaceVectorField N = (Sf/mesh().magSf()).ref();

// Spatial normals
surfaceVectorField n = N;


Info << "Creating linear momentum fields ..." << endl;

// Cell linear momentum
volVectorField lm
(
    IOobject
    (
        "lm",
        runTime.timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh(),
    dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
);

// Nodal linear momentum
pointVectorField lmN
(
    IOobject
    (
        "lmN",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh
);


Info << "Creating strain measure fields ..." << endl;

// Deformation gradient tensor
volTensorField F
(
    IOobject("F", mesh()),
    mesh(),
    tensor::I
);

// Cofactor of deformation
volTensorField H
(
    IOobject("H", mesh()),
    det(F)*op.invT(F)
);

// Jacobian of deformation
volScalarField J
(
    IOobject("J", mesh()),
    det(F)
);


Info << "Creating constitutive model ..." << endl;

// // Solid model class
solidMaterialModel model(F, mechanicalProperties_);
model.correct();

// Density
const dimensionedScalar& rho_ = model.density();

// Pressure
volScalarField p_ = model.pressure();

// First Piola Kirchhoff stress tensor
volTensorField P = model.piola();
volVectorField Px = op.decomposeTensorX(P);
volVectorField Py = op.decomposeTensorY(P);
volVectorField Pz = op.decomposeTensorZ(P);


Info << "Creating fields for wave speeds ..." << endl;

// Continuum mechanics class
mechanics mech(F, controlDict);

// Longitudinal wave speed
volScalarField Up
(
    IOobject("Up", mesh()),
    mesh(),
    model.Up()/beta
);

// Shear wave speed
volScalarField Us
(
    IOobject ("Us", mesh()),
    mesh(),
    model.Us()*beta
);

mech.correct(Up, Us, F);



// Wave speed for time increment
volScalarField Up_time = (model.Up()/mech.stretch()).ref();


Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh());

// Gradient of cell linear momentum
volTensorField lmGrad = grad.gradient(lm);

// Gradient of first Piola Kirchhoff stress tensor
volTensorField PxGrad = grad.gradient(Px);
volTensorField PyGrad = grad.gradient(Py);
volTensorField PzGrad = grad.gradient(Pz);


Info << "Creating fields for reconstruction ..." << endl;

// Reconstruction of linear momentum
surfaceVectorField lm_M
(
    IOobject("lm_M", mesh()),
    mesh(),
    dimensionedVector("lm_M", lm.dimensions(), vector::zero)
);
surfaceVectorField lm_P = lm_M;
grad.reconstruct(lm, lmGrad, lm_M, lm_P);

// Reconstruction of PK1 stresses
surfaceTensorField P_M
(
    IOobject("P_M", mesh()),
    mesh(),
    dimensionedTensor("P_M", P.dimensions(), tensor::zero)
);
surfaceTensorField P_P = P_M;
grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P);

// Reconstruction of traction
surfaceVectorField t_M
(
    IOobject("t_M", mesh()),
    P_M & N
);
surfaceVectorField t_P = (P_P & N).ref();


Info << "Creating fields for riemann solver ..." << endl;

surfaceTensorField S_lm = mech.Smatrix_lm();
surfaceTensorField S_t = mech.Smatrix_t();
S_lm.oriented() = false;
S_t.oriented() = false;
t_M.oriented() = false;
t_P.oriented() = false;
lm_P.oriented() = false;
lm_M.oriented() = false;

// Contact traction
surfaceVectorField tC = (0.5*(t_M + t_P) + (0.5*S_lm & (lm_P - lm_M))).ref();

volVectorField t_b
(
    IOobject
    (
        "t_b",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh()
);

// Contact linear momentum
surfaceVectorField lmC(0.5*(lm_M + lm_P) + (0.5*S_t & (t_P - t_M)));

volVectorField lm_b
(
    IOobject
    (
        "lm_b",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh()
);

// // Riemann solver
#include "riemannSolver.H"


Info << "Creating fields for the constrained procedure ..." << endl;

// Constrained class
interpolationSchemes interpolate(mesh());

// Cell averaged linear momentum
volVectorField lmR = interpolate.surfaceToVol(lmC);

// Local gradient of cell averaged linear momentum
volTensorField lmRgrad = grad.localGradient(lmR, lmC);

// Constrained fluxes
interpolate.volToPoint(lmR, lmRgrad, lmN);
#include "strongBCs.H"
lmN.correctBoundaryConditions();
lmC = interpolate.pointToSurface(lmN);


Info << "Creating fields for angular momentum ..." << endl;

// Angular momentum class
angularMomentum am(mesh(), mechanicalProperties_);

// RHS of linear momentum equation
volVectorField rhsLm
(
    IOobject("rhsLm", mesh()),
    mesh(),
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

// RHS of linear momentum equation (stage 1)
volVectorField rhsLm1 = rhsLm;

// RHS of angular momentum equation
volVectorField rhsAm
(
    IOobject("rhsAm", mesh()),
    mesh(),
    dimensionedVector("rhsAm", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);


Info << "Creating fields for post-processing ..." << endl;

// Nodal displacements
pointVectorField uN
(
    IOobject
    (
        "uN",
        runTime.timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("uN", dimLength, vector::zero)
);
uN.write();


Info << "Creating variables for time ..." << endl;

// Time increment
dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());

// Runge-Kutta stage
scalarList RKstages(2);
RKstages[0] = 0;
RKstages[1] = 1;


Info << "Printing data ..." << endl;

// Print material properties
model.printMaterialProperties();

// Print global linear and angular momentum
am.printGlobalMomentum(lm,x);

// Print centroid of geometry
mech.printCentroid();
