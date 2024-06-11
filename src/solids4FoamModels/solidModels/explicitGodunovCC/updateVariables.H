// Strain measures
H = det(F)*op.invT(F);
J = det(F);

// Wave speeds
mech.correct(Up, Us, F);
Up_time = model.Up()/mech.stretch();

// Constitutive model
model.correct();
P = model.piola();

// Cell gradients
lmGrad = grad.gradient(lm);
grad.gradient(P, PxGrad, PyGrad, PzGrad);

// Reconstruction
grad.reconstruct(lm, lmGrad, lm_M, lm_P);
grad.reconstruct(P, PxGrad, PyGrad, PzGrad, P_M, P_P);
t_M = P_M & N;
t_P = P_P & N;

// Riemann solver
S_lm = mech.Smatrix_lm();
S_t = mech.Smatrix_t();
#include "riemannSolver.H"

// Nodal linear momentum
lmR = interpolate.surfaceToVol(lmC);
lmRgrad = grad.localGradient(lmR, lmC);
interpolate.volToPoint(lmR, lmRgrad, lmN);
#include "strongBCs.H"
lmN.correctBoundaryConditions();

// Constrained fluxes
lmC = interpolate.pointToSurface(lmN);
