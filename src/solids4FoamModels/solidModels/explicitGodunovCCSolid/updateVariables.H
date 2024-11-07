// Strain measures
H_ = det(F_)*op_.invT(F_);
J_ = det(F_);

// Wave speeds
mech_.correct(Up_, Us_, F_);
Up_time_ = model_.Up()/mech_.stretch();

// Constitutive model
model_.correct();
P_ = model_.piola();

// Cell gradients
lmGrad_ = grad_.gradient(lm_);
grad_.gradient(P_, PxGrad_, PyGrad_, PzGrad_);

// Reconstruction
grad_.reconstruct(lm_, lmGrad_, lm_M_, lm_P_);
grad_.reconstruct(P_, PxGrad_, PyGrad_, PzGrad_, P_M_, P_P_);
t_M_ = P_M_ & N_;
t_P_ = P_P_ & N_;

// Riemann solver
S_lm_ = mech_.Smatrix_lm();
S_t_ = mech_.Smatrix_t();
#include "riemannSolver.H"

// Nodal linear momentum
lmR_ = interpolate_.surfaceToVol(lmC_);
lmRgrad_ = grad_.localGradient(lmR_, lmC_);
interpolate_.volToPoint(lmR_, lmRgrad_, lmN_);
#include "strongBCs.H" 
lmN_.correctBoundaryConditions();

// Constrained fluxes
lmC_ = interpolate_.pointToSurface(lmN_);
