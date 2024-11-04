// Compute right hand sides
rhsLm_ = fvc::surfaceIntegrate(tC_*magSf_);

if (angularMomentumConservation_ == "yes")
{
    rhsAm_ = fvc::surfaceIntegrate((xF_ ^ tC_)*magSf_);
    am_.AMconservation(rhsLm_, rhsLm1_, rhsAm_, stage);
}

// Update coordinates
x_ += deltaT_*(lm_/rho_);
xF_ += deltaT_*(lmC_/rho_);
xN_ += deltaT_*(lmN_/rho_);

D() += deltaT_ *(lm_/rho_);// -  pDeltaT_* ( (D().prevIter() - D().oldTime())/(deltaT_));
D().correctBoundaryConditions();


lm_ += deltaT_*rhsLm_ - deltaT_*dampingCoeff_*lm_;

// Update deformation gradient tensor
F_ += deltaT_*fvc::surfaceIntegrate((lmC_/rho_)*Sf_);
