// Compute right hand sides
rhsLm_ = fvc::surfaceIntegrate(tC_*magSf_);

// if (angularMomentumConservation_ == "yes")
// {
//     rhsAm_ = fvc::surfaceIntegrate((xF_ ^ tC_)*magSf_);
//     am_.AMconservation(rhsLm_, rhsLm1_, rhsAm_, stage);
// }

if (runTime_.timeIndex() == 0)
{
    // Update coordinates
    // x_ += pDeltaT_*(lm_/rho_) - pDeltaT_*( (x_.prevIter() - x_.oldTime())/(deltaT_));
    // xF_ += pDeltaT_*(lmC_/rho_) - pDeltaT_*( (xF_.prevIter() - xF_.oldTime())/(deltaT_));
    xN_ += pDeltaT_*(lmN_/rho_)- pDeltaT_*( (xN_.prevIter() - xN_.oldTime())/(deltaT_));

    lm_ += pDeltaT_*rhsLm_ - pDeltaT_*dampingCoeff_*lm_ - pDeltaT_*( (lm_.prevIter() - lm_.oldTime())/(deltaT_));

    // Update deformation gradient tensor
    F_ += pDeltaT_*fvc::surfaceIntegrate((lmC_/rho_)*Sf_) - pDeltaT_*( (F_.prevIter() - F_.oldTime())/(deltaT_)) ;
}

else
{
    // Update coordinates
    // x_ += pDeltaT_*(lm_/rho_)   - pDeltaT_ * ( (3*x_.prevIter() - 4* x_.oldTime()   + x_.oldTime().oldTime())/(2*deltaT_));
    // xF_ += pDeltaT_*(lmC_/rho_) - pDeltaT_ * ( (3*xF_.prevIter() - 4* xF_.oldTime() + xF_.oldTime().oldTime())/(2*deltaT_));
    xN_ += pDeltaT_*(lmN_/rho_) - pDeltaT_ * ( (3*xN_.prevIter() - 4* xN_.oldTime() + xN_.oldTime().oldTime())/(2*deltaT_));

    lm_ += pDeltaT_*rhsLm_ - pDeltaT_*dampingCoeff_*lm_ - pDeltaT_*( (3*lm_.prevIter() - 4* lm_.oldTime() + lm_.oldTime().oldTime())/(2*deltaT_));

    // Update deformation gradient tensor
    F_ += pDeltaT_*fvc::surfaceIntegrate((lmC_/rho_)*Sf_) - pDeltaT_*( (3*F_.prevIter() - 4* F_.oldTime() + F_.oldTime().oldTime())/(2*deltaT_));
}