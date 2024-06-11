// Compute right hand sides
rhsLm = fvc::surfaceIntegrate(tC*magSf);

if (angularMomentumConservation == "yes")
{
    rhsAm = fvc::surfaceIntegrate((xF ^ tC)*magSf);
    am.AMconservation(rhsLm, rhsLm1, rhsAm, stage);
}

// Update coordinates
x += deltaT*(lm/rho);
xF += deltaT*(lmC/rho);
xN += deltaT*(lmN/rho);

// Update linear momentum
lm += deltaT*rhsLm - deltaT*dampingCoeff*lm;

// Update deformation gradient tensor
F += deltaT*fvc::surfaceIntegrate((lmC/rho)*Sf);
