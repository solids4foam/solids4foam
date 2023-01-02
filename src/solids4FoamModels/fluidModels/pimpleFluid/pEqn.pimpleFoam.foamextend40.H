{
    p.boundaryField().updateCoeffs();

    // Prepare clean Ap without time derivative contribution and
    // without contribution from under-relaxation
    // HJ, 26/Oct/2015
    aU = HUEqn.A();

    // Store velocity under-relaxation point before using U for the flux
    // precursor
    U.storePrevIter();

    U = HUEqn.H()/aU;
    phi = (fvc::interpolate(U) & mesh.Sf());

    // Update flux on the FSI interface for Robin BCs
    robin_.updateInterface(phi, rAU);

    // Global flux balance
    adjustPhi(phi, U, p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1/aU, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve
        (
            mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

    // Correct flux due for Robin BCs according to temporal discretisation
    robin_.correctInterface(p, phi, rAU);

    // Explicitly relax pressure for momentum corrector except for last corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

#   include "movingMeshContinuityErrs.H"

    U = UUrf*
        (
            1.0/(aU + ddtUEqn.A())*
            (
                U*aU - fvc::grad(p) + ddtUEqn.H()
            )
        )
      + (1 - UUrf)*U.prevIter();
    U.correctBoundaryConditions();
}
