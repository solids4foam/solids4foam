{
    p.boundaryField().updateCoeffs();

    // Prepare clean 1/a_p without time derivative and under-relaxation
    // contribution
    rAU = 1.0/HUEqn.A();

    // Calculate U from convection-diffusion matrix
    U = rAU*HUEqn.H();

    // Consistently calculate flux
    pimple.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

    // Update flux on the FSI interface for Robin BCs
    robin_.updateInterface(phi, rAU);

    // Global flux balance
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                fvc::interpolate(rAU)/pimple.aCoeff(U.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
         ==
            fvc::div(phi)
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

    // Explicitly relax pressure for momentum corrector except for last
    // corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

#   include "movingMeshContinuityErrs.H"

    // Consistently reconstruct velocity after pressure equation. Note: flux is
    // made relative inside the function
    pimple.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);
}
