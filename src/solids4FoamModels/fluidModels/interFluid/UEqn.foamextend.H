    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    // Named operands: a tmp operand with the type of the result would be reused
    // as the storage of the inner product, which aliases an input of the loop
    // and gives wrong results with some -O3 builds (see
    // pimpleFluid::patchViscousForce)
    const volTensorField gradUvol(fvc::grad(U));
    const volVectorField gradMuEffVol(fvc::grad(muEff));

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (gradUvol & gradMuEffVol)
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(pd)
                )*mesh.magSf()
            )
        );
    }
