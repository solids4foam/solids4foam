    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U));

    // Get under-relaxation factor
    scalar UUrf =
        mesh.solutionDict().equationRelaxationFactor(U.select(pimple.finalIter()));

    if (pimple.momentumPredictor())
    {
        // Solve momentum predictor
        solve
        (
            ddtUEqn
          + relax(HUEqn, UUrf)
         ==
          - fvc::grad(p),
            mesh.solutionDict().solver((U.select(pimple.finalIter())))
        );
    }
