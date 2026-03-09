if (nonLinear_)
// if (false)
{
    surfaceTensorField Ff = I + gradDf_.T();

    surfaceTensorField Finvf = inv(Ff);

    surfaceScalarField Jf = det(Ff);

    if (gMax(nuCoeff_.internalField()) > 0.5-SMALL)
    {
        // Incompressible case
        Jf = 1.0;
    }

    surfaceVectorField curSf = Jf*(Finvf.T() & mesh().Sf());

    surfaceVectorField nonLinForce = (curSf & sigmaf_);

    force_.oldTime();
    force_ = nonLinForce;

    // forAll(DD().boundaryField(), patchI)
    // {
    //     if
    //     (
    //         isA<tractionExtrapolatedDisplacementFvPatchVectorField>
    //         (
    //             DD().boundaryField()[patchI]
    //         )
    //     )
    //     {
    //         tractionExtrapolatedDisplacementFvPatchVectorField ted =
    //             refCast<tractionExtrapolatedDisplacementFvPatchVectorField>
    //             (
    //                 DD().boundaryField()[patchI]
    //             );

    //         force_.boundaryField()[patchI] =
    //         (
    //             ted.traction()*mag(curSf.boundaryField()[patchI])
    //           - ted.pressure()*curSf.boundaryField()[patchI]
    //         );
    //     }
    // }
}
else
{
    // Info << "Calculate linear face forces" << endl;

    surfaceVectorField nf = mesh().Sf()/mesh().magSf();
    // surfaceVectorField tGradDDn = ((I-sqr(nf)) & (gradDDf_ & nf));
    surfaceVectorField nGradDDn = ((nf*nf) & fvc::snGrad(DD()));

    surfaceVectorField linearForceIncrement =
        impKf_*fvc::snGrad(DD())*mesh().magSf()
      + impKf_*nGradDDn*mesh().magSf()
      + impKf_*tGradDDn_*mesh().magSf()
      - fvc::interpolate(nuCoeff_)*Dpf_*mesh().Sf();

    force_ = force_.oldTime() + linearForceIncrement;

    forAll(DD().boundaryField(), patchI)
    {
        if
        (
            isA<tractionPressureDisplacementFvPatchVectorField>
            (
                DD().boundaryField()[patchI]
            )
        )
        {
            tractionPressureDisplacementFvPatchVectorField ted =
                refCast<tractionPressureDisplacementFvPatchVectorField>
                (
                    DD().boundaryField()[patchI]
                );

            force_.boundaryField()[patchI] =
            (
                ted.traction()*mesh().magSf().boundaryField()[patchI]
              - ted.pressure()*mesh().Sf().boundaryField()[patchI]
            );
        }
    }

    // surfaceVectorField nf = mesh().Sf()/mesh().magSf();
    // surfaceVectorField tGradDn = ((I-sqr(nf)) & (gradDf_ & nf));

    // surfaceVectorField nGradDn = ((nf*nf)&fvc::snGrad(D()));

    // surfaceVectorField linearForce =
    //     impKf_*fvc::snGrad(D())*mesh().magSf()
    //   + impKf_*nGradDn*mesh().magSf()
    //   + impKf_*tGradDn*mesh().magSf()
    //   - fvc::interpolate(nuCoeff_)*pf_*mesh().Sf();

    // force_ = linearForce;
}
