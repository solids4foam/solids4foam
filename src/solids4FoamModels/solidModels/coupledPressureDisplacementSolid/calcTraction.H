if (nonLinear_)
{
    surfaceVectorField linearForceIncrement =
        impKf_*fvc::snGrad(DD())*mesh().magSf()
      + impKf_*nGradDDn_*mesh().magSf()
      + impKf_*tGradDDn_*mesh().magSf()
      - fvc::interpolate(nuCoeff_)*Dpf_*curSf_.oldTime();

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

    // Non-linear force increment
    surfaceVectorField nonLinForceIncrement =
        nonLinForce - force_.oldTime();

    nonLinForceCorrection_ =
        nonLinForceIncrement
      - linearForceIncrement;

    nonLinForceCorrection_.relax();

    // nonLinForceCorrection_ =
    //     dimensionedVector("zero", dimForce, vector::zero);
}
else
{
    nonLinForceCorrection_ =
        dimensionedVector("zero", dimForce, vector::zero);
}
