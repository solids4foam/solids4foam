{
    volTensorField tmpGradDD
    (
        IOobject
        (
            "tmpGradDD",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        // zeroGradientFvPatchTensorField::typeName
        extrapolatedFvPatchTensorField::typeName
    );

    tmpGradDD = fvc::grad(DD());
    tmpGradDD.correctBoundaryConditions();
    fv::gradScheme<vector>::correctBoundaryConditions
    (
        DD(),
        tmpGradDD
    );

    gradDD() = tmpGradDD;
    gradDDf_ = fvc::interpolate(gradDD());
    const surfaceVectorField n = mesh().Sf()/mesh().magSf();
    gradDDf_ += n*fvc::snGrad(DD()) - (sqr(n) & gradDDf_);

    // gradDDf_ += n*fvc::snGrad(DD(), "snGrad(DD)ho") - (sqr(n) & gradDDf_);
}
