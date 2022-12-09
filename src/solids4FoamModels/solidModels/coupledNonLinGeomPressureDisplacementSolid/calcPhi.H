{
    surfaceVectorField HbyAf =
        linear<vector>(mesh()).interpolate(HbyA_);

    const skewCorrectionVectors& scv =
        skewCorrectionVectors::New(mesh());

    if (scv.skew())
    {
        volTensorField gradHbyA = fvc::grad(HbyA_);
        HbyAf += (scv() & linear<tensor>(mesh()).interpolate(gradHbyA));
    }

    HbyAf -= linear<vector>(mesh()).interpolate(D());

    phi_ = (HbyAf & mesh().Sf());
}
