{
    const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh());

	// Interpolates the cell pressure to the face
	// by using pf = pf' + grad(p) dot dist (correction with the gradient)?
	// pf' is the interpolated pressure (first line above)
    if (scv.skew())
    {
        pf_ +=
        (
            scv()
          & linear<vector>(mesh()).interpolate
            (
                gradp_
            )
        );
    }

    gradp_ = fv::gaussGrad<scalar>(mesh()).gradf(pf_, "grad(p)");
    gradp_.correctBoundaryConditions();

    // gradp = fvc::grad(p);
    // gradp.correctBoundaryConditions();

    fv::gaussGrad<scalar>(mesh()).correctBoundaryConditions(p_, gradp_);
}
