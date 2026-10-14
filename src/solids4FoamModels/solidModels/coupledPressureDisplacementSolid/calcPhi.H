{
    surfaceVectorField gradDpf =
        linear<vector>(mesh()).interpolate(gradDp_);

    phi_ = rADf*((delta_ & gradDpf)/mag(delta_))*mesh().magSf();

    surfaceVectorField HbyAf
    (
        IOobject
        (
            "HbyAf",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    );

    // Add consistency term
    if (consistentRhieChow_ && !steadyState)
    {
        Info << "Applying Rhie-Chow consistency correction" << endl;

        // Euler (and trapezoidal)
        // if (runTime().timeIndex()>0)
        {
            HbyAf +=
                (1-backward)*
                linear<scalar>(mesh()).interpolate(rAD*ADt_)*
                runTime().deltaT()*
                (
                    Uf_.oldTime()
                  - interpUf().oldTime()
                );
        }

        // Backward
        // if (runTime().timeIndex()>2)
        {
            HbyAf -=
                backward*
                linear<scalar>(mesh()).interpolate(rAD*ADt_)*
                ((Cn+Co)/Cn)*
                (
                    Df_.oldTime()
                  - interpDf().oldTime()
                );

            HbyAf -=
                backward*
                linear<scalar>(mesh()).interpolate(rAD*ADt_)*
                (Coo/Cn)*
                (
                    Df_.oldTime().oldTime()
                  - interpDf().oldTime().oldTime()
                );

            HbyAf -=
                backward*
                linear<scalar>(mesh()).interpolate(rAD*ADt_)*
                (Co/sqr(Cn))*runTime().deltaT()*
                (
                    Uf_.oldTime()
                  - interpUf().oldTime()
                );

            HbyAf -=
                backward*linear<scalar>(mesh()).interpolate(rAD*ADt_)*
                (Coo/sqr(Cn))*runTime().deltaT()*
                (
                    Uf_.oldTime().oldTime()
                  - interpUf().oldTime().oldTime()
                );
        }
    }

    // Initialize field to zero befor adding skewness correction
    DDf_ = dimensionedVector("0", dimLength, vector::zero);

    // if (skewCorr_)
    {
        const skewCorrectionVectors& scv =
            skewCorrectionVectors::New(mesh());

        if (scv.skew())
        {
            surfaceTensorField gradDDf =
                linear<tensor>(mesh()).interpolate
                (
                    fvc::grad(DD())
                );

            // Add skewness correction
            DDf_ = (scv() & gradDDf);

            HbyAf += DDf_;
        }
    }

    phi_ += (HbyAf & mesh().Sf());

    if (true)
    {
        forAll(DD().boundaryField(), patchI)
        {
            if
            (
                DD().boundaryField()[patchI].fixesValue()
            )
            {
                phi_.boundaryField()[patchI] = 0; // added through div(DD) therm
                // (
                //     DD().boundaryField()[patchI]
                //   & mesh().Sf().boundaryField()[patchI]
                // );
            }
            else if
            (
                isA<basicSymmetryFvPatchVectorField>
                (
                    DD().boundaryField()[patchI]
                )
            )
            {
                // Info << mesh().boundary()[patchI].name() << endl;

                phi_.boundaryField()[patchI] = 0;
                // phi0_.boundaryField()[patchI] = 0;
            }
        }
    }
}
