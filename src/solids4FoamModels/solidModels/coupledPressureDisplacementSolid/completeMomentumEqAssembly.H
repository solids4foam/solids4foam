// Add laplacian of normal displacement component
{
    // Get matrix diagonal
    tensor4Field& blockDiag = DDDpEqn.diag().asSquare();

    // Grab off-diagonal coeffs
    tensor4Field& blockUpper = DDDpEqn.upper().asSquare();

    // Grab off-diagonal coeffs
    tensor4Field& blockLower = DDDpEqn.lower().asSquare();

    // Grab source
    vector4Field& source = DDDpEqn.source();

    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    const surfaceVectorField& Sf = mesh().Sf();
    const surfaceScalarField& magSf = mesh().magSf();
    const surfaceScalarField& deltaCoeffs = mesh().deltaCoeffs();

    const surfaceScalarField& impKf = impKf_;
    // surfaceScalarField impKf =
    //     harmonic<scalar>(mesh()).interpolate(impK_);

    surfaceTensorField fCoeff =
        impKf*(Sf*Sf)*deltaCoeffs/magSf;

    surfaceVectorField correction
    (
        IOobject
        (
            "correction",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            impKf.dimensions()*dimArea,
            vector::zero
        )
    );

    if (!mesh().orthogonal())
    {
        surfaceVectorField corr =
            fv::snGradScheme<vector>::New
            (
                mesh(),
                mesh().schemesDict().snGradScheme("snGrad(" + DD().name() + ")")
            )().correction(DD());

        correction = impKf*Sf*(Sf & corr)/magSf;
    }

    forAll(own, faceI)
    {
        blockUpper[faceI](0,0) -= fCoeff[faceI].xx();
        blockUpper[faceI](0,1) -= fCoeff[faceI].xy();
        blockUpper[faceI](0,2) -= fCoeff[faceI].xz();

        blockUpper[faceI](1,0) -= fCoeff[faceI].yx();
        blockUpper[faceI](1,1) -= fCoeff[faceI].yy();
        blockUpper[faceI](1,2) -= fCoeff[faceI].yz();

        blockUpper[faceI](2,0) -= fCoeff[faceI].zx();
        blockUpper[faceI](2,1) -= fCoeff[faceI].zy();
        blockUpper[faceI](2,2) -= fCoeff[faceI].zz();


        blockLower[faceI](0,0) -= fCoeff[faceI].xx();
        blockLower[faceI](0,1) -= fCoeff[faceI].xy();
        blockLower[faceI](0,2) -= fCoeff[faceI].xz();

        blockLower[faceI](1,0) -= fCoeff[faceI].yx();
        blockLower[faceI](1,1) -= fCoeff[faceI].yy();
        blockLower[faceI](1,2) -= fCoeff[faceI].yz();

        blockLower[faceI](2,0) -= fCoeff[faceI].zx();
        blockLower[faceI](2,1) -= fCoeff[faceI].zy();
        blockLower[faceI](2,2) -= fCoeff[faceI].zz();


        blockDiag[own[faceI]](0,0) += fCoeff[faceI].xx();
        blockDiag[own[faceI]](0,1) += fCoeff[faceI].xy();
        blockDiag[own[faceI]](0,2) += fCoeff[faceI].xz();

        blockDiag[own[faceI]](1,0) += fCoeff[faceI].yx();
        blockDiag[own[faceI]](1,1) += fCoeff[faceI].yy();
        blockDiag[own[faceI]](1,2) += fCoeff[faceI].yz();

        blockDiag[own[faceI]](2,0) += fCoeff[faceI].zx();
        blockDiag[own[faceI]](2,1) += fCoeff[faceI].zy();
        blockDiag[own[faceI]](2,2) += fCoeff[faceI].zz();


        blockDiag[nei[faceI]](0,0) += fCoeff[faceI].xx();
        blockDiag[nei[faceI]](0,1) += fCoeff[faceI].xy();
        blockDiag[nei[faceI]](0,2) += fCoeff[faceI].xz();

        blockDiag[nei[faceI]](1,0) += fCoeff[faceI].yx();
        blockDiag[nei[faceI]](1,1) += fCoeff[faceI].yy();
        blockDiag[nei[faceI]](1,2) += fCoeff[faceI].yz();

        blockDiag[nei[faceI]](2,0) += fCoeff[faceI].zx();
        blockDiag[nei[faceI]](2,1) += fCoeff[faceI].zy();
        blockDiag[nei[faceI]](2,2) += fCoeff[faceI].zz();


        source[own[faceI]](0) += correction[faceI].x();
        source[own[faceI]](1) += correction[faceI].y();
        source[own[faceI]](2) += correction[faceI].z();

        source[nei[faceI]](0) -= correction[faceI].x();
        source[nei[faceI]](1) -= correction[faceI].y();
        source[nei[faceI]](2) -= correction[faceI].z();
    }

    forAll(correction.boundaryField(), patchI)
    {
        const labelList& faceCells =
            mesh().boundary()[patchI].faceCells();

        const vectorField nf =
            Sf.boundaryField()[patchI]
           /magSf.boundaryField()[patchI];

        const scalarField& pMagSf =
            magSf.boundaryField()[patchI];

        const vectorField& pCorrection =
            correction.boundaryField()[patchI];

        vectorField bouCoeffsSourceContrib(pCorrection.size(), vector::zero);
        tensorField internalCoeffs(pCorrection.size(), tensor::zero);

        if
        (
            DD().boundaryField()[patchI].coupled()
        )
        {
            tensor4Field& cu = DDDpEqn.coupleUpper()[patchI].asSquare();
            tensor4Field& cl = DDDpEqn.coupleLower()[patchI].asSquare();

            tensorField boundaryCoeffs =
               -impKf.boundaryField()[patchI]*pMagSf*
                (
                    (nf*nf)
                  & diag(DD().boundaryField()[patchI].gradientBoundaryCoeffs())
                );

            bouCoeffsSourceContrib = vector::zero;

            internalCoeffs =
                impKf.boundaryField()[patchI]*pMagSf*
                (
                    (nf*nf)
                  & diag(DD().boundaryField()[patchI].gradientInternalCoeffs())
                );

            forAll(cu, faceI)
            {
                cu[faceI](0,0) -= boundaryCoeffs[faceI].xx();
                cu[faceI](0,1) -= boundaryCoeffs[faceI].xy();
                cu[faceI](0,2) -= boundaryCoeffs[faceI].xz();

                cu[faceI](1,0) -= boundaryCoeffs[faceI].yx();
                cu[faceI](1,1) -= boundaryCoeffs[faceI].yy();
                cu[faceI](1,2) -= boundaryCoeffs[faceI].yz();

                cu[faceI](2,0) -= boundaryCoeffs[faceI].zx();
                cu[faceI](2,1) -= boundaryCoeffs[faceI].zy();
                cu[faceI](2,2) -= boundaryCoeffs[faceI].zz();


                cl[faceI](0,0) -= internalCoeffs[faceI].xx();
                cl[faceI](0,1) -= internalCoeffs[faceI].xy();
                cl[faceI](0,2) -= internalCoeffs[faceI].xz();

                cl[faceI](1,0) -= internalCoeffs[faceI].yx();
                cl[faceI](1,1) -= internalCoeffs[faceI].yy();
                cl[faceI](1,2) -= internalCoeffs[faceI].yz();

                cl[faceI](2,0) -= internalCoeffs[faceI].zx();
                cl[faceI](2,1) -= internalCoeffs[faceI].zy();
                cl[faceI](2,2) -= internalCoeffs[faceI].zz();
            }
        }
        else if
        (
            isA<tractionPressureDisplacementFvPatchVectorField>
            (
                DD().boundaryField()[patchI]
            )
            // isA<tractionExtrapolatedDisplacementFvPatchVectorField>
            // (
            //     DD().boundaryField()[patchI]
            // )
        )
        {
            // Already added by laplacian term (?)
            bouCoeffsSourceContrib = //vector::zero;
               -impKf.boundaryField()[patchI]*pMagSf
               *nGradDDn_.boundaryField()[patchI];
               // *(
               //     (nf*nf)
               //   & DD().boundaryField()[patchI].gradientBoundaryCoeffs()
               //  );
            internalCoeffs = tensor::zero;

            // Info << mesh().boundary()[patchI].name()
            //      << ", " << average(nGradDDn_.boundaryField()[patchI]) << endl;
        }
        else if
        (
            isA<fixedValueFvPatchVectorField>
            (
                DD().boundaryField()[patchI]
            )
        )
        {
            bouCoeffsSourceContrib =
               -impKf.boundaryField()[patchI]*pMagSf
               *(
                   (nf*nf)
                 & DD().boundaryField()[patchI].gradientBoundaryCoeffs()
                );

            internalCoeffs =
                impKf.boundaryField()[patchI]*pMagSf*
                (
                    (nf*nf)
                  & diag(DD().boundaryField()[patchI].gradientInternalCoeffs())
                );
        }
        else if
        (
            isA<wedgeFvPatchVectorField>
            (
                DD().boundaryField()[patchI]
            )
        )
        {
            internalCoeffs =
                impKf.boundaryField()[patchI]*pMagSf*
                (
                    (nf*nf)
                  & diag(DD().boundaryField()[patchI].gradientInternalCoeffs())
                );

            bouCoeffsSourceContrib =
               -impKf.boundaryField()[patchI]*pMagSf
               *(
                   (nf*nf)
                 & DD().boundaryField()[patchI].gradientBoundaryCoeffs()
                );
        }
        else if
        (
            isA<symmetryFvPatchVectorField>
            (
                DD().boundaryField()[patchI]
            )
        )
        {
            bouCoeffsSourceContrib =
               -impKf.boundaryField()[patchI]*pMagSf
               *(
                   (nf*nf)
                 & DD().boundaryField()[patchI].gradientBoundaryCoeffs()
                );

            internalCoeffs =
                impKf.boundaryField()[patchI]*pMagSf*
                (
                    (nf*nf)
                  & diag(DD().boundaryField()[patchI].gradientInternalCoeffs())
                );
        }
        else if
        (
            DD().boundaryField()[patchI].blockCoupled()
        )
        {
            bouCoeffsSourceContrib =
               -impKf.boundaryField()[patchI]*pMagSf
               *(
                   (nf*nf)
                 & DD().boundaryField()[patchI].blockGradientBoundaryCoeffs()
                );

            internalCoeffs =
                impKf.boundaryField()[patchI]*pMagSf*
                (
                    (nf*nf)
                  & DD().boundaryField()[patchI].blockGradientInternalCoeffs()()
                   .asSquare()
                );
        }
        else if
        (
            isA<emptyFvPatchVectorField>
            (
                DD().boundaryField()[patchI]
            )
        )
        {
            bouCoeffsSourceContrib = vector::zero;
            internalCoeffs = tensor::zero;
        }
        else
        {
            Info << "****** BC not implemented ******" << endl;
        }

        forAll(faceCells, faceI)
        {
            const label& curCell = faceCells[faceI];

            blockDiag[curCell](0,0) -= internalCoeffs[faceI].xx();
            blockDiag[curCell](0,1) -= internalCoeffs[faceI].xy();
            blockDiag[curCell](0,2) -= internalCoeffs[faceI].xz();

            blockDiag[curCell](1,0) -= internalCoeffs[faceI].yx();
            blockDiag[curCell](1,1) -= internalCoeffs[faceI].yy();
            blockDiag[curCell](1,2) -= internalCoeffs[faceI].yz();

            blockDiag[curCell](2,0) -= internalCoeffs[faceI].zx();
            blockDiag[curCell](2,1) -= internalCoeffs[faceI].zy();
            blockDiag[curCell](2,2) -= internalCoeffs[faceI].zz();

            source[curCell](0) -= bouCoeffsSourceContrib[faceI].x();
            source[curCell](1) -= bouCoeffsSourceContrib[faceI].y();
            source[curCell](2) -= bouCoeffsSourceContrib[faceI].z();

            source[curCell](0) += pCorrection[faceI].x();
            source[curCell](1) += pCorrection[faceI].y();
            source[curCell](2) += pCorrection[faceI].z();
        }
    }
}
