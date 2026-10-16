// Add laplacian of normal displacement component
{
    // Get matrix diagonal
    tensor4Field& blockDiag = UpEqn.diag().asSquare();

    // Grab off-diagonal coeffs
    tensor4Field& blockUpper = UpEqn.upper().asSquare();

    // Grab off-diagonal coeffs
    tensor4Field& blockLower = UpEqn.lower().asSquare();

    // Grab source
    vector4Field& source = UpEqn.source();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

    surfaceTensorField fCoeff =
        mu.value()*(Sf*Sf)*deltaCoeffs/magSf;

    surfaceVectorField correction
    (
        IOobject
        (
            "correction",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            mu.dimensions()*dimArea,
            vector::zero
        )
    );

    if (!mesh.orthogonal())
    {
        correction =
            mu*Sf
           *(
               Sf & (mesh.correctionVectors() & fvc::interpolate(gradU))
            )/magSf;
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

    forAll(fCoeff.boundaryField(), patchI)
    {
        const labelList& faceCells = mesh.boundary()[patchI].faceCells();

        const vectorField nf =
            Sf.boundaryField()[patchI]
           /magSf.boundaryField()[patchI];

        const scalarField& pMagSf = magSf.boundaryField()[patchI];

        tensorField fCoeff =
            mu.value()*pMagSf*(nf*nf)
           *(
                U.boundaryField()[patchI].gradientInternalCoeffs()
              & vector::one
            );

        vectorField sourceContrib =
            mu.value()*pMagSf
           *(
                (nf*nf)
              & U.boundaryField()[patchI].gradientBoundaryCoeffs()
            );

        forAll(faceCells, faceI)
        {
            const label& curCell = faceCells[faceI];

            blockDiag[curCell](0,0) -= fCoeff[faceI].xx();
            blockDiag[curCell](0,1) -= fCoeff[faceI].xy();
            blockDiag[curCell](0,2) -= fCoeff[faceI].xz();

            blockDiag[curCell](1,0) -= fCoeff[faceI].yx();
            blockDiag[curCell](1,1) -= fCoeff[faceI].yy();
            blockDiag[curCell](1,2) -= fCoeff[faceI].yz();

            blockDiag[curCell](2,0) -= fCoeff[faceI].zx();
            blockDiag[curCell](2,1) -= fCoeff[faceI].zy();
            blockDiag[curCell](2,2) -= fCoeff[faceI].zz();

            source[curCell](0) += sourceContrib[faceI].x();
            source[curCell](1) += sourceContrib[faceI].y();
            source[curCell](2) += sourceContrib[faceI].z();
        }
    }
}
