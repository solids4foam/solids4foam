{
    // Get matrix diagonal
    tensor4Field& blockDiag = DpEqn.diag().asSquare();

    // Grab off-diagonal coeffs
    tensor4Field& blockUpper = DpEqn.upper().asSquare();

    // Grab off-diagonal coeffs
    tensor4Field& blockLower = DpEqn.lower().asSquare();

    // Grab source
    vector4Field& source = DpEqn.source();

    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    const surfaceVectorField& Sf = mesh().Sf();
    const surfaceScalarField& magSf = mesh().magSf();
    const surfaceScalarField& deltaCoeffs = mesh().deltaCoeffs();

    surfaceTensorField fCoeff = impKf_*(Sf*Sf)*deltaCoeffs/magSf;

    if (!mesh().orthogonal())
    {
        correction_ =
            impKf_*Sf
           *(
               Sf & (mesh().correctionVectors() & gradDf_)
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

     
        source[own[faceI]](0) += correction_[faceI].x();
        source[own[faceI]](1) += correction_[faceI].y();
        source[own[faceI]](2) += correction_[faceI].z();

        source[nei[faceI]](0) -= correction_[faceI].x();
        source[nei[faceI]](1) -= correction_[faceI].y();
        source[nei[faceI]](2) -= correction_[faceI].z();
    }

    forAll(fCoeff.boundaryField(), patchI)
    {
        const labelList& faceCells = mesh().boundary()[patchI].faceCells();

        const vectorField nf =
            Sf.boundaryField()[patchI]
           /magSf.boundaryField()[patchI];

        const scalarField& pMagSf = magSf.boundaryField()[patchI];

        tensorField fCoeff =
            impK_.boundaryField()[patchI]*pMagSf*(nf*nf)
           *(
                D().boundaryField()[patchI].gradientInternalCoeffs()
              & vector::one
            );

        vectorField sourceContrib =
            impK_.boundaryField()[patchI]*pMagSf
           *(
                (nf*nf)
              & D().boundaryField()[patchI].gradientBoundaryCoeffs()
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
