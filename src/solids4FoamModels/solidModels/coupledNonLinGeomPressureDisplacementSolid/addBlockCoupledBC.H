// Hack block-coupled boundary conditions: due for rewrite

forAll (D().boundaryField(), patchI)
{
    if (D().boundaryField()[patchI].blockCoupled())
    {
        Info<< "Correcting coupled BC" << endl;

        // Insert correcting fully implicit coupling coefficient

        const labelList fc = mesh().boundary()[patchI].faceCells();

        const fvPatchVectorField& Dp = D().boundaryField()[patchI];

        // Warning: hacked for nuEff in viscosity
        const scalarField nutpMagSf =
            impK_.boundaryField()[patchI]
           *mesh().magSf().boundaryField()[patchI];

        // Get boundary condition contribution to matrix diagonal
        tensorField patchDiag =
            -Dp.blockGradientInternalCoeffs()().asSquare()*nutpMagSf;

        // Get matrix diagonal
        CoeffField<vector4>::squareTypeField& blockDiag =
            DpEqn.diag().asSquare();

        forAll (fc, faceI)
        {
            blockDiag[fc[faceI]](0, 0) += patchDiag[faceI].xx();
            blockDiag[fc[faceI]](0, 1) += patchDiag[faceI].xy();
            blockDiag[fc[faceI]](0, 2) += patchDiag[faceI].xz();

            blockDiag[fc[faceI]](1, 0) += patchDiag[faceI].yx();
            blockDiag[fc[faceI]](1, 1) += patchDiag[faceI].yy();
            blockDiag[fc[faceI]](1, 2) += patchDiag[faceI].yz();

            blockDiag[fc[faceI]](2, 0) += patchDiag[faceI].zx();
            blockDiag[fc[faceI]](3, 1) += patchDiag[faceI].zy();
            blockDiag[fc[faceI]](3, 2) += patchDiag[faceI].zz();
        }

        // Get boundary condition contribution to matrix source
        vectorField patchSource =
            -Dp.blockGradientBoundaryCoeffs()*nutpMagSf;

        // Get matrix source
        Field<vector4>& blockSource = DpEqn.source();

        forAll (fc, faceI)
        {
            blockSource[fc[faceI]](0) -= patchSource[faceI](0);
            blockSource[fc[faceI]](1) -= patchSource[faceI](1);
            blockSource[fc[faceI]](2) -= patchSource[faceI](2);
        }
    }
}
