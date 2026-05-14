// Calculate AD
{
    const surfaceScalarField& deltaCoeffs = mesh().deltaCoeffs();
    const scalarField& deltaCoeffsI = deltaCoeffs.internalField();

    const surfaceScalarField& magS = mesh().magSf();
    const scalarField& magSI = magS.internalField();

    const unallocLabelList& own = mesh().owner();
    const unallocLabelList& nei = mesh().neighbour();

    scalarField& ADI = AD_.internalField();
    ADI = 0;

    const surfaceScalarField& impKf = impKf_;
    // surfaceScalarField impKf =
    //     harmonic<scalar>(mesh()).interpolate(mechanical().impK());

    const scalarField& muI = impKf.internalField();

    scalar c = 1;

    forAll(deltaCoeffsI, faceI)
    {
        ADI[own[faceI]] += c*muI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
        ADI[nei[faceI]] += c*muI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
    }

    forAll(deltaCoeffs.boundaryField(), patchI)
    {
        {
            const scalarField& deltaCoeffsI =
                deltaCoeffs.boundaryField()[patchI];

            const scalarField& magSI =
                magS.boundaryField()[patchI];

            const labelList& faceCells =
                mesh().boundary()[patchI].faceCells();

            const scalarField& muI = impKf.boundaryField()[patchI];

            if
            (
                deltaCoeffs.boundaryField()[patchI].coupled() ||
                isA<wedgeFvPatchVectorField>
                (
                    DD().boundaryField()[patchI]
                )
            )
            {
                forAll(deltaCoeffsI, faceI)
                {
                    ADI[faceCells[faceI]] +=
                        c*muI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
                }
            }
            else
            {
                forAll(deltaCoeffsI, faceI)
                {
                    ADI[faceCells[faceI]] +=
                        // c*muI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
                        0.5*c*muI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
                }
            }
        }
    }

    scalar coeffAD_ = 1;

    ADI /= mesh().V().field()/coeffAD_;

    AD_.correctBoundaryConditions();
}

// Calculate ADt
// {
//     scalarField& ADtI = ADt_.internalField();
//     ADtI = 0;

//     // Add unsteady term contribution
//     if (!steadyState)
//     {
//         ADtI =
//             (1.0-backward.value())*sqr(coeff_.value())*
//             rho().internalField()*rDeltaT.value()*
//             rDeltaT.value()*mesh().V().field()

//           + backward.value()*rho().internalField()*
//             rDeltaT.value()*Cn.value()*
//             Cn.value()*rDeltaT.value()*mesh().V().field();

//         // AD_.internalField() /= sqr(coeff_.value());
//         // AD_.correctBoundaryConditions();
//     }

//     ADtI /= mesh().V().field();

//     ADt_.correctBoundaryConditions();
// }
