// Code by Z. Tukovic
// Correct central coefficient near boundaries
// in order to reduce jump due to laplacian term discretisation
volScalarField AUcorr("AUcorr", AU);

if (true)
{
    forAll(mesh.boundary(), patchI)
    {
        if
        (
            !U().boundaryField()[patchI].coupled()
         && !isA<wedgeFvPatchVectorField>(U().boundaryField()[patchI])
        )
        {
            const scalarField& deltaCoeffs =
                mesh.deltaCoeffs().boundaryField()[patchI];

            //const vectorField& Sf =
            //    mesh.Sf().boundaryField()[patchI];

            const scalarField& magSf =
                mesh.magSf().boundaryField()[patchI];

            //vectorField nf = Sf/magSf;

            //vectorField intCoeffsV =
            //    U().boundaryField()[patchI]
            //   .gradientInternalCoeffs();

            vectorField intCoeffsV
            (
                nu_.value()*magSf
               *U().boundaryField()[patchI].gradientInternalCoeffs()
            );

            scalarField intCoeffs(cmptAv(intCoeffsV));

            forAll(intCoeffs, faceI)
            {
                label curCell =
                    mesh.boundary()[patchI].faceCells()[faceI];

#ifdef OPENFOAMESIORFOUNDATION
                AUcorr.primitiveFieldRef()[curCell] +=
#else
                AUcorr.internalField()[curCell] +=
#endif
                    intCoeffs[faceI]/mesh.V().field()[curCell];

                //AUcorr.internalField()[curCell] +=
                //    nu_.value()*magSf[faceI]*intCoeffs[faceI]
                //   /mesh.V().field()[curCell];

#ifdef OPENFOAMESIORFOUNDATION
                AUcorr.primitiveFieldRef()[curCell] +=
#else
                AUcorr.internalField()[curCell] +=
#endif
                    nu_.value()*deltaCoeffs[faceI]*magSf[faceI]
                   /mesh.V().field()[curCell]/2.0;
            }
        }
    }
}

AUcorr.correctBoundaryConditions();
