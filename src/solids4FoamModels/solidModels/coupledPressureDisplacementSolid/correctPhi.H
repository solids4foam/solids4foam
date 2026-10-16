{
    forAll(D().boundaryField(), patchI)
    {
        if
        (
            D().boundaryField()[patchI].fixesValue()
        )
        {
            // Info << "Setting flux (phi) for patch: "
            //      << mesh().boundary()[patchI].name() << endl;

            phi_.boundaryField()[patchI] =
            (
                (
                    D().boundaryField()[patchI]
                )
              & mesh().Sf().boundaryField()[patchI]
            );
            // phi_.boundaryField()[patchI] =
            // (
            //     (
            //         D().boundaryField()[patchI]
            //       - D().oldTime().boundaryField()[patchI]
            //     )
            //   & mesh().Sf().boundaryField()[patchI]
            // );

            // if (nonLinear_)
            // {
            //     phi_.boundaryField()[patchI] =
            //     (
            //        -D().oldTime().boundaryField()[patchI]
            //       & mesh().Sf().boundaryField()[patchI]
            //     );
            // }
            // else
            // {
            //     phi_.boundaryField()[patchI] = 0;
            // }
        }
        else if
        (
            isA<basicSymmetryFvPatchVectorField>
            (
                D().boundaryField()[patchI]
            )
        )
        {
            // Info << "Setting zero flux (phi) for patch: "
            //      << mesh().boundary()[patchI].name() << endl;

            phi_.boundaryField()[patchI] = 0;
        }
        else if
        (
            isA<wedgeFvPatchVectorField>
            (
                D().boundaryField()[patchI]
            )
        )
        {
            // Info << "Setting zero flux (phi) for patch: "
            //      << mesh().boundary()[patchI].name() << endl;

            phi_.boundaryField()[patchI] = 0;
        }
    }
}
