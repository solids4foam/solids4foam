// Calculate delta
{
    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();
    const vectorField& C = mesh().C().internalField();

    forAll(delta_.internalField(), faceI)
    {
        delta_.internalField()[faceI] =
            (C[nei[faceI]] - C[own[faceI]]);
    }

    forAll(delta_.boundaryField(), patchI)
    {
        delta_.boundaryField()[patchI] =
            mesh().boundary()[patchI].delta();
        // if (delta_.boundaryField()[patchI].coupled())
        // {
        //     delta_.boundaryField()[patchI] =
        //         mesh().boundary()[patchI].delta();
        // }
        // else
        // {
        //     delta_.boundaryField()[patchI] =
        //         mesh().boundary()[patchI].nf()/
        //         mesh().deltaCoeffs().boundaryField()[patchI];
        // }
    }
}
