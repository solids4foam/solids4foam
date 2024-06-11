// Strong boundary conditions for nodal linear momentum

// Symmetric boundary patch
if (symmetricPatchID != -1)
{
    forAll(mesh.boundary()[symmetricPatchID], facei)
    {
        const label& faceID =
            mesh.boundary()[symmetricPatchID].start() + facei;

        forAll(mesh.faces()[faceID], node)
        {
            const label& nodeID = mesh.faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N.boundaryField()[symmetricPatchID][facei]
               *N.boundaryField()[symmetricPatchID][facei]);

            lmN[nodeID] = iMnCn_ & lmN[nodeID];
        }
    }
}

// SymmetricX boundary patch
if (symmetricXpatchID != -1)
{
    forAll(mesh.boundary()[symmetricXpatchID], facei)
    {
        const label& faceID =
            mesh.boundary()[symmetricXpatchID].start() + facei;

        forAll(mesh.faces()[faceID], node)
        {
            const label& nodeID = mesh.faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N.boundaryField()[symmetricXpatchID][facei]
               *N.boundaryField()[symmetricXpatchID][facei]);

            lmN[nodeID] = iMnCn_ & lmN[nodeID];
        }
    }
}

// SymmetricY boundary patch
if (symmetricYpatchID != -1)
{
    forAll(mesh.boundary()[symmetricYpatchID], facei)
    {
        const label& faceID =
            mesh.boundary()[symmetricYpatchID].start() + facei;

        forAll(mesh.faces()[faceID], node)
        {
            const label& nodeID = mesh.faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N.boundaryField()[symmetricYpatchID][facei]
               *N.boundaryField()[symmetricYpatchID][facei]);

            lmN[nodeID] = iMnCn_ & lmN[nodeID];
        }
    }
}

// SymmetricZ boundary patch
if (symmetricZpatchID != -1)
{
    forAll(mesh.boundary()[symmetricZpatchID], facei)
    {
        const label& faceID =
            mesh.boundary()[symmetricZpatchID].start() + facei;

        forAll(mesh.faces()[faceID], node)
        {
            const label& nodeID = mesh.faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N.boundaryField()[symmetricZpatchID][facei]
               *N.boundaryField()[symmetricZpatchID][facei]);

            lmN[nodeID] = iMnCn_ & lmN[nodeID];
        }
    }
}