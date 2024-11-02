// Strong boundary conditions for nodal linear momentum


// Boundary patches
const polyBoundaryMesh& bm = mesh().boundaryMesh();
const label& symmetricPatchID_ = bm.findPatchID("symmetric");
const label& symmetricXpatchID_ = bm.findPatchID("symmetricX");
const label& symmetricYpatchID_ = bm.findPatchID("symmetricY");
const label& symmetricZpatchID_ = bm.findPatchID("symmetricZ");

// Symmetric boundary patch
if (symmetricPatchID_ != -1)
{
    forAll(mesh().boundary()[symmetricPatchID_], facei)
    {
        const label& faceID =
            mesh().boundary()[symmetricPatchID_].start() + facei;

        forAll(mesh().faces()[faceID], node)
        {
            const label& nodeID = mesh().faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N_.boundaryField()[symmetricPatchID_][facei]
               *N_.boundaryField()[symmetricPatchID_][facei]);

            lmN_[nodeID] = iMnCn_ & lmN_[nodeID];
        }
    }
}

// SymmetricX boundary patch
if (symmetricXpatchID_ != -1)
{
    forAll(mesh().boundary()[symmetricXpatchID_], facei)
    {
        const label& faceID =
            mesh().boundary()[symmetricXpatchID_].start() + facei;

        forAll(mesh().faces()[faceID], node)
        {
            const label& nodeID = mesh().faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N_.boundaryField()[symmetricXpatchID_][facei]
               *N_.boundaryField()[symmetricXpatchID_][facei]);

            lmN_[nodeID] = iMnCn_ & lmN_[nodeID];
        }
    }
}

// SymmetricY boundary patch
if (symmetricYpatchID_ != -1)
{
    forAll(mesh().boundary()[symmetricYpatchID_], facei)
    {
        const label& faceID =
            mesh().boundary()[symmetricYpatchID_].start() + facei;

        forAll(mesh().faces()[faceID], node)
        {
            const label& nodeID = mesh().faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N_.boundaryField()[symmetricYpatchID_][facei]
               *N_.boundaryField()[symmetricYpatchID_][facei]);

            lmN_[nodeID] = iMnCn_ & lmN_[nodeID];
        }
    }
}

// SymmetricZ boundary patch
if (symmetricZpatchID_ != -1)
{
    forAll(mesh().boundary()[symmetricZpatchID_], facei)
    {
        const label& faceID =
            mesh().boundary()[symmetricZpatchID_].start() + facei;

        forAll(mesh().faces()[faceID], node)
        {
            const label& nodeID = mesh().faces()[faceID][node];

            const tensor& iMnCn_ =
                tensor::I
              - (N_.boundaryField()[symmetricZpatchID_][facei]
               *N_.boundaryField()[symmetricZpatchID_][facei]);

            lmN_[nodeID] = iMnCn_ & lmN_[nodeID];
        }
    }
}