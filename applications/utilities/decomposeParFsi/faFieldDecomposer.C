/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "faFieldDecomposer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const label sizeBeforeMapping,
    const unallocLabelList& addressingSlice,
    const label addressingOffset
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    directAddressing_(addressingSlice)
{
    forAll (directAddressing_, i)
    {
        // Subtract one to align addressing.
        // directAddressing_[i] -= addressingOffset + 1;
        // ZT, 12/Nov/2010
        directAddressing_[i] -= addressingOffset;
    }
}


faFieldDecomposer::processorAreaPatchFieldDecomposer::
processorAreaPatchFieldDecomposer
(
    const faMesh& mesh,
    const unallocLabelList& addressingSlice
)
:
    sizeBeforeMapping_(mesh.nFaces()),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    const scalarField& weights = mesh.weights().internalField();
    const labelList& own = mesh.edgeOwner();
    const labelList& neighb = mesh.edgeNeighbour();

    forAll (addressing_, i)
    {
        // Subtract one to align addressing.
        label ai = addressingSlice[i];
//         label ai = mag(addressingSlice[i]) - 1;

        if (ai < neighb.size())
        {
            // This is a regular edge. it has been an internal edge
            // of the original mesh and now it has become a edge
            // on the parallel boundary
            addressing_[i].setSize(2);
            weights_[i].setSize(2);

            addressing_[i][0] = own[ai];
            addressing_[i][1] = neighb[ai];

            weights_[i][0] = weights[ai];
            weights_[i][1] = 1.0 - weights[ai];
        }
        else
        {
            // This is a edge that used to be on a cyclic boundary
            // but has now become a parallel patch edge. I cannot
            // do the interpolation properly (I would need to look
            // up the different (edge) list of data), so I will
            // just grab the value from the owner face
            //
            addressing_[i].setSize(1);
            weights_[i].setSize(1);

            addressing_[i][0] = own[ai];

            weights_[i][0] = 1.0;
        }
    }
}


faFieldDecomposer::processorEdgePatchFieldDecomposer::
processorEdgePatchFieldDecomposer
(
    label sizeBeforeMapping,
    const unallocLabelList& addressingSlice
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    forAll (addressing_, i)
    {
        addressing_[i].setSize(1);
        weights_[i].setSize(1);

        addressing_[i][0] = mag(addressingSlice[i]) - 1;
        weights_[i][0] = sign(addressingSlice[i]);
    }
}


faFieldDecomposer::faFieldDecomposer
(
    const faMesh& completeMesh,
    const faMesh& procMesh,
    const labelList& edgeAddressing,
    const labelList& faceAddressing,
    const labelList& boundaryAddressing
)
:
    completeMesh_(completeMesh),
    procMesh_(procMesh),
    edgeAddressing_(edgeAddressing),
    faceAddressing_(faceAddressing),
    boundaryAddressing_(boundaryAddressing),
    patchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<patchFieldDecomposer*>(NULL)
    ),
    processorAreaPatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorAreaPatchFieldDecomposer*>(NULL)
    ),
    processorEdgePatchFieldDecomposerPtrs_
    (
        procMesh_.boundary().size(),
        static_cast<processorEdgePatchFieldDecomposer*>(NULL)
    )
{
    forAll (boundaryAddressing_, patchi)
    {
        if (boundaryAddressing_[patchi] >= 0)
        {
            patchFieldDecomposerPtrs_[patchi] = new patchFieldDecomposer
            (
                completeMesh_.boundary()[boundaryAddressing_[patchi]].size(),
                procMesh_.boundary()[patchi].patchSlice(edgeAddressing_),
//                 completeMesh_.boundaryMesh()
                completeMesh_.boundary()
                [
                    boundaryAddressing_[patchi]
                ].start()
            );
        }
        else
        {
            processorAreaPatchFieldDecomposerPtrs_[patchi] =
                new processorAreaPatchFieldDecomposer
                (
                    completeMesh_,
                    procMesh_.boundary()[patchi].patchSlice(edgeAddressing_)
                );

            processorEdgePatchFieldDecomposerPtrs_[patchi] =
                new processorEdgePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    static_cast<const unallocLabelList&>
                    (
                        procMesh_.boundary()[patchi].patchSlice
                        (
                            edgeAddressing_
                        )
                    )
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faFieldDecomposer::~faFieldDecomposer()
{
    forAll (patchFieldDecomposerPtrs_, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            delete patchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll (processorAreaPatchFieldDecomposerPtrs_, patchi)
    {
        if (processorAreaPatchFieldDecomposerPtrs_[patchi])
        {
            delete processorAreaPatchFieldDecomposerPtrs_[patchi];
        }
    }

    forAll (processorEdgePatchFieldDecomposerPtrs_, patchi)
    {
        if (processorEdgePatchFieldDecomposerPtrs_[patchi])
        {
            delete processorEdgePatchFieldDecomposerPtrs_[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
