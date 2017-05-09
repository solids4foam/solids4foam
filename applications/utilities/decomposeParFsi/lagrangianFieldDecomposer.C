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

Description
    Lagrangian field decomposer.

\*---------------------------------------------------------------------------*/

#include "lagrangianFieldDecomposer.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
lagrangianFieldDecomposer::lagrangianFieldDecomposer
(
    const polyMesh& mesh,
    const polyMesh& procMesh,
    const labelList& cellProcAddressing,
    const word& cloudName,
    const Cloud<indexedParticle>& lagrangianPositions,
    const List<SLList<indexedParticle*>*>& cellParticles
)
:
    procMesh_(procMesh),
    positions_(procMesh, cloudName, false),
    particleIndices_(lagrangianPositions.size())
{
    label pi = 0;

    forAll(cellProcAddressing, procCelli)
    {
        label celli = cellProcAddressing[procCelli];

        if (cellParticles[celli])
        {
            SLList<indexedParticle*>& particlePtrs = *cellParticles[celli];

            forAllIter(SLList<indexedParticle*>, particlePtrs, iter)
            {
                const indexedParticle& ppi = *iter();
                particleIndices_[pi++] = ppi.index();

                positions_.append
                (
                    new passiveParticle
                    (
                        positions_,
                        ppi.position(),
                        procCelli
                    )
                );
            }
        }
    }

    particleIndices_.setSize(pi);

    IOPosition<passiveParticle>(positions_).write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
