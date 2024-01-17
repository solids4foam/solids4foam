/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blockGlobalFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockGlobalFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, blockGlobalFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockGlobalFvPatch::~blockGlobalFvPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make patch weighting factors
void Foam::blockGlobalFvPatch::makeWeights(scalarField& w) const
{}


// Make patch face - neighbour cell distances
void Foam::blockGlobalFvPatch::makeDeltaCoeffs(scalarField& dc) const
{}


// Make patch face non-orthogonality correction vectors
void Foam::blockGlobalFvPatch::makeCorrVecs(vectorField& cv) const
{}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::blockGlobalFvPatch::delta() const
{
    if (patch().size())
    {
        FatalErrorIn("blockGlobalFvPatch::delta()")
            << "This patch should not have faces!" << abort(FatalError);
    }

    return tmp<vectorField>(new vectorField(0));
}

Foam::tmp<Foam::labelField> Foam::blockGlobalFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    if (patch().size())
    {
        FatalErrorIn("blockGlobalFvPatch::interfaceInternalField()")
            << "This patch should not have faces!" << abort(FatalError);
    }

    return tmp<labelField>(new labelField(0));
}


void Foam::blockGlobalFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{}


Foam::tmp<Foam::labelField> Foam::blockGlobalFvPatch::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    if (patch().size())
    {
        FatalErrorIn("blockGlobalFvPatch::transfer()")
            << "This patch should not have faces!" << abort(FatalError);
    }

    return tmp<labelField>(new labelField(0));
}


void Foam::blockGlobalFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{}


Foam::tmp<Foam::labelField> Foam::blockGlobalFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList& iF
) const
{
    if (patch().size())
    {
        FatalErrorIn("blockGlobalFvPatch::internalFieldTransfer()")
            << "This patch should not have faces!" << abort(FatalError);
    }

    return tmp<labelField>(new labelField(0));
}



// ************************************************************************* //
