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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "blockGlobalPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "foamTime.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        blockGlobalPolyPatch,
        0
    );

    addToRunTimeSelectionTable(polyPatch, blockGlobalPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, blockGlobalPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockGlobalPolyPatch::blockGlobalPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm)
{}


Foam::blockGlobalPolyPatch::blockGlobalPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm)
{}


Foam::blockGlobalPolyPatch::blockGlobalPolyPatch
(
    const blockGlobalPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm)
{}


//- Construct as copy, resetting the face list and boundary mesh data
Foam::blockGlobalPolyPatch::blockGlobalPolyPatch
(
    const blockGlobalPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockGlobalPolyPatch::~blockGlobalPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::blockGlobalPolyPatch::initAddressing()
{
    polyPatch::initAddressing();
}


void Foam::blockGlobalPolyPatch::calcAddressing()
{
    polyPatch::calcAddressing();
}


void Foam::blockGlobalPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}


void Foam::blockGlobalPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}


void Foam::blockGlobalPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}


void Foam::blockGlobalPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
}


void Foam::blockGlobalPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::blockGlobalPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
}


void Foam::blockGlobalPolyPatch::initOrder(const primitivePatch&) const
{}


bool Foam::blockGlobalPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return false;
}


void Foam::blockGlobalPolyPatch::syncOrder() const
{}


void Foam::blockGlobalPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
}


// ************************************************************************* //
