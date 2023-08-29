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

InClass
    frictionContactModel

\*---------------------------------------------------------------------------*/

#include "frictionContactModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(frictionContactModel, 0);
defineRunTimeSelectionTable(frictionContactModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frictionContactModel::frictionContactModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID
)
:
    name_(name),
    patch_(patch),
    masterPatchID_(masterPatchID),
    slavePatchID_(slavePatchID),
    stickSlipFaces_(patch.boundaryMesh()[slavePatchID].size(), 0.0)
{}


frictionContactModel::frictionContactModel(const frictionContactModel& fm)
:
    name_(fm.name_),
    patch_(fm.patch_),
    masterPatchID_(fm.masterPatchID_),
    slavePatchID_(fm.slavePatchID_),
    stickSlipFaces_(fm.stickSlipFaces_)
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //


void frictionContactModel::autoMap(const fvPatchFieldMapper& m)
{
#ifdef OPENFOAMFOUNDATION
    m(stickSlipFaces_, stickSlipFaces_);
#else
    stickSlipFaces_.autoMap(m);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
