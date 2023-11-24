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
    pointNormalContactModel

\*---------------------------------------------------------------------------*/

#include "pointNormalContactModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointNormalContactModel, 0);
defineRunTimeSelectionTable(pointNormalContactModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointNormalContactModel::pointNormalContactModel
(
    const word& name,
    const fvMesh& mesh, //const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
:
    name_(name),
    mesh_(mesh),
    masterPatchID_(masterPatchID),
    slavePatchID_(slavePatchID) //,
    // slaveContactPointGap_
    // (
    //     patch.boundaryMesh().mesh().boundaryMesh()[slavePatchID].nPoints(), 0.0
    // )
{
    if (slavePatchID == -1)
    {
        FatalErrorIn("pointNormalContactModel::pointNormalContactModel")
            << "slavePatchID cannot be -1" << abort(FatalError);
    }
}


pointNormalContactModel::pointNormalContactModel(const pointNormalContactModel& nm)
:
    name_(nm.name_),
    mesh_(nm.mesh_),
    masterPatchID_(nm.masterPatchID_),
    slavePatchID_(nm.slavePatchID_) //,
    // slaveContactPointGap_(nm.slaveContactPointGap_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
