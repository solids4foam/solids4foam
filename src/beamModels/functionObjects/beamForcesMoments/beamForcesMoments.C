/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamForcesMoments.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(beamForcesMoments, 0);
    addToRunTimeSelectionTable(functionObject, beamForcesMoments, dictionary);
}
}
// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::functionObjects::beamForcesMoments::writeData()
{
    if (patchFound_)
    {
        const fvMesh& mesh =
            time_.lookupObject<fvMesh>(regionName_);

        if (mesh.foundObject<surfaceVectorField>("Q"))
        {
            // Beam force field
            const vectorField& pQ =
                mesh.lookupObject<surfaceVectorField>
                (
                    "Q"
                ).boundaryField()[historyPatchID_];

            // Beam moment field
            const vectorField& pM =
                mesh.lookupObject<surfaceVectorField>
                (
                    "M"
                ).boundaryField()[historyPatchID_];

            const vector avgQ = gAverage(pQ);
            const vector avgM = gAverage(pM);

            if (Pstream::master())
            {
                historyFilePtr_()
                << mesh.time().value() << " "
                << avgQ.x() << " "
                << avgQ.y() << " "
                << avgQ.z() << " "
                << avgM.x() << " "
                << avgM.y() << " "
                << avgM.z()
                << endl;
            }
        }
        else
        {
            Info<< this->name() + " function object constructor"
                << "Q not found" << endl;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamForcesMoments::beamForcesMoments
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    name_(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(),
    historyPatchName_("notSpecified"),
    patchFound_(false),
    historyPatchID_(-1)
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch") >> historyPatchName_;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "solidDisplacements: historyPatch not specified" << endl;
    }

    // Lookup the beam mesh
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName_);

    if (historyPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName_ << " not found"
            << endl;
    }
    else
    {
        patchFound_ = true;
    }

    // Create history file if not already created
    if (!historyFilePtr_.valid() && patchFound_)
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir/"beamForcesMoments_"+historyPatchName_+".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "Time" << " "
                    << "Q_x" << " "
                    << "Q_y" << " "
                    << "Q_z" << " "
                    << "M_x" << " "
                    << "M_y" << " "
                    << "M_z"
                    << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamForcesMoments::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::beamForcesMoments::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamForcesMoments::end()
{
    return true;
}


bool Foam::functionObjects::beamForcesMoments::write()
{
    return true;
}


// ************************************************************************* //
