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

#include "beamDisplacements.H"
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
    defineTypeNameAndDebug(beamDisplacements, 0);
    addToRunTimeSelectionTable(functionObject, beamDisplacements, dictionary);
}
}
// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::functionObjects::beamDisplacements::writeData()
{
    if (patchFound_)
    {
        const fvMesh& mesh =
            time_.lookupObject<fvMesh>(regionName_);

        if (mesh.foundObject<volVectorField>("W"))
        {
            // Beam displacement field
            const vectorField& pW =
                mesh.lookupObject<volVectorField>
                (
                    "W"
                ).boundaryField()[historyPatchID_];

            // // Value of displacement at patch
            // const vectorField& patchW =
            //     W.boundaryField()[historyPatchID_];

            const vector avgW = gAverage(pW);
            if (Pstream::master())
            {
                historyFilePtr_()
                << mesh.time().value() << " "
                << avgW.x() << " "
                << avgW.y() << " "
                << avgW.z() << " "
                << mag(avgW)
                << endl;
            }
        }
        else
        {
             Info<< this->name() + " function object constructor"
                << "W not found" << endl;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamDisplacements::beamDisplacements
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

    // word historyPatchName("notSpecified");

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
                    historyDir/"beamDisplacements_"+historyPatchName_+".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "Time" << " "
                    << "w_x" << " "
                    << "w_y" << " "
                    << "w_z" << " "
                    << "mag(w)"
                    << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamDisplacements::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::beamDisplacements::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamDisplacements::end()
{
    return true;
}


bool Foam::functionObjects::beamDisplacements::write()
{
    return true;
}


// ************************************************************************* //
