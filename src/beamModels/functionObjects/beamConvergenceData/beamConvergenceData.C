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

#include "beamConvergenceData.H"
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
    defineTypeNameAndDebug(beamConvergenceData, 0);
    addToRunTimeSelectionTable(functionObject, beamConvergenceData, dictionary);
}
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

bool Foam::functionObjects::beamConvergenceData::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    if (Pstream::master())
    {
        historyFilePtr_()
        << time_.time().value() << " "
        << beam.iOuterCorr()
        << endl;
    }

    return true;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamConvergenceData::beamConvergenceData
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
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Create history file if not already created
    if (!historyFilePtr_.valid())
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
            OStringStream FileName;
            FileName() << "beamConvergenceData.dat";

            historyFilePtr_.reset
            (
                new OFstream(historyDir/word(FileName.str()))
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                << "Time" << " "
                << "nCorretors"
                << endl;
            }
        }
    }

    // read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamConvergenceData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::functionObjects::beamConvergenceData::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamConvergenceData::end()
{
    return true;
}


bool Foam::functionObjects::beamConvergenceData::write()
{
    return true;
}


// ************************************************************************* //
