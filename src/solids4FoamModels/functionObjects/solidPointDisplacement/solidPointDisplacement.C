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

\*----------------------------------------------------------------------------*/

#include "solidPointDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "volPointInterpolation.H"
#else
    #include "newLeastSquaresVolPointInterpolation.H"
#endif
#ifdef OPENFOAMFOUNDATION
    #include "OSspecific.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPointDisplacement, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPointDisplacement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPointDisplacement::writeData()
{
    // Lookup the solid mesh
    const fvMesh& mesh = time_.lookupObject<fvMesh>(region_);

    if (mesh.foundObject<pointVectorField>("pointD"))
    {
        // Lookup the point displacement field
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        vector pointDValue = vector::zero;
        if (pointID_ > -1)
        {
            pointDValue = pointD[pointID_];
        }
        reduce(pointDValue, sumOp<vector>());

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value()
                << " " << pointDValue.x()
                << " " << pointDValue.y()
                << " " << pointDValue.z()
                << " " << mag(pointDValue)
                << endl;
        }
    }
    else
    {
        InfoIn(this->name() + " function object constructor")
            << "pointVectorField pointD not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPointDisplacement::solidPointDisplacement
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    region_(dict.lookupOrDefault<word>("region", "UNDEFINED")),
    pointID_(-1),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;

    // Set region if it is undefined
    if (region_ == "UNDEFINED")
    {
        if (time_.foundObject<fvMesh>("solid"))
        {
            region_ = "solid";
        }
        else
        {
            region_ = fvMesh::defaultRegion;
        }
    }
    Info<< "    region = " << region_ << endl;

    // Lookup the point
    const vector point(dict.lookup("point"));

    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // Find the closest point
        scalar minDist = GREAT;

        forAll(mesh.points(), pI)
        {
            scalar dist = mag(mesh.points()[pI] - point);

            if (dist < minDist)
            {
                minDist = dist;
                pointID_ = pI;
            }
        }

        // Find global closest point
        const scalar globalMinDist = returnReduce(minDist, minOp<scalar>());
        int procNo = -1;
        if (mag(globalMinDist - minDist) < SMALL)
        {
            procNo = Pstream::myProcNo();
        }
        else
        {
            pointID_ = -1;
        }

        // More than one processor can have the point so we will take the proc
        // with the highest processor number
        const int globalMaxProc = returnReduce(procNo, maxOp<int>());
        if (mag(globalMaxProc - procNo) > SMALL)
        {
            pointID_ = -1;
        }

        if (pointID_ > -1)
        {
            Pout<< "    distance from specified point is " << minDist
                << endl;
        }

        if (returnReduce(pointID_, maxOp<int>()) == -1)
        {
            FatalErrorIn("solidPointDisplacement::solidPointDisplacement")
                << "Something went wrong: no proc found a point!"
                << abort(FatalError);
        }

        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            const word startTimeName =
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
                    historyDir/"solidPointDisplacement_" + name + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time Dx Dy Dz magD" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPointDisplacement::start()
{
    return false;
}


#if FOAMEXTEND
bool Foam::solidPointDisplacement::execute(const bool forceWrite)
#else
bool Foam::solidPointDisplacement::execute()
#endif
{
    return writeData();
}


bool Foam::solidPointDisplacement::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::solidPointDisplacement::write()
{
    return false;
}
#endif

// ************************************************************************* //
