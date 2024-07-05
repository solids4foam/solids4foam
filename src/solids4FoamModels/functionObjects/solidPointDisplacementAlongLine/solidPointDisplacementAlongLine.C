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

#include "solidPointDisplacementAlongLine.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "volPointInterpolation.H"
#else
    #include "newLeastSquaresVolPointInterpolation.H"
#endif
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPointDisplacementAlongLine, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPointDisplacementAlongLine,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPointDisplacementAlongLine::writeData()
{

    // Lookup the solid mesh
    const fvMesh& mesh = time_.lookupObject<fvMesh>(region_);

    if (mesh.foundObject<pointVectorField>("pointD"))
    {
        // Lookup the point displacement field
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        //Obtain pointD for all values on the specified line
        forAll(pointID_, pI)
        {
            const vector pointDValue = pointD[pointID_[pI]];

            historyFilePtr_()
                << pointID_[pI]
                << " " << pointCoord_[pI]
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


void Foam::solidPointDisplacementAlongLine::sortByComp
(
    DynamicList<vector>& pointCoord,
    DynamicList<label>& pointID,
    const label cmptI
)
{
    for (scalar i = 0; i < pointCoord.size(); i++)
    {
        for (scalar k = i + 1; k < pointCoord.size(); k++)
        {
            // I think there should be no mag here?
            if
            (
                mag(pointCoord[i].component(cmptI))
              > mag(pointCoord[k].component(cmptI))
            )
            {
                vector tempCoord(pointCoord[i]);
                pointCoord[i] = pointCoord[k];
                pointCoord[k] = tempCoord;

                label tempID(pointID[i]);
                pointID[i] = pointID[k];
                pointID[k] = tempID;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPointDisplacementAlongLine::solidPointDisplacementAlongLine
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
    pointID_(0),
    pointCoord_(0),
    historyFilePtr_()
{
    if (Pstream::parRun())
    {
        notImplemented
        (
            "This function object is currently only implemented for serial runs"
        );
    }

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

    // Lookup the start and end point of the line
    const vector pointA(dict.lookup("startPoint"));
    const vector pointB(dict.lookup("endPoint"));

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

    // Set capacity of point lists to be a fraction of the total number of points
    pointID_.setCapacity(0.001*mesh.nPoints());
    pointCoord_.setCapacity(0.001*mesh.nPoints());

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // Find the closest point on the line
        const scalar minDist(dict.lookupOrDefault<scalar>("minDist", 1e-6));
        Info<< "    minDist: " << minDist << endl;

        // Define vector between points A and B
        const vector line_vector = pointB - pointA;

        forAll(mesh.points(), pI)
        {
            // Define vector between point A and mesh point P
            const vector point_vector = mesh.points()[pI] - pointA;

            // Create bounding box for the line from A to B
            boundBox bb(pointA, pointB);

            // Inflate the box in case it has zero volume
            bb.inflate(0.01);

            // Check whether point is within the region defined by the segment
            if (bb.contains(mesh.points()[pI]))
            {
                // Calculate coordinates of projection point on line
                const vector proj_pt
                (
                    (
                        (point_vector & line_vector)/mag(line_vector)
                    )*(line_vector/mag(line_vector))
                  + pointA
                );

                // Calculate distance between mesh point and projection point
                const scalar dist = mag(mesh.points()[pI] - proj_pt);

                // Check if mesh point is on the line
                if (dist < minDist)
                {
                    pointID_.append(pI);
                    pointCoord_.append(mesh.points()[pI]);
                }
            }
        }

        // Sort point coordinates by x, y or z-coordinates
        if
        (
            mag(pointCoord_[0].component(vector::Y))
         == mag(pointCoord_[1].component(vector::Y))
         && mag(pointCoord_[0].component(vector::Z))
         == mag(pointCoord_[1].component(vector::Z))
        )
        {
            // Sort point coordinates by x-coordinates
            sortByComp(pointCoord_, pointID_, vector::X);
        }
        else if
        (
            mag(pointCoord_[0].component(vector::X))
         == mag(pointCoord_[1].component(vector::X))
         && mag(pointCoord_[0].component(vector::Z))
         == mag(pointCoord_[1].component(vector::Z)))
        {
            // Sort point coordinates by y-coordinates
            sortByComp(pointCoord_, pointID_, vector::Y);
        }
        else if
        (
            mag(pointCoord_[0].component(vector::X))
         == mag(pointCoord_[1].component(vector::X))
         && mag(pointCoord_[0].component(vector::Y))
         == mag(pointCoord_[1].component(vector::Y))
        )
        {
            // Sort point coordinates by z-coordinates
            sortByComp(pointCoord_, pointID_, vector::Z);
        }
        else
        {
            // Sort point coordinates by x-coordinates
            sortByComp(pointCoord_, pointID_, vector::X);
            Info<< "SortByCompX since two or three coordinates are different"
                << endl;
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
                    historyDir/"solidPointDisplacementAlongLine_" + name + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# PointID PointCoord Dx Dy Dz magD" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPointDisplacementAlongLine::start()
{
    return false;
}

#if FOAMEXTEND
bool Foam::solidPointDisplacementAlongLine::execute(const bool forceWrite)
#else
bool Foam::solidPointDisplacementAlongLine::execute()
#endif
{
    return writeData();
}


bool Foam::solidPointDisplacementAlongLine::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidPointDisplacementAlongLine::write()
{
    return false;
}
#endif

// ************************************************************************* //
