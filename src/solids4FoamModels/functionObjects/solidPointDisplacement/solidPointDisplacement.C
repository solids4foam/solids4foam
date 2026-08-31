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
#include "lookupSolidModel.H"
#ifndef FOAMEXTEND
    #include "symmTensor3rdOrder.H"
#endif
#ifdef OPENFOAM_NOT_EXTEND
    #include "volPointInterpolation.H"
#else
    #include "newLeastSquaresVolPointInterpolation.H"
#endif
#ifdef OPENFOAM_ORG
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

#ifndef FOAMEXTEND
void Foam::solidPointDisplacement::extrapolatedPointDisplacement()
{
    // Lookup the solid mesh
    const fvMesh& mesh = time_.lookupObject<fvMesh>(region_);

    // Lookup solid model
    const solidModel& solMod = lookupSolidModel(mesh);

    // Initialise resulting displacement
    vector pointDValue = vector::zero;

    if (mesh.foundObject<volVectorField>("D"))
    {
        if (cellID_ > -1)
        {
            const volVectorField& D = mesh.lookupObject<volVectorField>("D");
            pointDValue = D[cellID_];

            const vector distance = point_ - mesh.C()[cellID_];

            if (mesh.foundObject<volTensorField>("grad(D)"))
            {
                const volTensorField& gradD =
                    mesh.lookupObject<volTensorField>("grad(D)");

                pointDValue += distance & gradD[cellID_];
            }
            else
            {
                InfoInFunction
                    << "volTensorField grad(D) not found" << endl;

                // Without gradient available we will write closest cell centre
                // displacement
                closestPointDisplacement();
            }

            // Extrapolate from cell using higher order gradient
            if
            (
                solMod.highOrderResidual()
             && solMod.displacementMLS().polynomialOrder() >= 2
            )
            {
                const bool thirdOrder =
                    solMod.displacementMLS().polynomialOrder() >= 3;

                // The higher derivatives of all three components share a
                // single halo exchange and a single traversal of the stencil.
                // Only the cell values are needed here, so no fields are built
                FixedList<symmTensorField, 3> secondGradD;
                FixedList<List<symmTensor3rdOrder>, 3> thirdGradD;

                solMod.displacementMLS().cellDerivatives
                (
                    D,
                    nullptr,
                    &secondGradD,
                    thirdOrder ? &thirdGradD : nullptr
                );

                pointDValue.x() +=
                    0.5*(distance & secondGradD[0][cellID_] & distance);
                pointDValue.y() +=
                    0.5*(distance & secondGradD[1][cellID_] & distance);
                pointDValue.z() +=
                    0.5*(distance & secondGradD[2][cellID_] & distance);

                if (thirdOrder)
                {
                    pointDValue.x() +=
                        (1.0/6.0)*cubicForm(thirdGradD[0][cellID_], distance);
                    pointDValue.y() +=
                        (1.0/6.0)*cubicForm(thirdGradD[1][cellID_], distance);
                    pointDValue.z() +=
                        (1.0/6.0)*cubicForm(thirdGradD[2][cellID_], distance);
                }
            }
        }
    }
    else
    {
        InfoInFunction
            << "volVectorField D not found" << endl;
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
#endif


void Foam::solidPointDisplacement::closestPointDisplacement()
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
}


bool Foam::solidPointDisplacement::writeData()
{
#ifndef FOAMEXTEND
    if (extrapolate_)
    {
        extrapolatedPointDisplacement();
    }
    else
    {
        closestPointDisplacement();
    }
#else
    closestPointDisplacement();
#endif

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
    point_(dict.lookup("point")),
    pointID_(-1),
    cellID_(-1),
#ifndef FOAMEXTEND
    extrapolate_(dict.lookupOrDefault<bool>("extrapolate", false)),
#endif
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
        scalar minPointDist = GREAT;

        forAll(mesh.points(), pI)
        {
            scalar dist = mag(mesh.points()[pI] - point_);

            if (dist < minPointDist)
            {
                minPointDist = dist;
                pointID_ = pI;
            }
        }

        // Find the closest cell
        scalar minCellDist = GREAT;

        forAll(mesh.C(), cellI)
        {
            const scalar dist = mag(mesh.C()[cellI] - point_);

            if (dist < minCellDist)
            {
                minCellDist = dist;
                cellID_ = cellI;
            }
        }

        // Find global closest point
        const scalar globalPointMinDist =
            returnReduce(minPointDist, minOp<scalar>());
        int pointProcNo = -1;
        if (mag(globalPointMinDist - minPointDist) < SMALL)
        {
            pointProcNo = Pstream::myProcNo();
        }
        else
        {
            pointID_ = -1;
        }

        // Find global closest cell
        const scalar globalCellMinDist =
            returnReduce(minCellDist, minOp<scalar>());
        int cellProcNo = -1;
        if (mag(globalCellMinDist - minCellDist) < SMALL)
        {
            cellProcNo = Pstream::myProcNo();
        }
        else
        {
            cellID_ = -1;
        }

        // More than one processor can have the point so we will take the proc
        // with the highest processor number
        const int globalMaxPointProc = returnReduce(pointProcNo, maxOp<int>());
        if (globalMaxPointProc != pointProcNo)
        {
            pointID_ = -1;
        }

        if (pointID_ > -1)
        {
            Pout<< "    distance from specified point is " << minPointDist
                << endl;
        }

        const int globalMaxCellProc = returnReduce(cellProcNo, maxOp<int>());
        if (globalMaxCellProc != cellProcNo)
        {
            cellID_ = -1;
        }

        if
	(
	    returnReduce(pointID_, maxOp<int>()) == -1
	 || returnReduce(cellID_, maxOp<int>()) == -1
	)
        {
            FatalErrorIn("solidPointDisplacement::solidPointDisplacement")
                << "Something went wrong: no proc found a closest point/cell!"
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


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::solidPointDisplacement::write()
{
    return false;
}
#endif

// ************************************************************************* //
