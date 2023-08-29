/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "cellPointLeastSquaresVectors.H"
#include "volFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellPointLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::cellPointLeastSquaresVectors::cellPointLeastSquaresVectors(const fvMesh& mesh)
:
#ifdef OPENFOAM_NOT_EXTEND
    MeshObject<fvMesh, MoveableMeshObject, cellPointLeastSquaresVectors>(mesh),
#else
    MeshObject<fvMesh, cellPointLeastSquaresVectors>(mesh),
#endif
    vectorsPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::cellPointLeastSquaresVectors::~cellPointLeastSquaresVectors()
{
    vectorsPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellPointLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "cellPointLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = mesh_;
#else
    const fvMesh& mesh = this->mesh();
#endif

    vectorsPtr_.set(new List<vectorList>(mesh.nCells()));
    List<vectorList>& ls = vectorsPtr_();

    // Set local references to mesh data
    const labelListList& cellPoints = mesh.cellPoints();
    const vectorField& points = mesh.points();

    // Calculate average position in each cell
    vectorField cellAvPoint(mesh.nCells(), vector::zero);
    forAll(cellAvPoint, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];

        forAll(curCellPoints, cpI)
        {
            const label pointID = curCellPoints[cpI];

            cellAvPoint[cellI] += points[pointID];
        }

        cellAvPoint[cellI] /= curCellPoints.size();
    }

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh.nCells(), symmTensor::zero);

    forAll(dd, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
        const vector& avPoint = cellAvPoint[cellI];

        forAll(curCellPoints, cpI)
        {
            const label pointID = curCellPoints[cpI];

            // Vector from the avPoint to the point
            const vector d = points[pointID] - avPoint;

            // Weighted contribution
            // There are many ways to choose weights
            // This needs to be consistent with the loop below
            //const scalar w = 1.0;
            const scalar w = 1.0/magSqr(d);
            //const scalar w = 1.0/Foam::pow(mag(d), 3);

            dd[cellI] += w*sqr(d);
        }
    }

    // Invert least squares matrix using Householder transformations to avoid
    // badly posed cells
#ifdef OPENFOAM_NOT_EXTEND
    const symmTensorField invDd(inv(dd));
#else
    const symmTensorField invDd(hinv(dd));
#endif

    // Revisit all cells and calculate the ls vectors
    forAll(ls, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
        const vector& avPoint = cellAvPoint[cellI];
        vectorList& curLs = ls[cellI];
        curLs.resize(curCellPoints.size(), vector::zero);

        forAll(curCellPoints, cpI)
        {
            const label pointID = curCellPoints[cpI];

            // Vector from the avPoint to the point
            const vector d = points[pointID] - avPoint;

            // Weighted contribution
            //const scalar w = 1.0;
            const scalar w = 1.0/magSqr(d);
            //const scalar w = 1.0/Foam::pow(mag(d), 3);

            curLs[cpI] = w*(invDd[cellI] & d);
        }
    }

    // Revisit all cells (again) and subtract the cell-averaged ls vectors from
    // each ls vector
    forAll(ls, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];
        vectorList& curLs = ls[cellI];

        // Calculate cell-averaged ls vector
        vector cellAvLsVec = vector::zero;
        forAll(curCellPoints, cpI)
        {
            cellAvLsVec += curLs[cpI];
        }
        cellAvLsVec /= curCellPoints.size();

        // Subtract cell-averaged ls vector from each ls vector
        forAll(curCellPoints, cpI)
        {
            curLs[cpI] -= cellAvLsVec;
        }
    }

    if (debug)
    {
        Info<< "cellPointLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::List<Foam::vectorList>&
Foam::cellPointLeastSquaresVectors::vectors() const
{
    if (vectorsPtr_.empty())
    {
        makeLeastSquaresVectors();
    }

    return vectorsPtr_();
}


#ifdef OPENFOAM_NOT_EXTEND
    bool Foam::cellPointLeastSquaresVectors::movePoints()
#else
    bool Foam::cellPointLeastSquaresVectors::movePoints() const
#endif
{
    if (debug)
    {
        InfoIn("bool cellPointLeastSquaresVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    vectorsPtr_.clear();

    return true;
}


bool Foam::cellPointLeastSquaresVectors::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn
        (
            "bool cellPointLeastSquaresVectors::updateMesh(const mapPolyMesh&) const"
        )   << "Clearing least square data" << endl;
    }

    vectorsPtr_.clear();

    return true;
}

// ************************************************************************* //
