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

#include "pointPointLeastSquaresVectors.H"
#include "volFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointPointLeastSquaresVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::pointPointLeastSquaresVectors::pointPointLeastSquaresVectors
(
    const fvMesh& mesh
)
:
#ifdef OPENFOAMESIORFOUNDATION
    MeshObject<fvMesh, MoveableMeshObject, pointPointLeastSquaresVectors>(mesh),
#else
    MeshObject<fvMesh, pointPointLeastSquaresVectors>(mesh),
#endif
    vectorsPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::pointPointLeastSquaresVectors::~pointPointLeastSquaresVectors()
{
    vectorsPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointPointLeastSquaresVectors::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "pointPointLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    if (Pstream::parRun())
    {
        notImplemented
        (
            "Not yet implemented for parallel runs:"
            " actually it will work, but the results at parallel boundary "
            " will not be consistent with serial runs"
        );
    }

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = mesh_;
#else
    const fvMesh& mesh = this->mesh();
#endif

    vectorsPtr_.set(new List<vectorList>(mesh.nPoints()));
    List<vectorList>& ls = vectorsPtr_();

    // Set local references to mesh data
    const labelListList& pointPoints = mesh.pointPoints();
    const vectorField& points = mesh.points();

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh.nPoints(), symmTensor::zero);

    forAll(dd, pointI)
    {
        const vector& curPoint = points[pointI];
        const labelList& curPointPoints = pointPoints[pointI];

        forAll(curPointPoints, ppI)
        {
            const label pointPointID = curPointPoints[ppI];

            // Vector from curPoint to the point-point
            const vector d = points[pointPointID] - curPoint;

            // Weighted contribution
            // There are many ways to choose weights
            // This needs to be consistent with the loop below
            //const scalar w = 1.0;
            const scalar w = 1.0/magSqr(d);
            //const scalar w = 1.0/Foam::pow(mag(d), 3);

            dd[pointI] += w*sqr(d);
        }
    }

    // Invert least squares matrix using Householder transformations to avoid
    // badly posed cells
#ifdef OPENFOAMESIORFOUNDATION
    const symmTensorField invDd(inv(dd));
#else
    const symmTensorField invDd(hinv(dd));
#endif

    // Revisit all points and calculate the ls vectors
    forAll(ls, pointI)
    {
        const labelList& curPointPoints = pointPoints[pointI];
        const vector& curPoint = points[pointI];
        vectorList& curLs = ls[pointI];
        curLs.resize(curPointPoints.size(), vector::zero);

        forAll(curPointPoints, ppI)
        {
            const label pointPointID = curPointPoints[ppI];

            // Vector from the curPoint to the point-point
            const vector d = points[pointPointID] - curPoint;

            // Weighted contribution
            //const scalar w = 1.0;
            const scalar w = 1.0/magSqr(d);
            //const scalar w = 1.0/Foam::pow(mag(d), 3);

            curLs[ppI] = w*(invDd[pointI] & d);
        }
    }

    if (debug)
    {
        Info<< "pointPointLeastSquaresVectors::makeLeastSquaresVectors() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::List<Foam::vectorList>&
Foam::pointPointLeastSquaresVectors::vectors() const
{
    if (vectorsPtr_.empty())
    {
        makeLeastSquaresVectors();
    }

    return vectorsPtr_();
}


#ifdef OPENFOAMESIORFOUNDATION
    bool Foam::pointPointLeastSquaresVectors::movePoints()
#else
    bool Foam::pointPointLeastSquaresVectors::movePoints() const
#endif
{
    if (debug)
    {
        InfoIn("bool pointPointLeastSquaresVectors::movePoints() const")
            << "Clearing least square data" << endl;
    }

    vectorsPtr_.clear();

    return true;
}


bool Foam::pointPointLeastSquaresVectors::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn
        (
            "bool pointPointLeastSquaresVectors::updateMesh(const mapPolyMesh&) const"
        )   << "Clearing least square data" << endl;
    }

    vectorsPtr_.clear();

    return true;
}

// ************************************************************************* //
