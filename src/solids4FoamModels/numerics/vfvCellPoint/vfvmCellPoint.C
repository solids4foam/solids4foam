/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "vfvmCellPoint.H"
#include "multiplyCoeff.H"
#include "sparseMatrixTools.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::vfvm::divSigma
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const Field<RectangularMatrix<scalar>>& materialTangentField,
    const boolList& fixedDofs,
    const symmTensorField& fixedDofDirections,
    const scalar fixedDofScale,
    const scalar zeta,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::divSigma(...): start" << endl;
    }

    // Take reference for clarity and efficiency
    const labelListList& cellPoints = mesh.cellPoints();
    const pointField& points = mesh.points();
    const labelList& dualOwn = dualMesh.owner();
    const labelList& dualNei = dualMesh.neighbour();
    const vectorField& dualSf = dualMesh.faceAreas();
    const cellPointLeastSquaresVectors& cellPointLeastSquaresVecs =
        cellPointLeastSquaresVectors::New(mesh);
    const List<vectorList>& leastSquaresVecs =
        cellPointLeastSquaresVecs.vectors();

    // Loop over all internal faces of the dual mesh
    forAll(dualOwn, dualFaceI)
    {
        // Primary mesh cell in which dualFaceI resides
        const label cellID = dualFaceToCell[dualFaceI];

        // Material tangent at the dual mesh face
        const RectangularMatrix<scalar>& materialTangent =
            materialTangentField[dualFaceI];

        // Points in cellID
        const labelList& curCellPoints = cellPoints[cellID];

        // Dual cell owner of dualFaceI
        const label dualOwnCellID = dualOwn[dualFaceI];

        // Dual cell neighbour of dualFaceI
        const label dualNeiCellID = dualNei[dualFaceI];

        // Primary mesh point at the centre of dualOwnCellID
        const label ownPointID = dualCellToPoint[dualOwnCellID];

        // Primary mesh point at the centre of dualNeiCellID
        const label neiPointID = dualCellToPoint[dualNeiCellID];

        // dualFaceI area vector
        const vector& curDualSf = dualSf[dualFaceI];

        // Least squares vectors for cellID
        const vectorList& curLeastSquaresVecs = leastSquaresVecs[cellID];

        // Unit edge vector from the own point to the nei point
        vector edgeDir = points[neiPointID] - points[ownPointID];
        const scalar edgeLength = mag(edgeDir);
        edgeDir /= edgeLength;

        // dualFaceI will contribute coefficients to the equation for each
        // primary mesh point in the dual own cell, and, if an internal
        // face, the dual neighbour cell

        forAll(curCellPoints, cpI)
        {
            // Primary point index
            const label pointID = curCellPoints[cpI];

            // Take a copy of the least squares vector from the centre of
            // cellID to pointI
            vector lsVec = curLeastSquaresVecs[cpI];

            // Replace the component in the primary mesh edge direction with
            // a compact central-differencing calculation
            // We remove the edge direction component by multiplying the
            // least squares vectors by (I - sqr(edgeDir))
            // Note that the compact edge direction component is added below
            lsVec = ((I - zeta*sqr(edgeDir)) & lsVec);

            // Calculate the coefficient for this point coming from dualFaceI
            tensor coeff;
            multiplyCoeff(coeff, curDualSf, materialTangent, lsVec);

            // Add the coefficient to the ownPointID equation coming from
            // pointID
            matrix(ownPointID, pointID) += coeff;

            // Add the coefficient to the neiPointID equation coming from
            // pointID
            matrix(neiPointID, pointID) -= coeff;
        }

        // Add compact central-differencing component in the edge direction
        // This is the gradient in the direction of the edge

        // Edge unit direction vector divided by the edge length, so that we
        // can reuse the multiplyCoeff function
        const vector eOverLength = edgeDir/edgeLength;

        // Compact edge direction coefficient
        tensor edgeDirCoeff;
        multiplyCoeff
        (
            edgeDirCoeff, curDualSf, materialTangent, eOverLength
        );
        edgeDirCoeff *= zeta;

        // Insert coefficients for the ownPoint
        matrix(ownPointID, ownPointID) -= edgeDirCoeff;
        matrix(ownPointID, neiPointID) += edgeDirCoeff;

        // Insert coefficients for the neiPoint
        matrix(neiPointID, neiPointID) -= edgeDirCoeff;
        matrix(neiPointID, ownPointID) += edgeDirCoeff;
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::divSigma(...): end" << endl;
    }
}


void Foam::vfvm::d2dt2
(
    ITstream& d2dt2Scheme,
    const scalar& deltaT,
    const word& pointDname,
    sparseMatrix& matrix,
    const scalarField& pointRhoI,
    const scalarField& pointVolI,
    const int debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::d2dt2(...): start" << endl;
    }

    // Read time scheme
    const word d2dt2SchemeName(d2dt2Scheme);

    // Second order identity as a 9 component tensor
    const tensor I2(I);

    // Add transient term coefficients
    if (d2dt2SchemeName == "steadyState")
    {
        // Do nothing
    }
    else if (d2dt2SchemeName == "Euler")
    {
        forAll(pointRhoI, pointI)
        {
            matrix(pointI, pointI) -=
                I2*pointVolI[pointI]*pointRhoI[pointI]/sqr(deltaT);
        }
    }
    else if (d2dt2SchemeName == "backward")
    {
        forAll(pointRhoI, pointI)
        {
            matrix(pointI, pointI) -=
                (9.0/4.0)*I2*pointVolI[pointI]*pointRhoI[pointI]/sqr(deltaT);
        }
    }
    else if (d2dt2SchemeName == "NewmarkBeta")
    {
        const scalar beta(readScalar(d2dt2Scheme));
        forAll(pointRhoI, pointI)
        {
            matrix(pointI, pointI) -=
                I2*pointVolI[pointI]*pointRhoI[pointI]/(beta*sqr(deltaT));
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::d2dt2(...): end" << endl;
    }
}

// ************************************************************************* //
