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

\*---------------------------------------------------------------------------*/

#include "vfvmCellPointExtended.H"
#include "multiplyCoeffExtended.H"
#include "sparseMatrixTools.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::vfvm::divSigmaExtended
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const Field<RectangularMatrix<scalar>>& materialTangentField,
	const Field<RectangularMatrix<scalar>>& geometricStiffnessField,
	const symmTensorField& sigmaField,
	const tensorField& dualGradDField,
    const boolList& fixedDofs,
    const symmTensorField& fixedDofDirections,
    const scalar fixedDofScale,
    const scalar zeta,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::divSigmaExtended(...): start" << endl;
    }

//    Info << "sigmaField = " << sigmaField << endl;
//    Info << sigmaField.size() << endl;
//    Info << materialTangentField.size() << endl;
//    Info << geometricStiffnessField.size() << endl;
//    Info << dualGradDField.size() << endl;

//    for (int i = 0; i < geometricStiffnessField.size(); i++)
//    {
//        Info << "(" << i << ", 0) : " << geometricStiffnessField[i] << endl;
//
//    }
//    Info << endl;


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

        // Sensitivity term at the dual mesh face
        const RectangularMatrix<scalar>& geometricStiffness =
        	geometricStiffnessField[dualFaceI];

       	// Sigma at the dual mesh face
       	const symmTensor sigma = sigmaField[dualFaceI];

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

        // gradD at the dual mesh face
        const tensor& dualGradD = dualGradDField[dualFaceI];

		// Calculate F for dual faces
		const tensor& dualF = I + dualGradD.T();

		// Calculate invF for dual faces
		const tensor& dualInvF = inv(dualF);

		// Calculate J for dual faces
		const scalar& dualJ = det(dualF);

		// Calculate deformed Sf
		const vector& curDualSfDef = (dualJ*dualInvF.T()) & curDualSf;

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
            // cellID to pointID
            vector lsVec = curLeastSquaresVecs[cpI];

            // Replace the component in the primary mesh edge direction with
            // a compact central-differencing calculation
            // We remove the edge direction component by multiplying the
            // least squares vectors by (I - sqr(edgeDir))
            // Note that the compact edge direction component is added below
            lsVec = ((I - zeta*sqr(edgeDir)) & lsVec);

            // Calculate the coefficient for this point coming from dualFaceI
            tensor coeff;
            multiplyCoeffExtended
            (
            	coeff,
            	curDualSfDef,
            	materialTangent,
            	geometricStiffness,
            	sigma,
            	lsVec
            );

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
        multiplyCoeffExtended
        (
        	edgeDirCoeff,
        	curDualSfDef,
        	materialTangent,
        	geometricStiffness,
        	sigma,
        	eOverLength
        );
        edgeDirCoeff *= zeta;

        // Insert coefficients for the ownPoint
        matrix(ownPointID, ownPointID) -= edgeDirCoeff;
        matrix(ownPointID, neiPointID) += edgeDirCoeff;

        // Insert coefficients for the neiPoint
        matrix(neiPointID, neiPointID) -= edgeDirCoeff;
        matrix(neiPointID, ownPointID) += edgeDirCoeff;

    }

    // Enforce fixed degrees of freedom
    // This should probably be in its own function but it is ok here for now
    // Note: traction boundary conditions do not contribute to matrix
    // coefficients and instead contribute to just the source
    // sparseMatrixTools::addFixedDofToMatrix
    // (
    //     matrix, fixedDofs, fixedDofDirections, fixedDofScale
    // );

    if (debug)
    {
        Info<< "void Foam::vfvm::divSigmaExtended(...): end" << endl;
    }
}


void Foam::vfvm::d2dt2Extended
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
        Info<< "void Foam::vfvm::d2dt2Extended(...): start" << endl;
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
        Info<< "void Foam::vfvm::d2dt2Extended(...): end" << endl;
    }
}

// ************************************************************************* //
