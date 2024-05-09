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
#include "multiplyCoeffRectMat.H"
#include "sparseMatrixTools.H"
#include "sparseMatrixExtendedTools.H"
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
//			Info << "dualFaceI: " << dualFaceI << endl;
//			Info << "matrix ( " << ownPointID << ", " << pointID << " ) " << matrix(ownPointID, pointID) << endl;
//			Info << "coeff ( " << ownPointID << ", " << pointID << " ) " << coeff << endl;

            // Add the coefficient to the neiPointID equation coming from
            // pointID
            matrix(neiPointID, pointID) -= coeff;
//			Info << "matrix ( " << neiPointID << ", " << pointID << " ) " << matrix(neiPointID, pointID) << endl;
//			Info << "coeff ( " << neiPointID << ", " << pointID << " ) " << coeff << endl;
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
//		Info << "dualFaceI: " << dualFaceI << endl;
//		Info << "matrix ( " << ownPointID << ", " << ownPointID << " ) " << matrix(ownPointID, ownPointID) << endl;
//		Info << "coeff ( " << ownPointID << ", " << ownPointID << " ) " << edgeDirCoeff << endl;
//		Info << "matrix ( " << ownPointID << ", " << neiPointID << " ) " << matrix(ownPointID, neiPointID) << endl;
//		Info << "coeff ( " << ownPointID << ", " << neiPointID << " ) " << edgeDirCoeff << endl;

        // Insert coefficients for the neiPoint
        matrix(neiPointID, neiPointID) -= edgeDirCoeff;
        matrix(neiPointID, ownPointID) += edgeDirCoeff;
//		Info << "matrix ( " << neiPointID << ", " << neiPointID << " ) " << matrix(neiPointID, neiPointID) << endl;
//		Info << "coeff ( " << neiPointID << ", " << neiPointID << " ) " << edgeDirCoeff << endl;
//		Info << "matrix ( " << neiPointID << ", " << ownPointID << " ) " << matrix(neiPointID, ownPointID) << endl;
//		Info << "coeff ( " << neiPointID << ", " << ownPointID << " ) " << edgeDirCoeff << endl;

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
        Info<< "void Foam::vfvm::divSigma(...): end" << endl;
    }
}


void Foam::vfvm::divSigma
(
    sparseMatrixExtended& matrix,
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

        // Sensitivity term at the dual mesh face
        const RectangularMatrix<scalar>& geometricStiffness =
                geometricStiffnessField[dualFaceI];

        // Sigma at the dual mesh face
        const symmTensor sigmaf = sigmaField[dualFaceI];

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
        const tensor& dualGradDf = dualGradDField[dualFaceI];

        // Calculate F for dual faces
        const tensor& dualFf = I + dualGradDf.T();

        // Calculate invF for dual faces
        const tensor& dualInvFf = inv(dualFf);

        // Calculate J for dual faces
        const scalar& dualJf = det(dualFf);

        // Calculate deformed Sf
        const vector& curDualSfDef = (dualJf*dualInvFf.T()) & curDualSf;

        // Least squares vectors for cellID
        const vectorList& curLeastSquaresVecs = leastSquaresVecs[cellID];

        // Unit edge vector from the own point to the nei point
        vector edgeDir = points[neiPointID] - points[ownPointID];
        const scalar edgeLength = mag(edgeDir);
        edgeDir /= edgeLength;

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

            // Calculate the displacement coefficient for this point coming
            // from dualFaceI
            tensor coeff;
            multiplyCoeffExtended
                (
                    coeff,
                    curDualSfDef,
                    materialTangent,
                    geometricStiffness,
                    sigmaf,
                    lsVec
                );

            // Insert pressure coefficients
            for (int i = 0; i < 3; i++)
            {
                // Add the coefficient to the ownPointID equation
                matrix(ownPointID, pointID)(i,3) +=
                    -curDualSfDef.component(i)/curCellPoints.size();

                // Add the coefficient to the neiPointID equation coming from ownPointID
                matrix(neiPointID, pointID)(i,3) -=
                    -curDualSfDef.component(i)/curCellPoints.size();
            }

            // Insert displacement coefficients
            label cmptI = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    // Add the coefficient to the ownPointID equation coming from
                    // pointID
                    matrix(ownPointID, pointID)(i,k) += coeff.component(cmptI);

                    // Add the coefficient to the neiPointID equation coming from
                    // pointID
                    matrix(neiPointID, pointID)(i,k) -= coeff.component(cmptI);
                    cmptI++;
                }
            }
        }

        // Add compact central-differencing component in the edge direction
        // This is the gradient in the direction of the edge

        // Edge unit direction vector divided by the edge length, so that we
        // can reuse the multiplyCoeff function
        const vector eOverLength = edgeDir/edgeLength;

        // Compact edge direction displacement coefficient
        tensor edgeDirCoeff;
        multiplyCoeffExtended
        (
            edgeDirCoeff,
            curDualSfDef,
            materialTangent,
            geometricStiffness,
            sigmaf,
            eOverLength
        );
        edgeDirCoeff *= zeta;

        // Insert displacement edgeDir coefficients
        label cmptI = 0;
        for (int i = 0; i < 3; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                // Insert coefficients for the ownPoint
                matrix(ownPointID, ownPointID)(i,k) -=
                    edgeDirCoeff.component(cmptI);
                matrix(ownPointID, neiPointID)(i,k) +=
                    edgeDirCoeff.component(cmptI);

                // Insert coefficients for the neiPoint
                matrix(neiPointID, neiPointID)(i,k) -= edgeDirCoeff.component(cmptI);
                matrix(neiPointID, ownPointID)(i,k) += edgeDirCoeff.component(cmptI);

                cmptI++;
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::divSigma(...): end" << endl;
    }
}


void Foam::vfvm::divSigma
(
    sparseMatrixExtended& matrix,
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

            // Calculate the displacement coefficient for this point coming
            // from dualFaceI
            tensor coeff;
            multiplyCoeffRectMat
            (
                coeff,
                curDualSf,
                materialTangent,
                lsVec
            );

            // Insert pressure coefficients
            for (int i = 0; i < 3; i++)
            {
                // Add the coefficient to the ownPointID equation
                matrix(ownPointID, pointID)(i,3) +=
                    -curDualSf.component(i)/curCellPoints.size();

                // Add the coefficient to the neiPointID equation coming from ownPointID
                matrix(neiPointID, pointID)(i,3) -=
                    -curDualSf.component(i)/curCellPoints.size();
            }

            // Insert displacement coefficients
            label cmptI = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    // Add the coefficient to the ownPointID equation coming from
                    // pointID
                    matrix(ownPointID, pointID)(i,k) += coeff.component(cmptI);

                    // Add the coefficient to the neiPointID equation coming from
                    // pointID
                    matrix(neiPointID, pointID)(i,k) -= coeff.component(cmptI);

                    cmptI++;
                }
            }
        }

        // Add compact central-differencing component in the edge direction
        // This is the gradient in the direction of the edge

        // Edge unit direction vector divided by the edge length, so that we
        // can reuse the multiplyCoeff function
        const vector eOverLength = edgeDir/edgeLength;

        // Compact edge direction displacement coefficient
        tensor edgeDirCoeff;
        multiplyCoeffRectMat
        (
            edgeDirCoeff,
            curDualSf,
            materialTangent,
            eOverLength
        );
        edgeDirCoeff *= zeta;

        // Insert displacement edgeDir coefficients
        label cmptI = 0;
        for (int i = 0; i < 3; i++)
        {
            for (int k = 0; k < 3; k++)
            {
                // Insert coefficients for the ownPoint
                matrix(ownPointID, ownPointID)(i,k) -=
                    edgeDirCoeff.component(cmptI);
                matrix(ownPointID, neiPointID)(i,k) +=
                    edgeDirCoeff.component(cmptI);

                // Insert coefficients for the neiPoint
                matrix(neiPointID, neiPointID)(i,k) -=
                    edgeDirCoeff.component(cmptI);
                matrix(neiPointID, ownPointID)(i,k) +=
                    edgeDirCoeff.component(cmptI);

                cmptI++;
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::divSigma(...): end" << endl;
    }
}


void Foam::vfvm::laplacian
(
    sparseMatrixExtended& matrix,
    const Switch compactStencil,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const tensorField& dualGradDField,
    const scalar& diffusivity,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::laplacian(...): start" << endl;
    }

    // Take references for clarity and efficiency
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

        if (compactStencil)
        {
            // Compact central-differencing Laplacian

            // Dual face unit normal
            const scalar curDualMagSf = mag(curDualSfDef);
            const vector curDualN = curDualSfDef/curDualMagSf;

            // Delta coefficient
            const scalar deltaCoeff =
                1.0/(curDualN & (points[neiPointID] - points[ownPointID]));

            // Compact edge direction coefficient
            const scalar coeff = diffusivity*curDualMagSf*deltaCoeff;

            // Insert coefficients for the ownPoint
            matrix(ownPointID, ownPointID)(3,3) -= coeff;
            matrix(ownPointID, neiPointID)(3,3) += coeff;

            // Insert coefficients for the neiPoint
            matrix(neiPointID, neiPointID)(3,3) -= coeff;
            matrix(neiPointID, ownPointID)(3,3) += coeff;
        }
        else
        {
            // dualFaceI will contribute coefficients to the equation for each
            // primary mesh point in the dual own cell, and, if an internal
            // face, the dual neighbour cell

            // Least squares vectors for cellID
            const vectorList& curLeastSquaresVecs = leastSquaresVecs[cellID];

            forAll(curCellPoints, cpI)
            {
                // Primary point index
                const label pointID = curCellPoints[cpI];

                // Take a copy of the least squares vector from the centre of
                // cellID to pointI
                const vector lsVec = curLeastSquaresVecs[cpI];

                // Calculate the coefficient for this point coming from dualFaceI
                const scalar coeff = diffusivity*curDualSfDef & lsVec;

                // Add the coefficient to the ownPointID equation coming from
                // pointID
                matrix(ownPointID, pointID)(3,3) += coeff;

                // Add the coefficient to the neiPointID equation coming from
                // pointID
                matrix(neiPointID, pointID)(3,3) -= coeff;
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::laplacian(...): end" << endl;
    }
}


void Foam::vfvm::laplacian
(
    sparseMatrixExtended& matrix,
    const Switch compactStencil,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const scalar& diffusivity,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::laplacian(...): start" << endl;
    }

    // Take references for clarity and efficiency
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

        if (compactStencil)
        {
            // Compact central-differencing Laplacian

            // Dual face unit normal
            const scalar curDualMagSf = mag(curDualSf);
            const vector curDualN = curDualSf/curDualMagSf;

            // Delta coefficient
            const scalar deltaCoeff =
                1.0/(curDualN & (points[neiPointID] - points[ownPointID]));

            // Compact edge direction coefficient
            const scalar coeff = diffusivity*curDualMagSf*deltaCoeff;

            // Insert coefficients for the ownPoint
            matrix(ownPointID, ownPointID)(3,3) -= coeff;
            matrix(ownPointID, neiPointID)(3,3) += coeff;

            // Insert coefficients for the neiPoint
            matrix(neiPointID, neiPointID)(3,3) -= coeff;
            matrix(neiPointID, ownPointID)(3,3) += coeff;
        }
        else
        {
            // dualFaceI will contribute coefficients to the equation for each
            // primary mesh point in the dual own cell, and, if an internal
            // face, the dual neighbour cell

            // Least squares vectors for cellID
            const vectorList& curLeastSquaresVecs = leastSquaresVecs[cellID];

            forAll(curCellPoints, cpI)
            {
                // Primary point index
                const label pointID = curCellPoints[cpI];

                // Take a copy of the least squares vector from the centre of
                // cellID to pointI
                const vector lsVec = curLeastSquaresVecs[cpI];

                // Calculate the coefficient for this point coming from dualFaceI
                const scalar coeff = diffusivity*curDualSf & lsVec;

                // Add the coefficient to the ownPointID equation coming from
                // pointID
                matrix(ownPointID, pointID)(3,3) += coeff;

                // Add the coefficient to the neiPointID equation coming from
                // pointID
                matrix(neiPointID, pointID)(3,3) -= coeff;
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::laplacian(...): end" << endl;
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


void Foam::vfvm::Sp
(
    sparseMatrixExtended& matrix,
    const scalar s,
    const fvMesh& mesh,
    const labelList& dualCellToPoint,
    const scalarField& pointVolI,
    const tensorField& pBarSensitivity,
    const int debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::Sp(...): start" << endl;
    }

    // Insert pressure coefficient of pressure equation
    forAll(pointVolI, pointI)
    {
        matrix(pointI, pointI)(3,3) += s*pointVolI[pointI];
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::Sp(...): end" << endl;
    }
}


void Foam::vfvm::divU
(
    sparseMatrixExtended& matrix,
    const fvMesh& mesh,
    const labelList& dualCellToPoint,
    const scalarField& pointVolI,
    const tensorField& pBarSensitivity,
    const int debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::divU(...): start" << endl;
    }

    // Take references for clarity and efficiency
    const labelListList& pointPoints = mesh.pointPoints();
    const pointPointLeastSquaresVectors& pointPointLeastSquaresVecs =
        pointPointLeastSquaresVectors::New(mesh);
    const List<vectorList>& leastSquaresVecs =
        pointPointLeastSquaresVecs.vectors();

    forAll(pBarSensitivity, pointI)
    {
        // pBarSensitivity for pointI
        const tensor& pBarSensitivityI = pBarSensitivity[pointI];

        // Point-point neighbours
        const labelList& curPointPoints = pointPoints[pointI];

        // Least squares vectors for pointI
        const vectorList& curLeastSquaresVecs = leastSquaresVecs[pointI];

        // Accumulate contribution to the point from each point
        // neighbour
        forAll(curPointPoints, ppI)
        {
            const label pointPointID = curPointPoints[ppI];

            // Least squares vector from the average point to this neighbour
            // point
            const vector& lsVec = curLeastSquaresVecs[ppI];

            // Calculate the coefficient for this pointPoint
            const vector coeff = (pBarSensitivityI & lsVec)*pointVolI[pointI];

            // Insert the displacement coefficient of pressure equation
            for (int i = 0; i < 3; i++)
            {
                matrix(pointI, pointPointID)(3,i) += coeff.component(i);
                matrix(pointI, pointI)(3,i) -= coeff.component(i);
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::divU(...): end" << endl;
    }
}


void Foam::vfvm::grad
(
    sparseMatrixExtended& matrix,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::grad(...): start" << endl;
    }

    // Take reference for clarity and efficiency
    const labelListList& cellPoints = mesh.cellPoints();
    const labelList& dualOwn = dualMesh.owner();
    const labelList& dualNei = dualMesh.neighbour();
    const vectorField& dualSf = dualMesh.faceAreas();

    // Loop over all internal faces of the dual mesh
    forAll(dualOwn, dualFaceI)
    {
        // Primary mesh cell in which dualFaceI resides
        const label cellID = dualFaceToCell[dualFaceI];

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

        forAll(curCellPoints, cpI)
        {
            // Primary point index
            const label pointID = curCellPoints[cpI];

            // Insert pressure coefficients
            for (int i = 0; i < 3; i++)
            {
                // Add the coefficient to the ownPointID equation
                matrix(ownPointID, pointID)(i, 3) +=
                    -curDualSf.component(i)/curCellPoints.size();

                // Add the coefficient to the neiPointID equation coming from
                // ownPointID
                matrix(neiPointID, pointID)(i, 3) -=
                    -curDualSf.component(i)/curCellPoints.size();
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::grad(...): end" << endl;
    }
}


// ************************************************************************* //
