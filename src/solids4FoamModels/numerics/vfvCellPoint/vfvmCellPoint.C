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

#include "vfvmCellPoint.H"
#include "multiplyCoeff.H"
#include "multiplyCoeffExtended.H"
#include "sparseMatrixTools.H"
#include "cellPointLeastSquaresVectors.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::vfvm::divSigma
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const Field<scalarSquareMatrix>& materialTangentField,
    const scalar zeta,
    const bool debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::divSigma(...): start" << endl;
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

    // Check the material tangents are the correct shape
    forAll(materialTangentField, faceI)
    {
        if (materialTangentField[faceI].m() != 6)
        {
            FatalErrorIn("void Foam::vfvm::divSigma(...)")
                << "The materialTangent for face " << faceI << " has "
                << materialTangentField[faceI].m() << " rows "
                << "but it should have 6!" << abort(FatalError);
        }
        else if (materialTangentField[faceI].n() != 6)
        {
            FatalErrorIn("void Foam::vfvm::divSigma(...)")
                << "The materialTangent for face " << faceI << " has "
                << materialTangentField[faceI].m() << " rows "
                << "but it should have 6!" << abort(FatalError);
        }
    }

    // Loop over all internal faces of the dual mesh
    forAll(dualOwn, dualFaceI)
    {
        // Primary mesh cell in which dualFaceI resides
        const label cellID = dualFaceToCell[dualFaceI];

        // Material tangent at the dual mesh face
        const scalarSquareMatrix& materialTangent =
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
            // multiplyCoeffRectMat(coeff, curDualSf, materialTangent, lsVec);
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
        // multiplyCoeffRectMat
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


void Foam::vfvm::divSigma
(
    sparseMatrix& matrix,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const Field<scalarSquareMatrix>& materialTangentField,
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

    // Check the material tangents are the correct shape
    forAll(materialTangentField, faceI)
    {
        if (materialTangentField[faceI].m() != 6)
        {
            FatalErrorIn("void Foam::vfvm::divSigma(...)")
                << "The materialTangent for face " << faceI << " has "
                << materialTangentField[faceI].m() << " rows "
                << "but it should have 6!" << abort(FatalError);
        }
        else if (materialTangentField[faceI].n() != 6)
        {
            FatalErrorIn("void Foam::vfvm::divSigma(...)")
                << "The materialTangent for face " << faceI << " has "
                << materialTangentField[faceI].m() << " rows "
                << "but it should have 6!" << abort(FatalError);
        }
    }

    // Loop over all internal faces of the dual mesh
    forAll(dualOwn, dualFaceI)
    {
        // Primary mesh cell in which dualFaceI resides
        const label cellID = dualFaceToCell[dualFaceI];

        // Material tangent at the dual mesh face
        const scalarSquareMatrix& materialTangent =
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

    if (debug)
    {
        Info<< "void Foam::vfvm::divSigma(...): end" << endl;
    }
}


void Foam::vfvm::laplacian
(
    sparseScalarMatrix& matrix,
    const Switch compactStencil,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const scalarField& diffusivity,
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

        // Diffusivity in the primary cells
        const scalar diffCellID = diffusivity[cellID];

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
            const scalar coeff = diffCellID*curDualMagSf*deltaCoeff;

            // Insert coefficients for the ownPoint
            matrix(ownPointID, ownPointID) -= coeff;
            matrix(ownPointID, neiPointID) += coeff;

            // Insert coefficients for the neiPoint
            matrix(neiPointID, neiPointID) -= coeff;
            matrix(neiPointID, ownPointID) += coeff;
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
                const scalar coeff = diffCellID*curDualSf & lsVec;

                // Add the coefficient to the ownPointID equation coming from
                // pointID
                matrix(ownPointID, pointID) += coeff;

                // Add the coefficient to the neiPointID equation coming from
                // pointID
                matrix(neiPointID, pointID) -= coeff;
            }
        }
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::laplacian(...): end" << endl;
    }
}


Foam::tmp<Foam::sparseMatrix> Foam::vfvm::laplacian
(
    const pointVectorField& pointP,
    const meshDual& dualMesh,
    const scalar zeta,
    const bool debug
)
{
    if (debug)
    {
        Info<< "tmp<sparseMatrix> Foam::vfvm::laplacian(...): start" << endl;
    }

    // Prepare the result field
    tmp<sparseMatrix> tmatrix(new sparseMatrix(20*pointP.size()));
    sparseMatrix& matrix = tmatrix.ref();

    // Take references for clarity and efficiency
    const fvMesh& mesh = dualMesh.mesh();
    const labelListList& cellPoints = mesh.cellPoints();
    const pointField& points = mesh.points();
    const labelList& dualOwn = dualMesh.owner();
    const labelList& dualNei = dualMesh.neighbour();
    const vectorField& dualSf = dualMesh.faceAreas();
    const cellPointLeastSquaresVectors& cellPointLeastSquaresVecs =
        cellPointLeastSquaresVectors::New(mesh);
    const List<vectorList>& leastSquaresVecs =
        cellPointLeastSquaresVecs.vectors();
    const labelList& dualFaceToCell = dualMesh.dualMeshMap().dualFaceToCell();
    const labelList& dualCellToPoint = dualMesh.dualMeshMap().dualCellToPoint();

    // Second order identity as a 9 component tensor
    const tensor I2(I);

    // Loop over all internal faces of the dual mesh
    forAll(dualOwn, dualFaceI)
    {
        // Primary mesh cell in which dualFaceI resides
        const label cellID = dualFaceToCell[dualFaceI];

        // Diffusivity in the primary cells
        // const scalar diffCellID = diffusivity[cellID];

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

        // Edge unit direction vector divided by the edge length, so that we
        // can reuse the multiplyCoeff function
        const vector eOverLength = edgeDir/edgeLength;

        // Dual face unit normal
        // const scalar curDualMagSf = mag(curDualSf);
        // const vector curDualN = curDualSf/curDualMagSf;

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
            //lsVec = ((I - zeta*sqr(curDualN)) & lsVec);

            // Calculate the coefficient for this point coming from dualFaceI
            const tensor coeff((curDualSf & lsVec)*I2);

            // Add the coefficient to the ownPointID equation coming from
            // pointID
            matrix(ownPointID, pointID) += coeff;

            // Add the coefficient to the neiPointID equation coming from
            // pointID
            matrix(neiPointID, pointID) -= coeff;
        }

        // Add compact central-differencing component

        // Delta coefficient
        // const scalar deltaCoeff =
        //     1.0/(curDualN & (points[neiPointID] - points[ownPointID]));

        // Compact edge direction coefficient
        const tensor compactCoeff((curDualSf & eOverLength)*I2*zeta);
        //const tensor compactCoeff(curDualMagSf*deltaCoeff*I2*zeta);

        // Insert coefficients for the ownPoint
        matrix(ownPointID, ownPointID) -= compactCoeff;
        matrix(ownPointID, neiPointID) += compactCoeff;

        // Insert coefficients for the neiPoint
        matrix(neiPointID, neiPointID) -= compactCoeff;
        matrix(neiPointID, ownPointID) += compactCoeff;
    }

    return tmatrix;

    if (debug)
    {
        Info<< "tmp<sparseMatrix> Foam::vfvm::laplacian(...): end" << endl;
    }
}


Foam::tmp<Foam::sparseMatrix> Foam::vfvm::div
(
    //const scalar rho,
    const pointVectorField& pointV,
    const meshDual& dualMesh,
    const bool debug
)
{
    if (debug)
    {
        Info<< "tmp<sparseScalarMatrix> vfvm::div(...): start" << endl;
    }

    // Prepare the result field
    tmp<sparseMatrix> tmatrix(new sparseMatrix(20*pointV.size()));
    sparseMatrix& matrix = tmatrix.ref();

    // Take references for clarity and efficiency
    const fvMesh& mesh = dualMesh.mesh();
    const labelListList& cellPoints = mesh.cellPoints();
    const labelList& dualOwn = dualMesh.owner();
    const labelList& dualNei = dualMesh.neighbour();
    const vectorField& dualSf = dualMesh.faceAreas();
    const vectorField& pointVI = pointV;
    const labelList& dualFaceToCell = dualMesh.dualMeshMap().dualFaceToCell();
    const labelList& dualCellToPoint = dualMesh.dualMeshMap().dualCellToPoint();
    //const scalarField& dualRhoI = dualRho;
    // const cellPointLeastSquaresVectors& cellPointLeastSquaresVecs =
    //     cellPointLeastSquaresVectors::New(mesh);
    // const List<vectorList>& leastSquaresVecs =
    //     cellPointLeastSquaresVecs.vectors();

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

        // dualFaceI density
        // const scalar curDualRho = dualRhoI[dualFaceI];

        forAll(curCellPoints, cpI)
        {
            // Primary point index
            const label pID = curCellPoints[cpI];

            // Weight for pID
            const scalar wp = 1.0/curCellPoints.size();

            forAll(curCellPoints, cqI)
            {
                // Primary point index
                const label qID = curCellPoints[cqI];

                // Weight for qID
                const scalar wq = 1.0/curCellPoints.size();

                // Calculate the coefficient for point pID coming from dualFaceI
                // const tensor ddpCoeff(curDualRho*wp*wq*curDualSf*pointVI[qID]);
                const tensor ddpCoeff(wp*wq*curDualSf*pointVI[qID]);

                // Calculate the coefficient for point qID coming from dualFaceI
                const tensor ddqCoeff
                (
                    // curDualRho*wp*wq*(curDualSf & pointVI[pID])*I
                    wp*wq*(curDualSf & pointVI[pID])*I
                );

                // Add the coefficient to the ownPointID equation coming from
                // pID
                matrix(ownPointID, pID) += ddpCoeff;

                // Add the coefficient to the ownPointID equation coming from
                // qID
                matrix(ownPointID, qID) += ddqCoeff;

                // Subtract the coeffs from the neighbour point equation as the
                // dual face normal is flipped
                matrix(neiPointID, pID) -= ddpCoeff;
                matrix(neiPointID, qID) -= ddqCoeff;
            }
        }
    }

    if (debug)
    {
        Info<< "tmp<sparseScalarMatrix> vfvm::div(...): end" << endl;
    }
    return tmatrix;
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

    // TODO: why do we assume minus? For consistency, this should be a positive
    // contribution!

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


Foam::tmp<Foam::sparseMatrix> Foam::vfvm::ddt
(
    ITstream& ddtScheme,
    const pointVectorField& pointP,
    const int debug
)
{
    if (debug)
    {
        Info<< "void Foam::vfvm::ddt(...): start" << endl;
    }

    // Prepare the result field
    tmp<sparseMatrix> tmatrix(new sparseMatrix(pointP.size()));
    sparseMatrix& matrix = tmatrix.ref();

    // Read time scheme
    const word ddtSchemeName(ddtScheme);

    // Take a reference to the internal field for efficiency
    const vectorField& pointPI = pointP;

    // Second order identity as a 9 component tensor
    const tensor I2(I);

    // Time-step: assumed uniform
    const scalar deltaT = pointP.time().deltaTValue();

    // Add transient term coefficients
    if (ddtSchemeName == "steadyState")
    {
        // Do nothing
    }
    else if (ddtSchemeName == "Euler")
    {
        forAll(pointPI, pointI)
        {
            matrix(pointI, pointI) += I2/deltaT;
        }
    }
    else if (ddtSchemeName == "backward")
    {
        forAll(pointPI, pointI)
        {
            matrix(pointI, pointI) += (3.0/2.0)*I2/deltaT;
        }
    }
    else
    {
        FatalError
            << "Unknown ddtScheme: " << ddtScheme << abort(FatalError);
    }

    if (debug)
    {
        Info<< "void Foam::vfvm::ddt(...): end" << endl;
    }

    return tmatrix;
}

// ************************************************************************* //
