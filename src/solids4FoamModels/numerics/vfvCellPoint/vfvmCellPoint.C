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

#ifdef USE_PETSC

#include "vfvmCellPoint.H"
#include "multiplyCoeff.H"
#include "cellPointLeastSquaresVectors.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::vfvm::divSigma
(
    Mat jac,
    const pointVectorField& pointD,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const label nScalarEqns,
    const labelList& localToGlobalPointMap,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const List<mat66>& materialTangentField,
    const scalar zeta,
    const bool flipSign
)
{
    const scalar sign = flipSign ? -1.0 : 1.0;
    const label colOffset = 0;
    const label rowOffset = 0;

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

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise block coeff
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

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
        const mat66& materialTangent = materialTangentField[dualFaceI];

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

        // Calculate the global owner and neighbour point indices
        const label globalOwnPointID = localToGlobalPointMap[ownPointID];
        const label globalNeiPointID = localToGlobalPointMap[neiPointID];

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

            // Calculate the global owner point index
            const label globalPointID = localToGlobalPointMap[pointID];

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
            multiplyCoeff(coeff, curDualSf, materialTangent, lsVec);
            coeff *= sign;

            // Construct the block coeff
            for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
            {
                for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                {
                    values
                    [
                        (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                    ] = coeff[cmptI*3 + cmptJ];
                }
            }

            // Add the coefficient to the ownPointID equation coming from
            // pointID
            // matrix(ownPointID, pointID) += coeff;
            MatSetValuesBlocked
            (
                jac, 1, &globalOwnPointID, 1, &globalPointID,
                values.cdata(),
                ADD_VALUES
            );

            // Add the coefficient to the neiPointID equation coming from
            // pointID
            // matrix(neiPointID, pointID) -= coeff;
            forAll(values, cmptI)
            {
                values[cmptI] = -values[cmptI];
            }
            MatSetValuesBlocked
            (
                jac, 1, &globalNeiPointID, 1, &globalPointID,
                values.cdata(),
                ADD_VALUES
            );
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
        edgeDirCoeff *= sign;

        // Construct the block coeff
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
            {
                values
                [
                    (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                ] = edgeDirCoeff[cmptI*3 + cmptJ];
            }
        }

        // Insert coefficients for the ownPoint-neiPoint
        // matrix(ownPointID, neiPointID) += edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalOwnPointID, 1, &globalNeiPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Insert coefficients for the neiPoint-ownPointID
        // matrix(neiPointID, ownPointID) += edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalNeiPointID, 1, &globalOwnPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Flip the coefficient signs
        forAll(values, cmptI)
        {
            values[cmptI] = -values[cmptI];
        }

        // Insert coefficients for the ownPoint-ownPoint
        // matrix(ownPointID, ownPointID) -= edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalOwnPointID, 1, &globalOwnPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Insert coefficients for the neiPoint-neiPoint
        // matrix(neiPointID, neiPointID) -= edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalNeiPointID, 1, &globalNeiPointID,
            values.cdata(),
            ADD_VALUES
        );
    }
}


void Foam::vfvm::divSigma
(
    Mat jac,
    const pointVectorField& pointD,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const label nScalarEqns,
    const labelList& localToGlobalPointMap,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const List<mat66>& materialTangentField,
    const List<mat39>& geometricStiffnessField,
    const symmTensorField& sigmaField,
    const tensorField& dualGradDField,
    const scalar zeta,
    const bool flipSign
)
{
    const scalar sign = flipSign ? -1.0 : 1.0;
    const label colOffset = 0;
    const label rowOffset = 0;

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

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise block coeff
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

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
        const mat66& materialTangent = materialTangentField[dualFaceI];

        // Sensitivity term at the dual mesh face
        const mat39& geometricStiffness = geometricStiffnessField[dualFaceI];

        // Info<< "dualFaceI = " << dualFaceI << " ";
        // for (int i = 0; i < 3; ++i)
        // {
        //     for (int j = 0; j < 9; ++j)
        //     {
        //         if (mag(geometricStiffness(i, j)) > SMALL)
        //         {
        //             Info<< " " << geometricStiffness(i, j);
        //         }
        //     }
        // }
        // Info<< endl;

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

        // Calculate the global owner and neighbour point indices
        const label globalOwnPointID = localToGlobalPointMap[ownPointID];
        const label globalNeiPointID = localToGlobalPointMap[neiPointID];

        // dualFaceI area vector
        const vector& curDualSf = dualSf[dualFaceI];

        // gradD at the dual mesh face
        const tensor& dualGradD = dualGradDField[dualFaceI];

        // Calculate F for dual faces
        const tensor dualF(I + dualGradD.T());

        // Calculate invF for dual faces
        const tensor dualInvF(inv(dualF));

        // Calculate J for dual faces
        const scalar dualJ(det(dualF));

        // Calculate deformed Sf
        const vector curDualSfDef((dualJ*dualInvF.T()) & curDualSf);

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

            // Calculate the global owner point index
            const label globalPointID = localToGlobalPointMap[pointID];

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
            multiplyCoeff
            (
                coeff,
                curDualSfDef,
                materialTangent,
                geometricStiffness,
                sigma,
                lsVec
            );
            coeff *= sign;

            // Construct the block coeff
            for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
            {
                for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                {
                    values
                    [
                        (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                    ] = coeff[cmptI*3 + cmptJ];
                }
            }

            // Add the coefficient to the ownPointID equation coming from
            // pointID
            //matrix(ownPointID, pointID) += coeff;
            MatSetValuesBlocked
            (
                jac, 1, &globalOwnPointID, 1, &globalPointID,
                values.cdata(),
                ADD_VALUES
            );

            // Add the coefficient to the neiPointID equation coming from
            // pointID
            //matrix(neiPointID, pointID) -= coeff;
            forAll(values, cmptI)
            {
                values[cmptI] = -values[cmptI];
            }
            MatSetValuesBlocked
            (
                jac, 1, &globalNeiPointID, 1, &globalPointID,
                values.cdata(),
                ADD_VALUES
            );
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
            edgeDirCoeff,
            curDualSfDef,
            materialTangent,
            geometricStiffness,
            sigma,
            eOverLength
        );
        edgeDirCoeff *= zeta;
        edgeDirCoeff *= sign;

        // Construct the block coeff
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
            {
                values
                [
                    (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                ] = edgeDirCoeff[cmptI*3 + cmptJ];
            }
        }

        // Insert coefficients for the ownPoint-neiPoint
        // matrix(ownPointID, neiPointID) += edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalOwnPointID, 1, &globalNeiPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Insert coefficients for the neiPoint-ownPointID
        // matrix(neiPointID, ownPointID) += edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalNeiPointID, 1, &globalOwnPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Flip the coefficient signs
        forAll(values, cmptI)
        {
            values[cmptI] = -values[cmptI];
        }

        // Insert coefficients for the ownPoint-ownPoint
        // matrix(ownPointID, ownPointID) -= edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalOwnPointID, 1, &globalOwnPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Insert coefficients for the neiPoint-neiPoint
        // matrix(neiPointID, neiPointID) -= edgeDirCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalNeiPointID, 1, &globalNeiPointID,
            values.cdata(),
            ADD_VALUES
        );
    }
}


void Foam::vfvm::d2dt2
(
    Mat jac,
    const pointVectorField& pointD,
    const scalarField& pointRhoI,
    const scalarField& pointVolI,
    ITstream& d2dt2Scheme,
    const label nScalarEqns,
    const labelList& localToGlobalPointMap,
    const bool flipSign
)
{
    const scalar sign = flipSign ? -1.0 : 1.0;
    const scalar deltaT = pointD.mesh().time().deltaTValue();
    const label colOffset = 0;
    const label rowOffset = 0;

    // Read time scheme
    const word d2dt2SchemeName(d2dt2Scheme);

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Calculate the scalar coefficient field
    scalarField coeffs(pointD.size(), 0.0);

    // Add transient term coefficients
    if (d2dt2SchemeName == "steadyState")
    {
        // Do nothing
        return;
    }
    else if (d2dt2SchemeName == "Euler")
    {
        forAll(pointRhoI, pointI)
        {
            coeffs[pointI] =
                sign*pointVolI[pointI]*pointRhoI[pointI]/sqr(deltaT);
        }
    }
    else if (d2dt2SchemeName == "backward")
    {
        forAll(pointRhoI, pointI)
        {
            coeffs[pointI] =
                sign*(9.0/4.0)*pointVolI[pointI]*pointRhoI[pointI]/sqr(deltaT);
        }
    }
    else if (d2dt2SchemeName == "NewmarkBeta")
    {
        const scalar beta(readScalar(d2dt2Scheme));
        forAll(pointRhoI, pointI)
        {
            coeffs[pointI] =
                sign*pointVolI[pointI]*pointRhoI[pointI]/(beta*sqr(deltaT));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown d2dt2SchemeName scheme = " << d2dt2SchemeName
            << exit(FatalError);
    }

    // Initialise block coeff
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

    // Insert the coeffs into the PETSc matrix
    forAll(pointRhoI, pointI)
    {
        // Construct the block coeff
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
            {
                if (cmptI == cmptJ)
                {
                    values
                    [
                        (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                    ] = coeffs[pointI];
                }
            }
        }

        // Determine the global block row/column
        const label globalBlockRowI = localToGlobalPointMap[pointI];

        // Insert the block coefficient
        MatSetValuesBlocked
        (
            jac, 1, &globalBlockRowI, 1, &globalBlockRowI,
            values.cdata(),
            ADD_VALUES
        );
    }
}


void Foam::vfvm::laplacian
(
    Mat jac,
    const Switch compactStencil,
    const scalar zeta, // compact stencil coefficient
    const meshDual& dualMesh,
    const label nScalarEqns,
    const labelList& localToGlobalPointMap,
    const scalarField& diffusivity, // diffusivity in the primary cells
    const bool flipSign
)
{
    const scalar sign = flipSign ? -1.0 : 1.0;
    const label colOffset = 0;
    const label rowOffset = 0;

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

    // Get the blockSize
    label blockSize;
    MatGetBlockSize(jac, &blockSize);

    // Initialise block coeff
    const label nCoeffCmpts = blockSize*blockSize;
    List<PetscScalar> values(nCoeffCmpts, 0.0);

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

        // Calculate the global owner and neighbour point indices
        const label globalOwnPointID = localToGlobalPointMap[ownPointID];
        const label globalNeiPointID = localToGlobalPointMap[neiPointID];

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

            // Calculate the global owner point index
            const label globalPointID = localToGlobalPointMap[pointID];

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
            tensor coeff((curDualSf & lsVec)*I2);
            coeff *= sign;

            // Construct the block coeff
            for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
            {
                for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
                {
                    values
                    [
                        (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                    ] = coeff[cmptI*3 + cmptJ];
                }
            }

            // Add the coefficient to the ownPointID equation coming from
            // pointID
            // matrix(ownPointID, pointID) += coeff;
            MatSetValuesBlocked
            (
                jac, 1, &globalOwnPointID, 1, &globalPointID,
                values.cdata(),
                ADD_VALUES
            );

            // Add the coefficient to the neiPointID equation coming from
            // pointID
            // matrix(neiPointID, pointID) -= coeff;
            forAll(values, cmptI)
            {
                values[cmptI] = -values[cmptI];
            }
            MatSetValuesBlocked
            (
                jac, 1, &globalNeiPointID, 1, &globalPointID,
                values.cdata(),
                ADD_VALUES
            );
        }

        // Add compact central-differencing component

        // Delta coefficient
        // const scalar deltaCoeff =
        //     1.0/(curDualN & (points[neiPointID] - points[ownPointID]));

        // Compact edge direction coefficient
        tensor compactCoeff((curDualSf & eOverLength)*I2*zeta);
        //const tensor compactCoeff(curDualMagSf*deltaCoeff*I2*zeta);
        compactCoeff *= sign;

        // Construct the block coeff
        for (label cmptI = 0; cmptI < nScalarEqns; ++cmptI)
        {
            for (label cmptJ = 0; cmptJ < nScalarEqns; ++cmptJ)
            {
                values
                [
                    (cmptI + rowOffset)*blockSize + cmptJ + colOffset
                ] = compactCoeff[cmptI*3 + cmptJ];
            }
        }

        // Insert coefficients for the ownPoint-neiPoint
        // matrix(ownPointID, neiPointID) += compactCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalOwnPointID, 1, &globalNeiPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Insert coefficients for the neiPoint-ownPointID
        // matrix(neiPointID, ownPointID) += compactCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalNeiPointID, 1, &globalOwnPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Flip the coefficient signs
        forAll(values, cmptI)
        {
            values[cmptI] = -values[cmptI];
        }

        // Insert coefficients for the ownPoint-ownPoint
        // matrix(ownPointID, ownPointID) -= compactCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalOwnPointID, 1, &globalOwnPointID,
            values.cdata(),
            ADD_VALUES
        );

        // Insert coefficients for the neiPoint-neiPoint
        // matrix(neiPointID, neiPointID) -= compactCoeff;
        MatSetValuesBlocked
        (
            jac, 1, &globalNeiPointID, 1, &globalNeiPointID,
            values.cdata(),
            ADD_VALUES
        );
    }
}


// Foam::tmp<Foam::sparseMatrix> Foam::vfvm::div
// (
//     //const scalar rho,
//     const pointVectorField& pointV,
//     const meshDual& dualMesh,
//     const bool debug
// )
// {
//     if (debug)
//     {
//         Info<< "tmp<sparseScalarMatrix> vfvm::div(...): start" << endl;
//     }

//     // Prepare the result field
//     tmp<sparseMatrix> tmatrix(new sparseMatrix(20*pointV.size()));
//     sparseMatrix& matrix = tmatrix.ref();

//     // Take references for clarity and efficiency
//     const fvMesh& mesh = dualMesh.mesh();
//     const labelListList& cellPoints = mesh.cellPoints();
//     const labelList& dualOwn = dualMesh.owner();
//     const labelList& dualNei = dualMesh.neighbour();
//     const vectorField& dualSf = dualMesh.faceAreas();
//     const vectorField& pointVI = pointV;
//     const labelList& dualFaceToCell = dualMesh.dualMeshMap().dualFaceToCell();
//     const labelList& dualCellToPoint = dualMesh.dualMeshMap().dualCellToPoint();
//     //const scalarField& dualRhoI = dualRho;
//     // const cellPointLeastSquaresVectors& cellPointLeastSquaresVecs =
//     //     cellPointLeastSquaresVectors::New(mesh);
//     // const List<vectorList>& leastSquaresVecs =
//     //     cellPointLeastSquaresVecs.vectors();

//     // Loop over all internal faces of the dual mesh
//     forAll(dualOwn, dualFaceI)
//     {
//         // Primary mesh cell in which dualFaceI resides
//         const label cellID = dualFaceToCell[dualFaceI];

//         // Points in cellID
//         const labelList& curCellPoints = cellPoints[cellID];

//         // Dual cell owner of dualFaceI
//         const label dualOwnCellID = dualOwn[dualFaceI];

//         // Dual cell neighbour of dualFaceI
//         const label dualNeiCellID = dualNei[dualFaceI];

//         // Primary mesh point at the centre of dualOwnCellID
//         const label ownPointID = dualCellToPoint[dualOwnCellID];

//         // Primary mesh point at the centre of dualNeiCellID
//         const label neiPointID = dualCellToPoint[dualNeiCellID];

//         // dualFaceI area vector
//         const vector& curDualSf = dualSf[dualFaceI];

//         // dualFaceI density
//         // const scalar curDualRho = dualRhoI[dualFaceI];

//         forAll(curCellPoints, cpI)
//         {
//             // Primary point index
//             const label pID = curCellPoints[cpI];

//             // Weight for pID
//             const scalar wp = 1.0/curCellPoints.size();

//             forAll(curCellPoints, cqI)
//             {
//                 // Primary point index
//                 const label qID = curCellPoints[cqI];

//                 // Weight for qID
//                 const scalar wq = 1.0/curCellPoints.size();

//                 // Calculate the coefficient for point pID coming from dualFaceI
//                 // const tensor ddpCoeff(curDualRho*wp*wq*curDualSf*pointVI[qID]);
//                 const tensor ddpCoeff(wp*wq*curDualSf*pointVI[qID]);

//                 // Calculate the coefficient for point qID coming from dualFaceI
//                 const tensor ddqCoeff
//                 (
//                     // curDualRho*wp*wq*(curDualSf & pointVI[pID])*I
//                     wp*wq*(curDualSf & pointVI[pID])*I
//                 );

//                 // Add the coefficient to the ownPointID equation coming from
//                 // pID
//                 matrix(ownPointID, pID) += ddpCoeff;

//                 // Add the coefficient to the ownPointID equation coming from
//                 // qID
//                 matrix(ownPointID, qID) += ddqCoeff;

//                 // Subtract the coeffs from the neighbour point equation as the
//                 // dual face normal is flipped
//                 matrix(neiPointID, pID) -= ddpCoeff;
//                 matrix(neiPointID, qID) -= ddqCoeff;
//             }
//         }
//     }

//     if (debug)
//     {
//         Info<< "tmp<sparseScalarMatrix> vfvm::div(...): end" << endl;
//     }
//     return tmatrix;
// }

#endif // ifdef USE_PETSC

// ************************************************************************* //
