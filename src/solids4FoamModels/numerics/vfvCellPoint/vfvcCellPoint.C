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

#include "vfvcCellPoint.H"
#include "cellPointLeastSquaresVectors.H"
#include "pointPointLeastSquaresVectors.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace vfvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volTensorField> grad
(
    const pointVectorField& pointD,
    const fvMesh& mesh
)
{
    // Prepare the result field
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "grad("+ pointD.name() +")",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "zero", pointD.dimensions()/dimLength, tensor::zero
            ),
            "zeroGradient"
        )
    );
#ifdef OPENFOAMESIORFOUNDATION
    volTensorField& result = tresult.ref();
#else
    volTensorField& result = tresult();
#endif

    // Take references for clarity and efficiency

    tensorField& resultI = result;
    const labelListList& cellPoints = mesh.cellPoints();
    const vectorField& pointDI = pointD.internalField();
    const cellPointLeastSquaresVectors& cellPointLeastSquaresVecs =
        cellPointLeastSquaresVectors::New(mesh);
    const List<vectorList>& leastSquaresVecs =
        cellPointLeastSquaresVecs.vectors();

    // Calculate the gradient for each cell
    forAll(resultI, cellI)
    {
        // Points in the current cell
        const labelList& curCellPoints = cellPoints[cellI];

        // Least squares vectors for cellI
        const vectorList& curLeastSquaresVecs = leastSquaresVecs[cellI];

        // Accumulate contribution to the cell gradient from each point
        tensor& cellGrad = resultI[cellI];
        forAll(curCellPoints, cpI)
        {
            // Least squares vector from the centre of cellID to pointI
            const vector& lsVec = curLeastSquaresVecs[cpI];

            // Primary point index
            const label pointID = curCellPoints[cpI];

            // Add least squares contribution to the cell gradient
            cellGrad += lsVec*pointDI[pointID];
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


tmp<pointTensorField> pGrad
(
    const pointVectorField& pointD,
    const fvMesh& mesh
)
{
    // Get a reference to the pointMesh
    const pointMesh& pMesh = pointD.mesh();

    // Prepare the result field
    tmp<pointTensorField> tresult
    (
        new pointTensorField
        (
            IOobject
            (
                "pGrad("+ pointD.name() +")",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedTensor
            (
                "zero", pointD.dimensions()/dimLength, tensor::zero
            ),
            "calculated"
        )
    );
#ifdef OPENFOAMESIORFOUNDATION
    pointTensorField& result = tresult.ref();
#else
    pointTensorField& result = tresult();
#endif

    // Take references for clarity and efficiency

    tensorField& resultI = result;
    const labelListList& pointPoints = mesh.pointPoints();
    const vectorField& pointDI = pointD.internalField();
    const pointPointLeastSquaresVectors& pointPointLeastSquaresVecs =
        pointPointLeastSquaresVectors::New(mesh);
    const List<vectorList>& leastSquaresVecs =
        pointPointLeastSquaresVecs.vectors();

    // Calculate the gradient for each point
    forAll(resultI, pointI)
    {
        // Point-point neighbours
        const labelList& curPointPoints = pointPoints[pointI];

        // Least squares vectors for pointI
        const vectorList& curLeastSquaresVecs = leastSquaresVecs[pointI];

        // Accumulate contribution to the point gradient from each point
        // neighbour
        tensor& pointGrad = resultI[pointI];
        forAll(curPointPoints, ppI)
        {
            // Least squares vector from the average point to this neighbour
            // point
            const vector& lsVec = curLeastSquaresVecs[ppI];

            // Point index
            const label pointPointID = curPointPoints[ppI];

            // Delta pointD
            const vector deltaPointD = pointDI[pointPointID] - pointDI[pointI];

            // Add least squares contribution to the gradient
            pointGrad += lsVec*deltaPointD;
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


tmp<surfaceTensorField> fGrad
(
    const pointVectorField& pointD,
    const fvMesh& mesh,
    const fvMesh& dualMesh,
    const labelList& dualFaceToCell,
    const labelList& dualCellToPoint,
    const scalar zeta,
    const bool debug
)
{
    if (debug)
    {
        Info<< "surfaceTensorField fGrad(...): start" << endl;
    }

    // Prepare the result field
    tmp<surfaceTensorField> tresult
    (
        new surfaceTensorField
        (
            IOobject
            (
                "fGrad("+ pointD.name() +")",
                dualMesh.time().timeName(),
                dualMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dualMesh,
            dimensionedTensor
            (
                "zero", pointD.dimensions()/dimLength, tensor::zero
            )
        )
    );
#ifdef OPENFOAMESIORFOUNDATION
    surfaceTensorField& result = tresult.ref();
#else
    surfaceTensorField& result = tresult();
#endif

    // Take references for clarity and efficiency
    tensorField& resultI = result;
    const vectorField& pointDI = pointD.internalField();
    const pointField& points = mesh.points();
    const labelList& dualOwn = dualMesh.faceOwner();
    const labelList& dualNei = dualMesh.faceNeighbour();

    // Approach
    // Step 1: Calculate constant gradient in each primary mesh cell
    // Step 2: Set dual face gradient to primary mesh constant cell gradient and
    //         replace the component in the edge direction

    // Calculate constant gradient in each primary mesh cell
    const volTensorField gradD(vfvc::grad(pointD, mesh));
    const tensorField& gradDI = gradD.internalField();

    // Set dual face gradient to primary mesh constant cell gradient and
    // replace the component in the edge direction
    // We only replace the edge component for internal dual faces

    // For all faces - internal and boundary
    forAll(dualOwn, dualFaceI)
    {
        // Only calculate the gradient for internal faces
        if (dualMesh.isInternalFace(dualFaceI))
        {
            // Primary mesh cell in which dualFaceI resides
            const label cellID = dualFaceToCell[dualFaceI];

            // Dual cell owner of dualFaceI
            const label dualOwnCellID = dualOwn[dualFaceI];

            // Dual cell neighbour of dualFaceI
            const label dualNeiCellID = dualNei[dualFaceI];

            // Primary mesh point at the centre of dualOwnCellID
            const label ownPointID = dualCellToPoint[dualOwnCellID];

            // Primary mesh point at the centre of dualNeiCellID
            const label neiPointID = dualCellToPoint[dualNeiCellID];

            // Unit edge vector from the own point to the nei point
            vector edgeDir = points[neiPointID] - points[ownPointID];
            const scalar edgeLength = mag(edgeDir);
            edgeDir /= edgeLength;

            // Calculate the gradient component in the edge direction using
            // central-differencing and use the primary mesh cell value for the
            // tangential directions
            resultI[dualFaceI] =
              //   edgeDir*(pointDI[neiPointID] - pointDI[ownPointID])/edgeLength
              // + ((I - sqr(edgeDir)) & gradDI[cellID]);
                zeta*edgeDir
               *(
                   pointDI[neiPointID] - pointDI[ownPointID]
               )/edgeLength
              + ((I - zeta*sqr(edgeDir)) & gradDI[cellID]);
        }
        else // boundary face
        {
            // Dual patch which this dual face resides on
            const label dualPatchID =
                dualMesh.boundaryMesh().whichPatch(dualFaceI);

            if (dualMesh.boundaryMesh()[dualPatchID].type() != "empty")
            {
                // Find local face index
                const label localDualFaceID =
                    dualFaceI - dualMesh.boundaryMesh()[dualPatchID].start();

                // Primary mesh cell in which dualFaceI resides
                const label cellID = dualFaceToCell[dualFaceI];

                // Use the gradient in the adjacent primary cell-centre
                // This will result in inconsistent values at processor patches
                // Is this an issue?
                result.boundaryFieldRef()[dualPatchID][localDualFaceID] =
                    gradDI[cellID];
            }
        }
    }

    if (debug)
    {
        Info<< "surfaceTensorField fGrad(...): end" << endl;
    }

    return tresult;
}


tmp<vectorField> d2dt2
(
    ITstream& d2dt2Scheme,
    const pointVectorField& pointD, // displacement
    const pointVectorField& pointU, // velocity
    const pointVectorField& pointA, // acceleration
    const scalarField& pointRho,    // density
    const scalarField& pointVol,    // volumes
    const int debug // debug switch
)
{
    // Take a reference to the internal field
    const vectorField& pointDI = pointD.internalField();

    // Create result field
    tmp<vectorField> tresult(new vectorField(pointDI.size(), vector::zero));
#ifdef OPENFOAMESIORFOUNDATION
    vectorField& result = tresult.ref();
#else
    vectorField& result = tresult();
#endif

    // Read the time-scheme
    const word d2dt2SchemeName(d2dt2Scheme);

    // Time-step: assumed uniform
    const dimensionedScalar deltaT = pointD.time().deltaT();

    if (d2dt2SchemeName == "steadyState")
    {
        // Do nothing
    }
    else if (d2dt2SchemeName == "Euler")
    {
        result =
        (
            pointD - 2.0*pointD.oldTime() + pointD.oldTime().oldTime()
        )().internalField()*pointVol*pointRho/Foam::pow(deltaT.value(), 2.0);
    }
    else if (d2dt2SchemeName == "backward")
    {
        result =
        (
            1.5*
            (
                1.5*pointD
              - 2.0*pointD.oldTime()
              + 0.5*pointD.oldTime().oldTime()
            )/deltaT
          - 2.0*pointU.oldTime()
          + 0.5*pointU.oldTime().oldTime()
        )().internalField()*pointVol*pointRho/deltaT.value();
    }
    else if (d2dt2SchemeName == "NewmarkBeta")
    {
        const scalar beta(readScalar(d2dt2Scheme));

        const pointField pointDbar
        (
            pointD.oldTime()
          + deltaT*pointU.oldTime()
          + 0.5*sqr(deltaT)*(1.0 - 2.0*beta)*pointA.oldTime()
        );

        result =
            (pointDI - pointDbar)*pointVol*pointRho/(beta*sqr(deltaT.value()));
    }
    else
    {
        FatalErrorIn("tmp<vectorField> d2dt2(...)")
            << "Not implemented for d2dt2Scheme = " << d2dt2SchemeName << nl
            << "Available d2dt2Schemes are: " << nl
            << "    steadyState" << nl
            << "    Euler" << nl
            << "    backward" << nl
            << "    NewmarkBeta" << nl
            << abort(FatalError);
    }

    return tresult;
}


tmp<vectorField> ddt
(
    ITstream& ddtScheme,
    ITstream& d2dt2Scheme,
    const pointVectorField& pointP
)
{
    // Take a reference to the internal field
    const vectorField& pointPI = pointP.internalField();

    // Create result field
    tmp<vectorField> tresult(new vectorField(pointPI.size(), vector::zero));
#ifdef OPENFOAMESIORFOUNDATION
    vectorField& result = tresult.ref();
#else
    vectorField& result = tresult();
#endif

    // Read ddt time-scheme
    const word ddtSchemeName(ddtScheme);

    // Check that ddt and d2dt2 schemes are consistent for pointD and pointU
    word d2dt2SchemeName("none");
    if (pointP.name() == "pointD" || pointP.name() == "pointU")
    {
        d2dt2SchemeName = word(d2dt2Scheme);

        if (ddtSchemeName != d2dt2SchemeName)
        {
            FatalErrorIn("tmp<vectorField> ddt(...)")
                << "The ddtScheme and d2dt2Scheme for " << pointP.name()
                << " are not consistant"
                << abort(FatalError);
        }
    }

    // Time-step: assumed uniform
    const dimensionedScalar deltaT = pointP.time().deltaT();

    if (ddtSchemeName == "steadyState")
    {
        // Do nothing
    }
    else if (ddtSchemeName == "Euler")
    {
        result = (pointP - pointP.oldTime())/deltaT;
    }
    else if (ddtSchemeName == "backward")
    {
        result =
        (
            1.5*pointP - 2.0*pointP.oldTime() + 0.5*pointP.oldTime().oldTime()
        )/deltaT;
    }
    else if (ddtSchemeName == "NewmarkBeta")
    {
        // Read the beta and gamma parameters
        const scalar beta(readScalar(d2dt2Scheme));
        const scalar gamma(readScalar(d2dt2Scheme));

        if (pointP.name() == "pointU")
        {
            // Lookup the point displacement field
            const pointVectorField& pointD =
                pointP.mesh().db().lookupObject<pointVectorField>("pointD");

            // Lookup the point velocity field
            const pointVectorField& pointU = pointP;

            // Lookup the point acceleration field
            const pointVectorField& pointA =
                pointP.mesh().db().lookupObject<pointVectorField>("pointA");

            const pointVectorField pointDbar
            (
                pointD.oldTime()
              + deltaT*pointU.oldTime()
              + 0.5*sqr(deltaT)*(1.0 - 2.0*beta)*pointA.oldTime()
            );

            // Acceleration
            result = (pointD - pointDbar)/(beta*sqr(deltaT));
        }
        else if (pointP.name() == "pointD")
        {
            // Lookup the point velocity field
            const pointVectorField& pointU =
                pointP.mesh().db().lookupObject<pointVectorField>("pointU");

            // Lookup the point acceleration field
            const pointVectorField& pointA =
                pointP.mesh().db().lookupObject<pointVectorField>("pointA");

            const pointVectorField pointUbar
            (
                pointU.oldTime() + (1.0 - gamma)*deltaT*pointA.oldTime()
            );

            // Velocity
            result = pointUbar + gamma*deltaT*pointA;
        }
        else
        {
            FatalErrorIn("tmp<vectorField> ddt(...)")
                << "NewmarkBeta only implemented for pointD and pointU"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("tmp<vectorField> ddt(...)")
            << "Not implemented for ddtScheme = " << ddtSchemeName << nl
            << "Available ddtSchemes are: " << nl
            << "    steadyState" << nl
            << "    Euler" << nl
            << "    backward" << nl
            << "    NewmarkBeta" << nl
            << abort(FatalError);
    }

    return tresult;
}


// tmp<pointVectorField> laplacian
// (
//     const tensor& gamma,
//     const pointVectorField& pf,
//     const fvMesh& dualMesh,
//     const labelList& dualCellToPoint
// )
// {
//     // Take references to the meshes
//     const pointMesh& pMesh = pf.mesh();

//     // Prepare the result field
//     tmp<pointVectorField> tresult
//     (
//         new pointVectorField
//         (
//             IOobject
//             (
//                 "laplacian(gamma," + pf.name() + ")",
//                 pMesh.time().timeName(),
//                 pMesh.db(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             pMesh,
//             dimensionedVector("zero", dimForce, vector::zero),
//             "calculated"
//         )
//     );
// #ifdef OPENFOAMESIORFOUNDATION
//     pointVectorField& result = tresult.ref();
// #else
//     pointVectorField& result = tresult();
// #endif

//     // Take references
//     vectorField& resultI = result.internalField();
//     const vectorField& dualSfI = dualMesh.Sf();
//     const labelList& dualOwn = dualMesh.owner();
//     const labelList& dualNei = dualMesh.neighbour();
//     const vectorField& pfI = pf.internalField();
//     const pointField& points = pMesh.mesh().points();

//     // Loop over all dual faces in the dual mesh and add contributions to the
//     // primary mesh own and nei points associated with that dual face
//     forAll(dualOwn, dualFaceI)
//     {
//         // Dual cell owner of dualFaceI
//         const label dualOwnCellID = dualOwn[dualFaceI];

//         // Dual cell neighbour of dualFaceI
//         const label dualNeiCellID = dualNei[dualFaceI];

//         // Primary mesh point at the centre of dualOwnCellID
//         const label ownPointID = dualCellToPoint[dualOwnCellID];

//         // Primary mesh point at the centre of dualNeiCellID
//         const label neiPointID = dualCellToPoint[dualNeiCellID];

//         // dualFaceI area vector
//         const vector& curDualSf = dualSfI[dualFaceI];
//         const scalar curMagDualSf = mag(curDualSf);

//         // Calculate contribute from the dual face to the own and nei cells
//         const scalar deltaCoeff =
//             curMagDualSf
//            /(curDualSf & (points[neiPointID] - points[ownPointID]));
//         const vector contrib
//         (
//             gamma & curMagDualSf*(pfI[neiPointID] - pfI[ownPointID])*deltaCoeff
//         );

//         // Add contribution to own and nei cells i.e. to own and nei points
//         resultI[ownPointID] += contrib;
//         resultI[neiPointID] -= contrib;
//     }

//     return tresult;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vfvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
