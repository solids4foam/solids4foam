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

#ifdef OPENFOAM_COM

#include "vertexCentredFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrix.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "vfvmCellPointExtended.H"
#include "fvcDiv.H"
#include "sparseMatrixTools.H"
#include "sparseMatrixExtendedTools.H"
#include "fixedValuePointPatchFields.H"
#include "symmetryPointPatchFields.H"
#include "slipPointPatchFields.H"
#include "pointFieldFunctions.H"
#ifdef USE_PETSC
    #include <petscksp.h>
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredFluid, 0);
addToRunTimeSelectionTable(fluidModel, vertexCentredFluid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

tmp<vectorField> vertexCentredFluid::residualU
(
    const pointVectorField& pointU,
    const pointScalarField& pointP
)
{
    if (debug)
    {
        Info<< "void vertexCentredFluid::residualMomentum(...): start"
            << endl;
    }

    // Prepare the result field
    tmp<vectorField> tresidual(new vectorField(pointU.size(), vector::zero));
    vectorField& residual = tresidual.ref();

    // The residual F is defined as
    // F = rho*ddt(U) + rho*div(U*U) - div(sigma)
    //   = rho*ddt(U) + rho*div(U*U) - div(s - p*I)
    //   = rho*ddt(U) + rho*div(U*U) - div(mu*grad(U) - p*I)
    // where density is assumed constant and uniform.
    // Take care that the source of the linear system will be the negative of F

    // Take references to internal fields for clarity and brevity
    const scalarField& pointVolI = pointVol_.internalField();
    //const scalarField& pointRhoI = pointRho_.internalField();


    // Add transient term
    residual +=
        rho_.value()*pointVolI
       *vfvc::ddt(mesh().ddtScheme("ddt(pointU)"), pointU);

    // Calculate the advection term div(rho*U*U)

    // Interpolate pointU to the dual faces
    const surfaceVectorField dualUf
    (
        "dualUf",
        vfvc::interpolate(pointU, mesh(), dualMesh_)
    );

    // Calculate the flux through the dual faces
    surfaceScalarField dualFlux("dualFlux", dualMesh_.Sf() & dualUf);

    // Enforce the flux at the boundaries
    enforceFluxBoundaries
    (
        pointU, dualFlux, mesh(), dualMesh_.dualMeshMap().pointToDualFaces()
    );

    // Calculate advection divergence term for the dual cells
    const vectorField dualDivUU = fvc::div(dualFlux*dualUf);

    // Map dual cell field to the primary mesh point field
    vectorField pointDivUU(pointU.size(), vector::zero);
    const labelList& dualCellToPoint =
        dualMesh_.dualMeshMap().dualCellToPoint();
    forAll(dualDivUU, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivUU[pointID] = dualDivUU[dualCellI];
    }

    // Add the advection term to residual
    residual += rho_.value()*pointDivUU*pointVolI;


    // Calculate the surface forces term div(mu*gradU - p*I)

    // Calculate the velocity gradient at the dual faces
    const surfaceTensorField dualGradUf
    (
        "dualGradUf",
        vfvc::fGrad
        (
            pointU,
            mesh(),
            dualMesh_,
            dualMesh_.dualMeshMap().dualFaceToCell(),
            dualMesh_.dualMeshMap().dualCellToPoint(),
            zeta_,
            false // debug
        )
    );

    // Interpolate pointP to the dual faces
    const surfaceScalarField dualPf
    (
        "dualPf",
        vfvc::interpolate
        (
            pointP,
            mesh(),
            dualMesh_,
            dualMesh_.dualMeshMap().dualFaceToCell(),
            dualMesh_.dualMeshMap().dualCellToPoint(),
            false // debug
        )
    );

    // Dual unit normals
    const surfaceVectorField dualN(dualMesh_.Sf()/dualMesh_.magSf());

    // Calculate the tractions on the dual faces
    surfaceVectorField dualTraction
    (
        "dualTraction", mu_*(dualN & dualGradUf) - dualPf*dualN
    );

    // Enforce exact tractions on traction boundaries, e.g. outlets
    enforceTractionBoundaries
    (
        pointP,
        pointU,
        dualTraction,
        mesh(),
        dualMesh_.dualMeshMap().pointToDualFaces()
    );

    // Calculate divergence of stress for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh_.magSf());

    // Map dual cell field to primary mesh point field
    vectorField pointDivSigma(pointU.size(), vector::zero);
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigma[pointID] = dualDivSigma[dualCellI];
    }

    // Add the surface force term to residual
    residual -= pointDivSigma*pointVolI;

    if (debug)
    {
        Info<< "void vertexCentredFluid::residualMomentum(...): end"
            << endl;
    }

    return tresidual;
}


tmp<scalarField> vertexCentredFluid::residualP
(
    const pointVectorField& pointU,
    const pointScalarField& pointP
)
{
    if (debug)
    {
        Info<< "void vertexCentredFluid::residualMomentum(...): start"
            << endl;
    }

    // Prepare the result field
    tmp<scalarField> tresidual(new scalarField(pointU.size(), 0.0));
    scalarField& residual = tresidual.ref();

    // The residual for the pressure equation is:
    // F = p - gamma*laplacian(p) - k*div(U)

    // Take references to internal fields for clarity and brevity
    const scalarField& pointVolI = pointVol_.internalField();

    // Gradient of U at the points
    // We will calculate the div(U) term as tr(grad(U))
    const pointTensorField pointGradU(vfvc::pGrad(pointU_, mesh()));

    // Calculate laplacian(p) field
    const pointScalarField laplacianP
    (
        vfvc::laplacian
        (
            pointP,
            mesh(),
            dualMesh_,
            dualMesh_.dualMeshMap().dualFaceToCell(),
            dualMesh_.dualMeshMap().dualCellToPoint(),
            zeta_,
            debug
        )
    );

    // Calculate the K*div(U) field, where K presents a penalty bulk modulus
    const pointScalarField kDivU(-k_*tr(pointGradU));

    // Add pressure term
    residual += pointP*pointVolI;

    // Add gamma*laplacian(p) term
    residual -= pressureSmoothingCoeff_*laplacianP*pointVolI;

    // Add kDivU term
    residual += kDivU*pointVolI;

    return tresidual;
}


void vertexCentredFluid::setFixedUDofs
(
    const pointVectorField& pointU,
    boolList& fixedUDofs,
    pointField& fixedUDofValues,
    symmTensorField& fixedUDofDirections
) const
{
    // Flag all fixed DOFs

    forAll(pointU.boundaryField(), patchI)
    {
        if
        (
            isA<fixedValuePointPatchVectorField>
            (
                pointU.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointU.mesh().mesh().boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                const vector& val = pointU[pointID];

                // Check if this point has already been fixed
                if (fixedUDofs[pointID])
                {
                    // Check if the existing prescribed value is
                    // consistent with the new one
                    if
                    (
                        mag
                        (
                            fixedUDofDirections[pointID]
                          & (fixedUDofValues[pointID] - val)
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredFluid::setFixedDofs(...)"
                        )   << "Inconsistent values prescribed at point "
                            << "= " << pointU.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // Set all directions as fixed, just in case it was
                    // previously marked as a symmetry point
                    fixedUDofDirections[pointID] = symmTensor(I);
                }
                else
                {
                    fixedUDofs[pointID] = true;
                    fixedUDofValues[pointID] = val;
                    fixedUDofDirections[pointID] = symmTensor(I);
                }
            }
        }
        else if
        (
            isA<symmetryPointPatchVectorField>
            (
                pointU.boundaryField()[patchI]
            )
         || isA<slipPointPatchVectorField>
            (
                pointU.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointU.mesh().boundary()[patchI].meshPoints();
            const vectorField& pointNormals =
                pointU.mesh().boundary()[patchI].pointNormals();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                // Check if this point has already been fixed
                if (fixedUDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the current condition
                    if
                    (
                        mag
                        (
                            (pointNormals[pI] & fixedUDofValues[pointID])
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredFluid::setFixedDofs(...)"
                        )   << "Inconsistent values prescribed at point "
                            << "= " << pointU.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // If the point is not fully fixed then make sure the normal
                    // direction is fixed
                    if (mag(fixedUDofDirections[pointID] - symmTensor(I)) > 0)
                    {
                        // If the directions are orthogonal we can add them
                        const symmTensor curDir = sqr(pointNormals[pI]);
                        if (mag(fixedUDofDirections[pointID] & curDir) > 0)
                        {
                            FatalError
                                << "Point " << pointID << " is fixed in two "
                                << "directions: this is only implemented for "
                                << "Cartesian axis directions" << abort(FatalError);
                        }

                        fixedUDofDirections[pointID] += curDir;
                    }
                }
                else
                {
                    fixedUDofs[pointID] = true;
                    fixedUDofValues[pointID] = vector::zero;
                    fixedUDofDirections[pointID] = sqr(pointNormals[pI]);
                }
            }
        }
    }
}


void vertexCentredFluid::setFixedPDofs
(
    const pointScalarField& pointP,
    boolList& fixedPDofs,
    scalarField& fixedPDofValues
) const
{
    // Flag all fixed DOFs

    forAll(pointP.boundaryField(), patchI)
    {
        if
        (
            isA<fixedValuePointPatchScalarField>
            (
                pointP.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointP.mesh().mesh().boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                const scalar val = pointP[pointID];

                // Check if this point has already been fixed
                if (fixedPDofs[pointID])
                {
                    // Check if the existing prescribed value is
                    // consistent with the new one
                    if (mag(fixedPDofValues[pointID] - val)  > SMALL)
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredFluid::setFixedPDofs(...)"
                        )   << "Inconsistent values prescribed at point "
                            << "= " << pointP.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }
                }
                else
                {
                    fixedPDofs[pointID] = true;
                    fixedPDofValues[pointID] = val;
                }
            }
        }
    }
}


void vertexCentredFluid::enforceFluxBoundaries
(
    const pointVectorField& pointU,
    surfaceScalarField& dualFlux,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    const fvMesh& dualMesh = dualFlux.mesh();

    forAll(pointU.boundaryField(), patchI)
    {
        if
        (
            isA<fixedValuePointPatchVectorField>
            (
                pointU.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // Calculate the point patch flux field
            const scalarField patchPointFlux
            (
                mesh.boundaryMesh()[patchI].pointNormals()
              & pointU.boundaryField()[patchI].patchInternalField()
            );

            // Initialise dual mesh faces flux field
            scalarField dualFaceFlux
            (
                dualMesh.boundaryMesh()[patchI].size(), 0.0
            );

            // Multiple points map to each dual face so we will count them
            // and then divide the dualFaceFlux by this field so that it is
            // the average of all the points that map to it
            scalarField nPointsPerDualFace(dualFaceFlux.size(), 0.0);

            // Map from primary mesh point field to dual mesh face field using
            // the pointToDualFaces map
            forAll(patchPointFlux, pI)
            {
                const label pointID = meshPoints[pI];
                const labelList& curDualFaces = pointToDualFaces[pointID];

                forAll(curDualFaces, dfI)
                {
                    const label dualFaceID = curDualFaces[dfI];

                    if (!dualMesh.isInternalFace(dualFaceID))
                    {
                        // Check which patch this dual face belongs to
                        const label dualPatchID =
                            dualMesh.boundaryMesh().whichPatch(dualFaceID);

                        if (dualPatchID == patchI)
                        {
                            // Find local face index
                            const label localDualFaceID =
                                dualFaceID
                              - dualMesh.boundaryMesh()[dualPatchID].start();

                            // Set dual face flux
                            dualFaceFlux[localDualFaceID] += patchPointFlux[pI];

                            // Update the count for this face
                            nPointsPerDualFace[localDualFaceID]++;
                        }
                    }
                }
            }

            if (gMin(nPointsPerDualFace) < 1)
            {
                FatalErrorIn
                (
                    "void vertexCentredLinGeomSolid::"
                    "enforceFluxBoundaries(...)"
                )   << "Problem setting fluxs: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceFlux /= nPointsPerDualFace;

            // Overwrite the dual patch face flux
            dualFlux.boundaryFieldRef()[patchI] = dualFaceFlux;
        }
        else if
        (
            isA<symmetryPointPatchVectorField>(pointU.boundaryField()[patchI])
        )
        {
            // Set the dual patch flux to zero
            dualFlux.boundaryFieldRef()[patchI] = 0.0;
        }
    }

    // Set coupled boundary (e.g. processor) flux fields to zero: this
    // ensures their global contribution is zero
    forAll(dualFlux.boundaryField(), patchI)
    {
        if (dualFlux.boundaryField()[patchI].coupled())
        {
            dualFlux.boundaryFieldRef()[patchI] = 0.0;
        }
    }
}


void vertexCentredFluid::enforceTractionBoundaries
(
    const pointScalarField& pointP,
    const pointVectorField& pointU,
    surfaceVectorField& dualTraction,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    const pointMesh& pMesh = pointP.mesh();
    const fvMesh& dualMesh = dualTraction.mesh();

    forAll(pointP.boundaryField(), patchI)
    {
        if
        (
            isA<fixedValuePointPatchScalarField>
            (
                pointP.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // Primary mesh point normals
            const vectorField& n =
                pMesh.boundary()[patchI].pointNormals();

            // We assume the traction is equal to n*p on a boundary where
            // pressure is prescribed
            const vectorField patchPointTraction
            (
                -pointP.boundaryField()[patchI].patchInternalField()*n
            );

            // Create dual mesh faces traction field
            vectorField dualFaceTraction
            (
                dualMesh.boundaryMesh()[patchI].size(), vector::zero
            );

            // Multiple points map to each dual face so we will count them
            // and then divide the dualFaceTraction by this field so that it is
            // the average of all the points that map to it
            scalarField nPointsPerDualFace(dualFaceTraction.size(), 0.0);

            // Map from primary mesh point field to second mesh face field using
            // the pointToDualFaces map
            forAll(patchPointTraction, pI)
            {
                const label pointID = meshPoints[pI];
                const labelList& curDualFaces = pointToDualFaces[pointID];

                forAll(curDualFaces, dfI)
                {
                    const label dualFaceID = curDualFaces[dfI];

                    if (!dualMesh.isInternalFace(dualFaceID))
                    {
                        // Check which patch this dual face belongs to
                        const label dualPatchID =
                            dualMesh.boundaryMesh().whichPatch(dualFaceID);

                        if (dualPatchID == patchI)
                        {
                            // Find local face index
                            const label localDualFaceID =
                                dualFaceID
                              - dualMesh.boundaryMesh()[dualPatchID].start();

                            // Set dual face traction
                            dualFaceTraction[localDualFaceID] +=
                                patchPointTraction[pI];

                            // Update the count for this face
                            nPointsPerDualFace[localDualFaceID]++;
                        }
                    }
                }
            }

            if (gMin(nPointsPerDualFace) < 1)
            {
                FatalErrorIn
                (
                    "void vertexCentredFluid::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;

            // Overwrite the dual patch face traction
            dualTraction.boundaryFieldRef()[patchI] = dualFaceTraction;
        }
        else if
        (
            isA<symmetryPointPatchScalarField>(pointP.boundaryField()[patchI])
         || isA<slipPointPatchVectorField>
            (
                pointU.boundaryField()[patchI]
            )
        )
        {
            // Set the dual patch face shear traction to zero
            const vectorField n(dualMesh.boundary()[patchI].nf());
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
        }
    }

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
        }
    }
}


bool vertexCentredFluid::vertexCentredFluid::converged
(
    const label iCorr,
    scalar& initResidualU,
    scalar& initResidualP,
    const scalar res,
    const label nInterations,
    const pointVectorField& pointU,
    const pointScalarField& pointP,
    const Field<RectangularMatrix<scalar>>& pointUPcorr
) const
{
    // Calculate the residuals as the root mean square of the correction
    // Velocity residual
    scalar residualUAbs = 0;
    scalar residualPAbs = 0;
    scalar maxUCorr = 0.0;
    scalar maxPCorr = 0.0;
    forAll(pointUPcorr, pointI)
    {
        // Velocity residual
        const vector curUCorr
        (
            pointUPcorr[pointI](0,0),
            pointUPcorr[pointI](1,0),
            pointUPcorr[pointI](2,0)
        );

        residualUAbs += magSqr(curUCorr);

        maxUCorr = max(maxUCorr, mag(curUCorr));

        // Pressure residual
        residualPAbs += sqr(pointUPcorr[pointI](3,0));

        maxPCorr = max(maxPCorr, mag(pointUPcorr[pointI](3,0)));
    }

    residualUAbs /= sqrt(residualUAbs/pointU.size());
    residualPAbs /= sqrt(residualPAbs/pointU.size());

    // Store initial residual
    if (iCorr == 0)
    {
        initResidualU = residualUAbs;
        initResidualP = residualPAbs;

        // If the initial residual is small then convergence has been achieved
        if (initResidualU < SMALL && initResidualP < SMALL)
        {
            Info<< "    Both displacement and pressure residuals are less than " << VSMALL
                        << "	Converged" << endl;
                return true;
        }
        Info<< "	Initial displacement residual = " << initResidualU << nl
            << "	Initial pressure residual = " << initResidualP << endl;
    }

    // Define a normalised residual wrt the initial residual
    const scalar residualUNorm = residualUAbs/initResidualU;
    const scalar residualPNorm = residualPAbs/initResidualP;

    // Calculate the maximum displacement
    const scalar maxMagU = gMax(mag(pointU.primitiveField()));

    // Calculate the maximum pressure
    const scalar maxMagP = gMax(mag(pointP.primitiveField()));

    // Print information for the displacement
    Info<< "	Iter = " << iCorr
        << ", relRef = " << residualUNorm
        << ", resAbs = " << residualUAbs
        << ", nIters = " << nInterations
        << ", maxU = " << maxMagU
        << ", maxUCorr = " << maxUCorr << endl;

    // Print information for the pressure
    Info<< "	Iter = " << iCorr
        << ", relRef = " << residualPNorm
        << ", resAbs = " << residualPAbs
        << ", nIters = " << nInterations
        << ", maxP = " << maxMagP
        << ", maxPCorr = " << maxPCorr << endl;

    // Velocity tolerance
    const scalar UTol
    (
        readScalar(fluidProperties().lookup("solutionUTolerance"))
    );

    // Pressure tolerance
    const scalar PTol
    (
        readScalar(fluidProperties().lookup("solutionPTolerance"))
    );

    // Check for convergence
    if (residualUNorm < UTol && residualPNorm < PTol)
    {
        Info<< "	Converged" << endl;
        return true;
    }
    else if (iCorr >= nCorr_ - 1)
    {
        if (nCorr_ > 1)
        {
            Warning
                << "Max iterations reached within the momentum Newton-Raphson "
                   "loop" << endl;
        }

        return true;
    }

    // Convergence has not been reached
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredFluid::vertexCentredFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    pMesh_(pointMesh::New(mesh())),
    dualMesh_(mesh(), fluidProperties().subDict("dualMesh")),
    zeta_(readScalar(fluidProperties().lookup("zeta"))),
    compactStencil_("compactPressureStencil", fluidProperties()),
    pressureSmoothingCoeff_
    (
        readScalar(fluidProperties().lookup("pressureSmoothingCoeff"))
    ),
    mu_(fluidProperties().lookup("mu")),
    k_(fluidProperties().lookup("k")),
    rho_(fluidProperties().lookup("rho")),
    fullNewton_(fluidProperties().lookup("fullNewton")),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    twoDCorrector_(mesh()),
    fixedUDofs_(mesh().nPoints(), false),
    fixedPDofs_(mesh().nPoints(), false),
    fixedUDofValues_(fixedUDofs_.size(), vector::zero),
    fixedPDofValues_(fixedUDofs_.size(), 0.0),
    fixedUDofDirections_(fixedUDofs_.size(), symmTensor::zero),
    fixedUDofScale_
    (
        fluidProperties().lookupOrDefault<scalar>
        (
            "fixedUDofScale",
            (mu_*Foam::sqrt(gAverage(mesh().magSf()))).value()
        )
    ),
    fixedPDofScale_
    (
        fluidProperties().lookupOrDefault<scalar>
        (
            "fixedPDofScale",
            Foam::sqrt(gAverage(mesh().magSf()))
        )
    ),
    pointU_
    (
        IOobject
        (
            "pointU",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
    ),
    pointP_
    (
        IOobject
        (
            "pointP",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh_
    ),
    // pointRho_
    // (
    //     IOobject
    //     (
    //         "point(rho)",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     pMesh_,
    //     dimensionedScalar(fluidProperties().lookup("rho"))
    // ),
    pointVol_
    (
        IOobject
        (
            "pointVolumes",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    solutionTol_(readScalar(fluidProperties().lookup("solutionTolerance"))),
    nCorr_(readInt(fluidProperties().lookup("nCorr"))),
    infoFrequency_(readInt(fluidProperties().lookup("infoFrequency"))),
    // dualGradDf_
    // (
    //     IOobject
    //     (
    //         "grad(D)f",
    //         runTime.timeName(),
    //         dualMesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     dualMesh_,
    //     dimensionedTensor("zero", dimless, tensor::zero),
    //     "calculated"
    // ),
    // dualSigmaf_
    // (
    //     IOobject
    //     (
    //         "sigmaf",
    //         runTime.timeName(),
    //         dualMesh_,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     dualMesh_,
    //     dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
    //     "calculated"
    // ),
    globalPointIndices_(mesh()),
    pointVolInterp_(pointMesh::New(mesh()), mesh())
{
    // Set fixed degree of freedom list
    setFixedUDofs(pointU_, fixedUDofs_, fixedUDofValues_, fixedUDofDirections_);
    setFixedPDofs(pointP_, fixedPDofs_, fixedPDofValues_);

    // Set point density field
    // mechanical().volToPoint().interpolate(rho(), pointRho_);

    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
    scalarField& pointVolI = pointVol_.primitiveFieldRef();
    // scalarField& pointGlobalVolI = pointGlobalVol_;
    const scalarField& dualCellVol = dualMesh_.V();
    const labelList& dualCellToPoint = dualMesh_.dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
        // pointGlobalVolI[pointID] = dualCellVol[dualCellI];
    }

    // // Sum the shared point volumes to create the point global volumes
    // pointConstraints::syncUntransformedData
    // (
    //     mesh(), pointGlobalVol_, plusEqOp<scalar>()
    // );

    // Store old time fields
    pointU_.storeOldTime();
    //pointP_.storeOldTime();

    // Write parameter values
    Info<< "zeta: " << zeta_ << nl
        << "fixedUDofScale: " << fixedUDofScale_ << nl
        << "fixedPDofScale: " << fixedPDofScale_ << endl;

    // Disable the writing of the unused fields
    U().writeOpt() = IOobject::NO_WRITE;
    p().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

vertexCentredFluid::~vertexCentredFluid()
{
#ifdef USE_PETSC
    PetscFinalize();
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool vertexCentredFluid::evolve()
{
    Info<< "Evolving fluid solver" << endl;

    // Update boundary conditions
    pointU_.correctBoundaryConditions();
    pointP_.correctBoundaryConditions();

    // Initialise matrix
    sparseMatrixExtended matrix(sum(globalPointIndices_.stencilSize()));

    // A 6x6 material tangent is overkill for a Newtonian incompressible fluid
    // where the viscosity (mu) is the only parameter; however, it will be
    // useful when including turbulence models. So we will leave this comment
    // here for now to remind us
    // Store material tangent field for dual mesh faces
    // Field<scalarSquareMatrix> materialTangent
    // (
    //     dualMechanicalPtr_().materialTangentFaceField()
    // );

    // Point volumes
    //const scalarField& pointVolI = pointVol_.internalField();

    // Point density
    //const scalarField& pointRhoI = pointRho_.internalField();

    // Solution field: point velocity vector (0, 1, 2) and point pressure (3)
    Field<RectangularMatrix<scalar>> pointUPcorr
    (
        mesh().nPoints(), RectangularMatrix<scalar>(4, 1, 0.0)
    );

    // Newton-Raphson loop over momentum equation
    int iCorr = 0;
    scalar initResidualU = 0.0;
    scalar initResidualP = 0.0;
    SolverPerformance<vector> solverPerf;
    do
    {
        if (iCorr == 0 || fullNewton_)
        {
            // Assemble the matrix once per outer iteration
            matrix.clear();

            // If we had a material tangent, we would update it here
            // materialTangent = dualMechanicalPtr_().materialTangentFaceField();

            // Todo: vfvm operators are inconsisent in that some return the
            // results per unit volume and some do not
            // Add ddt coefficients
            sparseMatrix matrixU;
            // matrixU +=
            //     vfvm::ddt
            //     (
            //         mesh().ddtScheme("ddt(pointU)"),
            //         pointU_
            //     )()*(rho_.value()*pointVolI);

            // Add div(rho*U*U) coefficients
            //matrixU += vfvm::div(pointU_, dualMesh_)()*rho_.value();

            // Add div(s) coefficients
            // where div(s) == div(mu*grad(U)) == mu*laplacian(U)
            matrixU -= vfvm::laplacian(pointU_, dualMesh_, zeta_)()*mu_.value();

            // Insert matrixU into matrix
            const sparseMatrixData& dataU = matrixU.data();
            for
            (
                sparseMatrixData::const_iterator iter = dataU.begin();
                iter != dataU.end();
                ++iter
            )
            {
                // Convert the coefficient from a tensor to a RectangularMatrix
                // and insert it into the matrix
                const tensor& coeffU = iter();
                const label blockRowI = iter.key()[0];
                const label blockColI = iter.key()[1];

                if (twoD_)
                {
                    RectangularMatrix<scalar> coeff(3, 3, 0.0);
                    coeff(0, 0) = coeffU.xx();
                    coeff(0, 1) = coeffU.xy();
                    coeff(1, 0) = coeffU.yx();
                    coeff(1, 1) = coeffU.yy();

                    matrix(blockRowI, blockColI) += coeff;
                }
                else
                {
                    RectangularMatrix<scalar> coeff(4, 4, 0.0);
                    coeff(0, 0) = coeffU.xx();
                    coeff(0, 1) = coeffU.xy();
                    coeff(0, 2) = coeffU.xz();
                    coeff(1, 0) = coeffU.yx();
                    coeff(1, 1) = coeffU.yy();
                    coeff(1, 2) = coeffU.yz();
                    coeff(2, 0) = coeffU.zx();
                    coeff(2, 1) = coeffU.zy();
                    coeff(2, 2) = coeffU.zz();

                    matrix(blockRowI, blockColI) += coeff;
                }
            }

            // Add pressure terms
            // Ideally this would return a matrix
            // Add laplacian coefficient to the pressure equation
            // Note: we have not yet implemented operator- so we will
            // flip the sign of the pressureSmoothingCoeff as a hack to flip
            // the sign of the laplacian
            vfvm::laplacian
            (
                 matrix,
                 compactStencil_,
                 mesh(),
                 dualMesh_,
                 dualMesh_.dualMeshMap().dualFaceToCell(),
                 dualMesh_.dualMeshMap().dualCellToPoint(),
                 -pressureSmoothingCoeff_, // note negative sign!
                 debug
            );

            // Add diagonal and div(U) coefficients to the pressure equation
            // Ideally this would return a matrix and be split into two
            // operators (Sp name is misleading as it includes div(U)
            const tensorField pBarSensitivity(mesh().nPoints(), I);
            vfvm::Sp
            (
                matrix,
                mesh(),
                dualMesh_.dualMeshMap().dualCellToPoint(),
                pointVol_,
                pBarSensitivity,
                debug
            );
        }

        // Calculate the source as the negative of the residuals
        Field<RectangularMatrix<scalar>> source
        (
            mesh().nPoints(), RectangularMatrix<scalar>(4, 1, 0)
        );
        {
            vectorField sourceU(-residualU(pointU_, pointP_));
            scalarField sourceP(-residualP(pointU_, pointP_));
            forAll(source, pointI)
            {
                source[pointI](0, 0) = sourceU[pointI].x();
                source[pointI](1, 0) = sourceU[pointI].y();
                source[pointI](2, 0) = sourceU[pointI].z();
                source[pointI](3, 0) = sourceP[pointI];
            }
        }

        // Enforce fixed DOF on the linear system

        sparseMatrixExtendedTools::enforceFixedDisplacementDof
        (
            matrix,
            source,
            twoD_,
            fixedUDofs_,
            fixedUDofDirections_,
            fixedUDofValues_,
            fixedUDofScale_
        );

        sparseMatrixExtendedTools::enforceFixedPressureDof
        (
            matrix,
            source,
            twoD_,
            fixedPDofs_,
            fixedPDofValues_,
            fixedPDofScale_
        );

        // Solve linear system for displacement correction
        if (debug)
        {
            Info<< "bool vertexCentredFluid::evolve(): "
                << " solving linear system: start" << endl;
        }

        if (Switch(fluidProperties().lookup("usePETSc")))
        {
#ifdef USE_PETSC
            fileName optionsFile(fluidProperties().lookup("optionsFile"));
            solverPerf = sparseMatrixExtendedTools::solveLinearSystemPETSc
            (
                matrix,
                source,
                pointUPcorr,
                twoD_,
                optionsFile,
                mesh().points(),
                globalPointIndices_.ownedByThisProc(),
                globalPointIndices_.localToGlobalPointMap(),
                globalPointIndices_.stencilSizeOwned(),
                globalPointIndices_.stencilSizeNotOwned(),
                fluidProperties().lookupOrDefault<bool>("debugPETSc", false)
            );
#else
            FatalErrorIn("vertexCentredFluid::evolve()")
                << "PETSc not available. Please set the PETSC_DIR environment "
                << "variable and re-compile solids4foam" << abort(FatalError);
#endif
        }
        else
        {
            // Use Eigen SparseLU direct solver
            sparseMatrixExtendedTools::solveLinearSystemEigen
            (
                matrix, source, pointUPcorr, twoD_, false, debug
            );
        }

        if (debug)
        {
            Info<< "bool vertexCentredFluid::evolve(): "
                << " solving linear system: end" << endl;
        }

        // Retrieve the solution
        vectorField& pointUI = pointU_.primitiveFieldRef();
        scalarField& pointPI = pointP_.primitiveFieldRef();
        forAll(pointUPcorr, pointI)
        {
            pointUI[pointI][vector::X] += pointUPcorr[pointI](0, 0);
            pointUI[pointI][vector::Y] += pointUPcorr[pointI](1, 0);
            pointUI[pointI][vector::Z] += pointUPcorr[pointI](2, 0);
            pointPI[pointI] += pointUPcorr[pointI](3, 0);
        }

        // Enforce boundary conditions
        pointU_.correctBoundaryConditions();
        pointP_.correctBoundaryConditions();
    }
    while
    (
        !converged
        (
            iCorr,
            initResidualU,
            initResidualP,
            mag(solverPerf.finalResidual()),
            cmptMax(solverPerf.nIterations()),
            pointU_,
            pointP_,
            pointUPcorr
        ) && ++iCorr
    );

    // // Calculate gradD at dual faces
    // dualGradDf_ = vfvc::fGrad
    // (
    //     pointU(),
    //     mesh(),
    //     dualMesh_,
    //     dualMesh_.dualMeshMap().dualFaceToCell(),
    //     dualMesh_.dualMeshMap().dualCellToPoint(),
    //     zeta,
    //     debug
    // );

    // // Calculate cell gradient
    // // This assumes a constant gradient within each primary mesh cell
    // // This is a first-order approximation
    // gradD() = vfvc::grad(pointU(), mesh());

    // Interpolate pointU to D
    // This is useful for visualisation but it is also needed when using preCICE
    // pointVolInterp_.interpolate(pointU(), D());

    return true;
}


void vertexCentredFluid::writeFields(const Time& runTime)
{
    // // Calculate cell gradient
    // // This assumes a constant gradient within each primary mesh cell
    // // This is a first-order approximation
    // gradD() = vfvc::grad(pointU(), mesh());

    // // Calculate gradD at the primary points using least squares: this should
    // // be second-order accurate (... I think).
    // const pointTensorField pGradD(vfvc::pGrad(pointU(), mesh()));

    fluidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // OPENFOAM_COM

// ************************************************************************* //
