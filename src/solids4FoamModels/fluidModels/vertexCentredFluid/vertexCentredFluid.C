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
#include "fvcDiv.H"
#include "sparseMatrixTools.H"
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

tmp<vectorField> vertexCentredFluid::residualMomentum
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


void vertexCentredFluid::setFixedDofs
(
    const pointVectorField& pointU,
    boolList& fixedDofs,
    pointField& fixedDofValues,
    symmTensorField& fixedDofDirections
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
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed value is
                    // consistent with the new one
                    if
                    (
                        mag
                        (
                            fixedDofDirections[pointID]
                          & (fixedDofValues[pointID] - val)
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
                    fixedDofDirections[pointID] = symmTensor(I);
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = val;
                    fixedDofDirections[pointID] = symmTensor(I);
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
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the current condition
                    if
                    (
                        mag
                        (
                            (pointNormals[pI] & fixedDofValues[pointID])
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
                    if (mag(fixedDofDirections[pointID] - symmTensor(I)) > 0)
                    {
                        // If the directions are orthogonal we can add them
                        const symmTensor curDir = sqr(pointNormals[pI]);
                        if (mag(fixedDofDirections[pointID] & curDir) > 0)
                        {
                            FatalError
                                << "Point " << pointID << " is fixed in two "
                                << "directions: this is only implemented for "
                                << "Cartesian axis directions" << abort(FatalError);
                        }

                        fixedDofDirections[pointID] += curDir;
                    }
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = vector::zero;
                    fixedDofDirections[pointID] = sqr(pointNormals[pI]);
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
    scalar& initResidual,
    const scalar res,
    const label nInterations,
    const pointVectorField& pointU,
    const vectorField& pointUcorr
) const
{
    // Calculate the residual as the root mean square of the correction
    const scalar residualAbs = gSum(magSqr(pointUcorr));

    // Store initial residual
    if (iCorr == 0)
    {
        initResidual = residualAbs;

        // If the initial residual is small then convergence has been achieved
        if (false) //initResidual < SMALL)
        {
            Info<< "    Initial residual is less than " << SMALL << nl
                << "    Converged" << endl;
            return true;
        }
        Info<< "    Initial residual = " << initResidual << endl;
    }

    // Define a normalised residual wrt the initial residual
    const scalar residualNorm = residualAbs/initResidual;

    // Calculate the maximum displacement
    const scalar maxMagD = gMax(mag(pointU.primitiveField()));

    // Check for convergence
    bool converged = false;
    if (residualNorm < solutionTol_)
    {
        Info<< "    Converged" << endl;
        converged = true;
    }
    // else if (residualAbs < SMALL)
    // {
    //     Info<< "    Converged: absolute residual is less than " << SMALL
    //         << endl;
    //     converged = true;
    // }

    if (iCorr == 0)
    {
        Info<< "    Corr, res, relRes, resAbs, iters, maxMag" << pointU.name()
            << endl;
    }

    if (iCorr % infoFrequency_ == 0 || converged || iCorr >= nCorr_ - 1)
    {
        Info<< "    " << iCorr
            << ", " << res
            << ", " << residualNorm
            << ", " << residualAbs
            << ", " << nInterations
            << ", " << maxMagD << endl;

        if (iCorr >= nCorr_)
        {
            Warning
                << "Max iterations reached within the momentum loop"
                << endl;
            converged = true;
        }
    }

    return converged;
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
    mu_(fluidProperties().lookup("mu")),
    rho_(fluidProperties().lookup("rho")),
    fullNewton_(fluidProperties().lookup("fullNewton")),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    twoDCorrector_(mesh()),
    fixedDofs_(mesh().nPoints(), false),
    fixedDofValues_(fixedDofs_.size(), vector::zero),
    fixedDofDirections_(fixedDofs_.size(), symmTensor::zero),
    fixedDofScale_
    (
        fluidProperties().lookupOrDefault<scalar>
        (
            "fixedDofScale",
            (mu_*Foam::sqrt(gAverage(mesh().magSf()))).value()
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
    setFixedDofs(pointU_, fixedDofs_, fixedDofValues_, fixedDofDirections_);

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
    Info<< "zeta: " << zeta_ << endl;
    Info<< "fixedDofScale: " << fixedDofScale_ << endl;

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
    sparseMatrix matrix(sum(globalPointIndices_.stencilSize()));

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
    const scalarField& pointVolI = pointVol_.internalField();

    // Point density
    //const scalarField& pointRhoI = pointRho_.internalField();

    // Solution field: point velocity correction
    vectorField pointUcorr(pointU_.internalField().size(), vector::zero);

    // Newton-Raphson loop over momentum equation
    int iCorr = 0;
    scalar initResidual = 0.0;
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
            matrix +=
                vfvm::ddt
                (
                    mesh().ddtScheme("ddt(pointU)"),
                    pointU_
                )()*(rho_.value()*pointVolI);

            // Add div(rho*U*U) coefficients
            matrix += vfvm::div(pointU_, dualMesh_)()*rho_.value();

            // Add div(s) coefficients
            // where div(s) == div(mu*grad(U)) == mu*laplacian(U)
            matrix -= vfvm::laplacian(pointU_, dualMesh_, zeta_)()*mu_.value();
        }

        if (debug > 1)
        {
            matrix.print();
        }

        // Calculate the source as the negative of the momentum residual
        vectorField source(-residualMomentum(pointU_, pointP_));

        // Enforce fixed DOF on the linear system
        sparseMatrixTools::enforceFixedDof
        (
            matrix,
            source,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofValues_,
            fixedDofScale_
        );

        if (debug > 1)
        {
            matrix.print();
        }

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
            solverPerf = sparseMatrixTools::solveLinearSystemPETSc
            (
                matrix,
                source,
                pointUcorr,
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
            sparseMatrixTools::solveLinearSystemEigen
            (
                matrix, source, pointUcorr, twoD_, false, debug
            );
        }

        if (debug)
        {
            Info<< "bool vertexCentredFluid::evolve(): "
                << " solving linear system: end" << endl;
        }

        // Update point displacement field
        if (mesh().relaxField(pointU_.name()))
        {
            // Relaxing the correction can help convergence

            const scalar rf
            (
                mesh().fieldRelaxationFactor(pointU_.name())
            );
            pointU_.primitiveFieldRef() += rf*pointUcorr;
        }
        else
        {
            pointU_.primitiveFieldRef() += pointUcorr;
        }
        pointU_.correctBoundaryConditions();
    }
    while
    (
        !converged
        (
            iCorr,
            initResidual,
            mag(solverPerf.finalResidual()),
            cmptMax(solverPerf.nIterations()),
            pointU_,
            pointUcorr
        ) && ++iCorr
    );

    // // Calculate gradD at dual faces
    // dualGradDf_ = vfvc::fGrad
    // (
    //     pointD(),
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
    // gradD() = vfvc::grad(pointD(), mesh());

    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    // pointVolInterp_.interpolate(pointD(), D());

    return true;
}


void vertexCentredFluid::writeFields(const Time& runTime)
{
    // // Calculate cell gradient
    // // This assumes a constant gradient within each primary mesh cell
    // // This is a first-order approximation
    // gradD() = vfvc::grad(pointD(), mesh());

    // // Calculate gradD at the primary points using least squares: this should
    // // be second-order accurate (... I think).
    // const pointTensorField pGradD(vfvc::pGrad(pointD(), mesh()));

    fluidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // OPENFOAM_COM

// ************************************************************************* //
