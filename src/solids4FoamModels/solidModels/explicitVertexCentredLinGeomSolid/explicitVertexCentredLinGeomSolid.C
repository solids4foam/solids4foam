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

#include "explicitVertexCentredLinGeomSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "sparseMatrix.H"
#include "symmTensor4thOrder.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#include "pointConstraints.H"
#include "sparseMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitVertexCentredLinGeomSolid, 0);
addToRunTimeSelectionTable(solidModel, explicitVertexCentredLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void explicitVertexCentredLinGeomSolid::updatePointDivSigma
(
    const pointVectorField& pointD,
    surfaceTensorField& dualGradDf,
    surfaceSymmTensorField& dualSigmaf,
    pointVectorField& pointDivSigma
)
{
    if (debug)
    {
        Info<< "void explicitVertexCentredLinGeomSolid::"
            << "updatePointDivSigma(...): start" << endl;
    }


    // Calculate gradD at dual faces
    dualGradDf = vfvc::fGrad
    (
        pointD,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta_,
        debug
    );

    // Calculate stress at dual faces
    dualMechanicalPtr_().correct(dualSigmaf);

    // Calculate the tractions at the dual faces
    surfaceVectorField dualTraction
    (
        (dualMesh().Sf()/dualMesh().magSf()) & dualSigmaf
    );

    // Enforce extract tractions on traction boundaries
    enforceTractionBoundaries
    (
        pointD, dualTraction, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
#else
            dualTraction.boundaryField()[patchI] = vector::zero;
#endif
        }
    }

    // Calculate divergence of stress for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh().magSf());

    // Map dual cell field to primary mesh point field
    vectorField& pointDivSigmaI = pointDivSigma;
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigmaI[pointID] = dualDivSigma[dualCellI]/pointVol_[pointID];
    }

    // What happens at processor boundaries? pointDivSigma needs to be
    // calculated consistently at processor boundaries: to be fixed!
    // For points on proc boundaries, I need to sum pointDivSigma
    pointConstraints::syncUntransformedData
    (
        mesh(), pointDivSigma, plusEqOp<vector>()
    );

    forAll(pointDivSigmaI, pointI)
    {
        pointDivSigmaI[pointI] *= pointVol_[pointI];
    }
    // Something is still not right at these points: maybe the problem is with
    // volumes: I tried converting to force then syncing then converting back to
    // force per unit volume, but this also doesn't work.
    // To-do: run a simple case serial and parallel and print out the values
    // if (Pstream::parRun())
    // {
    //     FatalError
    //         << "To fixed for parallel running" << abort(FatalError);
    // }

    if (debug)
    {
        Info<< "void explicitVertexCentredLinGeomSolid::"
            << " updatePointDivSigma(...): end" << endl;
    }
}


void explicitVertexCentredLinGeomSolid::enforceTractionBoundaries
(
    const pointVectorField& pointD,
    surfaceVectorField& dualTraction,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    const pointMesh& pMesh = pointD.mesh();
    const fvMesh& dualMesh = dualTraction.mesh();

    forAll(pointD.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const solidTractionPointPatchVectorField& tracPatch =
                refCast<const solidTractionPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                );

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // Primary mesh point normals
            const vectorField& n =
                pMesh.boundary()[patchI].pointNormals();

            // Primary mesh point tractions
            const vectorField totalTraction
            (
                tracPatch.traction() - n*tracPatch.pressure()
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
            forAll(totalTraction, pI)
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
                                totalTraction[pI];

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
                    "void explicitVertexCentredLinGeomSolid::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = " << nPointsPerDualFace
                    << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;

            // Overwrite the dual patch face traction
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] = dualFaceTraction;
#else
            dualTraction.boundaryField()[patchI] = dualFaceTraction;
#endif
        }
        else if
        (
            isA<symmetryPointPatchVectorField>(pointD.boundaryField()[patchI])
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            // Set the dual patch face shear traction to zero
            const vectorField n(dualMesh.boundary()[patchI].nf());
#ifdef OPENFOAMESIORFOUNDATION
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#else
            dualTraction.boundaryField()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#endif
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitVertexCentredLinGeomSolid::explicitVertexCentredLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    zeta_(solidModelDict().lookupOrDefault<scalar>("zeta", 1.0)),
    dampingCoeff_(solidModelDict().lookup("dampingCoeff")),
    dualMechanicalPtr_
    (
        new dualMechanicalModel
        (
            dualMesh(),
            nonLinGeom(),
            incremental(),
            mechanical(),
            dualMeshMap().dualFaceToCell()
        )
    ),
    twoD_(sparseMatrixTools::checkTwoD(mesh())),
    pointU_
    (
        IOobject
        (
            "pointU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pointA_
    (
        IOobject
        (
            "pointA",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity/dimTime, vector::zero)
    ),
    pointRho_
    (
        IOobject
        (
            "point(rho)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    ),
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
        pMesh(),
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    pointDivSigma_
    (
        IOobject
        (
            "pointDivSigma",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimForce/dimVolume, vector::zero)
    ),
    dualGradDf_
    (
        IOobject
        (
            "grad(D)f",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualSigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        "calculated"
    ),
    globalPointIndices_(mesh())
#ifdef OPENFOAMESI
    ,
    pointVolInterp_(pMesh(), mesh())
#endif
{
    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // pointD field must be defined
    pointDisRequired();

    // Set point density field
    mechanical().volToPoint().interpolate(rho(), pointRho_);


    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
#ifdef OPENFOAMESIORFOUNDATION
    scalarField& pointVolI = pointVol_.primitiveFieldRef();
#else
    scalarField& pointVolI = pointVol_.internalField();
#endif
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
    }

    // Sum the shared point volumes for consistency with serial runs
    pointConstraints::syncUntransformedData
    (
        mesh(), pointVol_, plusEqOp<scalar>()
    );

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointU_.oldTime().storeOldTime();
    pointA_.storeOldTime();


    // Write the compact edge gradient factor
    Info<< "zeta: " << zeta_ << endl;

    // Disable the writing of the unused fields
    D().writeOpt() = IOobject::NO_WRITE;
    D().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    DD().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    U().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

explicitVertexCentredLinGeomSolid::~explicitVertexCentredLinGeomSolid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitVertexCentredLinGeomSolid::setDeltaT(Time& runTime)
{
    // Max wave speed in the domain
    const scalar waveSpeed = max
    (
        Foam::sqrt(mechanical().impK()/mechanical().rho())
    ).value();

    // deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e. deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // stableTimeStep/sqrt(2): we will default to this value
    const scalar requiredDeltaT =
        1.0/
        gMax
        (
#ifdef OPENFOAMESIORFOUNDATION
            DimensionedField<scalar, Foam::surfaceMesh>
#else
            Field<scalar>
#endif
            (
                dualMesh().surfaceInterpolation::deltaCoeffs().internalField()
               *waveSpeed
            )
        );

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    // Update print info
    physicsModel::printInfo() = bool
    (
        runTime.timeIndex() % infoFrequency() == 0
     || mag(runTime.value() - runTime.endTime().value()) < SMALL
    );

    physicsModel::printInfo() = false;

    if (time().timeIndex() == 1)
    {
        Info<< nl << "Setting deltaT = " << newDeltaT
            << ", maxCo = " << maxCo << endl;
    }

    runTime.setDeltaT(newDeltaT);
}


bool explicitVertexCentredLinGeomSolid::evolve()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the solid momentum equation for pointD" << nl
            << "Simulation Time, Clock Time, Max Stress" << endl;
    }

    physicsModel::printInfo() = bool
    (
        time().timeIndex() % infoFrequency() == 0
     || mag(time().value() - time().endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< time().value() << " " << time().elapsedClockTime()
            << " " << max(mag(dualSigmaf_)).value() << endl;

        physicsModel::printInfo() = false;
    }

    // Central difference scheme

    // Take a reference to the current and previous time-step
    const dimensionedScalar& deltaT = time().deltaT();
    //const dimensionedScalar& deltaT0 = time().deltaT0();

    // Compute the velocity
    // Note: this is the velocity at the middle of the time-step
    //pointU_ = pointU_.oldTime() + 0.5*(deltaT + deltaT0)*pointA_.oldTime();
    pointU_ = pointU_.oldTime() + deltaT*pointA_.oldTime();

    // Compute displacement
    pointD() = pointD().oldTime() + deltaT*pointU_;

    // Enforce boundary conditions on the displacement field
    pointD().correctBoundaryConditions();

    // Remove displacement in the empty directions
    forAll(mesh().geometricD(), dirI)
    {
        if (mesh().geometricD()[dirI] < 0)
        {
            pointD().primitiveFieldRef().replace(dirI, 0.0);
        }
    }

    // Update the divergence of stress based on the latest pointD field
    updatePointDivSigma(pointD(), dualGradDf_, dualSigmaf_, pointDivSigma_);

    // Compute acceleration
    pointA_ = pointDivSigma_/pointRho_ - dampingCoeff_*pointU_ + g();

    // Check energies
    // To-do

    return true;
}


void explicitVertexCentredLinGeomSolid::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    // Get point field on patch
    const vectorField traction
    (
        globalPatches()[interfaceI].globalPointToPatch
        (
            globalPatches()[interfaceI].interpolator().faceToPointInterpolate
            (
                faceZoneTraction
            )()
        )
    );

    // Lookup point patch field
#ifdef OPENFOAMESIORFOUNDATION
    pointPatchVectorField& ptPatch = pointD().boundaryFieldRef()[patchID];
#else
    pointPatchVectorField& ptPatch = pointD().boundaryField()[patchID];
#endif

    if (isA<solidTractionPointPatchVectorField>(ptPatch))
    {
        solidTractionPointPatchVectorField& patchD =
            refCast<solidTractionPointPatchVectorField>(ptPatch);

        patchD.traction() = traction;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::explicitVertexCentredLinGeomSolid::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << ptPatch.type()
            << " for point patch " << ptPatch.patch().name()
            << " should instead be type "
            << solidTractionPointPatchVectorField::typeName
            << abort(FatalError);
    }
}


void explicitVertexCentredLinGeomSolid::writeFields(const Time& runTime)
{
    Info<< nl << "Writing fields to " << runTime.timeName() << endl;

    // Update the increment of displacement
    pointDD() = pointD() - pointD().oldTime();

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1)
    {
        mechanical().mapGradToSubMeshes(gradD());
    }

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate
    mechanical().correct(sigma());

#ifdef OPENFOAMESI
    // Interpolate pointD to D
    // This is useful for visualisation but it is also needed when using preCICE
    pointVolInterp_.interpolate(pointD(), D());
#endif

//     // Calculate gradD at the primary points using least squares: this should
//     // be second-order accurate (... I think).
//     const pointTensorField pGradD(vfvc::pGrad(pointD(), mesh()));

//     // Calculate strain at the primary points based on pGradD
//     // Note: the symm operator is not defined for pointTensorFields so we will
//     // do it manually
//     // const pointSymmTensorField pEpsilon("pEpsilon", symm(pGradD));
//     pointSymmTensorField pEpsilon
//     (
//         IOobject
//         (
//             "pEpsilon",
//             runTime.timeName(),
//             runTime,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         pMesh(),
//         dimensionedSymmTensor("0", dimless, symmTensor::zero)
//     );

// #ifdef FOAMEXTEND
//     pEpsilon.internalField() = symm(pGradD.internalField());
// #else
//     pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
// #endif
//     pEpsilon.write();

//     // Equivalent strain at the points
//     pointScalarField pEpsilonEq
//     (
//         IOobject
//         (
//             "pEpsilonEq",
//             runTime.timeName(),
//             runTime,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         pMesh(),
//         dimensionedScalar("0", dimless, 0.0)
//     );

// #ifdef FOAMEXTEND
//     pEpsilonEq.internalField() =
//         sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
// #else
//     pEpsilonEq.primitiveFieldRef() =
//         sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
// #endif
//     pEpsilonEq.write();

//     Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << nl << endl;

    solidModel::writeFields(runTime);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
