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

#include "linGeomTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "compatibilityFunctions.H"
#ifndef FOAMEXTEND
    #include "hofvm.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomTotalDispSolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomTotalDispSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void linGeomTotalDispSolid::predict()
{
    Info<< "Applying linear predictor to D" << endl;

    // Predict D using previous time steps
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


void linGeomTotalDispSolid::enforceTractionBoundaries
(
    surfaceVectorField& traction,
    const volVectorField& D,
    const surfaceVectorField& n
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const vectorField& nPatch = n.boundaryField()[patchI];

            if (highOrderResidual())
            {
#ifndef FOAMEXTEND
                // Face quadrature points weights
                const CompactListList<scalar>& faceQuadWeights =
                    displacementLeastSquares().quadrature().faceQuadWeights();

                const surfaceScalarField& magSf = mesh().magSf();

                // Get value at patch faces quadrature points
                autoPtr<CompactListList<vector>> patchQuadraturePointsValue =
                    tracPatch.evaluateQuadrature();

                const CompactListList<vector>& quadratureValues =
                    patchQuadraturePointsValue();

                forAll(mesh().boundaryMesh()[patchI], faceI)
                {
                    const label start = mesh().boundaryMesh()[patchI].start();

                    // Get global face index
                    const label faceID = faceI + start;

                    // Get the number of quadrature points for this face
                    const label nPoints = faceQuadWeights[faceID].size();

                    // Loop over quadrature points and add their contribution
                    traction.boundaryFieldRef()[patchI][faceI] = vector::zero;
                    for (label pointI = 0; pointI < nPoints; ++pointI)
                    {
                        traction.boundaryFieldRef()[patchI][faceI] +=
                            quadratureValues[faceI][pointI]
                           *faceQuadWeights[faceID][pointI];
                    }
                    // Divide with area because we use physical weights
                    traction.boundaryFieldRef()[patchI][faceI] *=
                        (1.0/(magSf.boundaryField()[patchI][faceI]));
                }
#endif
            }
            else
            {
#ifdef OPENFOAM_NOT_EXTEND
                traction.boundaryFieldRef()[patchI] =
                    tracPatch.traction() - nPatch*tracPatch.pressure();
#else
                traction.boundaryField()[patchI] =
                    tracPatch.traction() - nPatch*tracPatch.pressure();
#endif
            }
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            // Unit normals
            const vectorField& nPatch = n.boundaryField()[patchI];

            // Set shear traction to zero
#ifdef OPENFOAM_NOT_EXTEND
            traction.boundaryFieldRef()[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
#else
            traction.boundaryField()[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
#endif
        }
    }
}


bool linGeomTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver using an implicit segregated approach"
        << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    if (predictor_ && newTimeStep())
    {
        predict();
    }

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    // Mesh update loop
    do
    {
        int iCorr = 0;
        scalar currentResidualNorm = 0;
        scalar initialResidualNorm = 0;
        scalar deltaXNorm = 0;
        scalar xNorm = 0;
        const convergenceParameters convParam =
            readConvergenceParameters(solidModelDict());

        Info<< "Solving the momentum equation for D" << endl;

        // Unit normal vectors at the faces
        const surfaceVectorField n(mesh().Sf()/mesh().magSf());

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Calculate raction vectors at the faces
            surfaceVectorField traction(n & fvc::interpolate(sigma()));

            // Add stabilisation to the traction
            // We add this before enforcing the traction condition as the
            // stabilisation is set to zero on traction boundaries
            momentumStabilisation().updateVector(D(), &gradD());
            traction += impKf_*momentumStabilisation().faceVector();

            // Enforce traction boundary conditions
            enforceTractionBoundaries(traction, D(), n);

            // Linear momentum equation total displacement form
            // Assemble the RHS in stages.
            tmp<fvVectorMatrix> tRhsEqn
            (
                fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
            );
            tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, D(), "laplacian(DD,D)");
            tmpRef(tRhsEqn) += fvc::div(mesh().magSf()*traction);
            tmpRef(tRhsEqn) += rho()*g();
#ifdef OPENFOAM_COM
            tmpRef(tRhsEqn) += fvOptions()(ds_, D());
#endif

            fvVectorMatrix DEqn
            (
                rho()*fvm::d2dt2(D()) == tRhsEqn
            );

            // Add damping
            if (dampingCoeff().value() > SMALL)
            {
                DEqn += dampingCoeff()*rho()*fvm::ddt(D());
            }

            // Under-relaxation the linear system
            DEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DEqn);

            // Solve the linear system and store the residual
            currentResidualNorm = mag(DEqn.solve().initialResidual());

            // Norm of the solution correction
            deltaXNorm =
                sqrt
                (
                    gSum
                    (
                        magSqr
                        (
#ifdef OPENFOAM_NOT_EXTEND
                            D().primitiveField()
                          - D().prevIter().primitiveField()
#else
                            D().internalField()
                          - D().prevIter().internalField()
#endif
                        )
                    )
                );

            // Norm of the solution
#ifdef OPENFOAM_NOT_EXTEND
            xNorm = sqrt(gSum(magSqr(D().primitiveField())));
#else
            xNorm = sqrt(gSum(magSqr(D().internalField())));
#endif

            // Store the initial residual
            if (iCorr == 0)
            {
                initialResidualNorm = currentResidualNorm;
                Info<< "Initial Residual Norm = " << initialResidualNorm << nl
                    << "Initial Solution Norm = " << xNorm << endl;
            }

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Update the momentum equation inverse diagonal field
            // This may be used by the mechanical law when calculating the
            // hydrostatic pressure
            const volScalarField DEqnA("DEqnA", DEqn.A());

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());
        }
        while
        (
            !checkConvergence
            (
                currentResidualNorm,
                initialResidualNorm,
                deltaXNorm,
                xNorm,
                ++iCorr,
                convParam
            )
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), gradD(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (solidModel::mesh().update());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


bool linGeomTotalDispSolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Map the D field to the SNES solution vector
        foamPetscSnesHelper::InsertFieldComponents<vector>
        (
#ifdef OPENFOAM_NOT_EXTEND
            D().primitiveFieldRef(),
#else
            D().internalField(),
#endif
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD()
          ? makeList<label>({0,1})
          : makeList<label>({0,1,2})
        );
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the D field
    vectorField& DI = D();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    D().correctBoundaryConditions();

    if (solvePressure())
    {
        // Map the PETSc solution to the p field
        // p is located in the 4th component
        scalarField& pI = p();
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            foamPetscSnesHelper::solution(),
            pI,
            blockSize_ - 1 // Location of p component
        );

        p().correctBoundaryConditions();
    }

    // Update gradient of displacement
    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        gradD() = displacementLeastSquares().grad(D());

        // Calculate the cell centre stress using run-time selectable
        // mechanical law
        mechanical().correct(sigma());
#endif
    }
    else
    {
        mechanical().grad(D(), gradD());
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


bool linGeomTotalDispSolid::evolveExplicit()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the solid momentum equation for D using an explicit "
            << "approach" << nl
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
            << " " << max(mag(sigma())).value() << endl;

        physicsModel::printInfo() = false;
    }

    // Take references for brevity and efficiency
    const fvMesh& mesh = solidModel::mesh();
    volVectorField& D = solidModel::D();
    volTensorField& gradD = solidModel::gradD();
    volVectorField& U = solidModel::U();
    volSymmTensorField& sigma = solidModel::sigma();
    const volScalarField& rho = solidModel::rho();

    // Central difference scheme

    // Take a reference to the current and previous time-step
    const dimensionedScalar& deltaT = time().deltaT();
    //const dimensionedScalar& deltaT0 = time().deltaT0();

    // Compute the velocity
    // Note: this is the velocity at the middle of the time-step
    //pointU_ = pointU_.oldTime() + 0.5*(deltaT + deltaT0)*pointA_.oldTime();
    U = U.oldTime() + deltaT*A_.oldTime();

    // Compute displacement
    D = D.oldTime() + deltaT*U;

    // Enforce boundary conditions on the displacement field
    D.correctBoundaryConditions();

    if (solidModel::twoD())
    {
        // Remove displacement in the empty directions
        forAll(mesh.geometricD(), dirI)
        {
            if (mesh.geometricD()[dirI] < 0)
            {
#ifdef OPENFOAM_NOT_EXTEND
                D.primitiveFieldRef().replace(dirI, 0.0);
#else
                D.internalField().replace(dirI, 0.0);
#endif
            }
        }
    }

    // Update gradient of displacement
    mechanical().grad(D, gradD);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma);

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Calculate the traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the
    // stabilisation is set to zero on traction boundaries
    momentumStabilisation().updateVector(D, &gradD);
    traction += impKf_*momentumStabilisation().faceVector();

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // Solve the momentum equation for acceleration
    A_ = fvc::div(mesh.magSf()*traction)/rho
       + g()
       - dampingCoeff()*fvc::ddt(D);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomTotalDispSolid::linGeomTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "D",
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh(),
        solutionLocation::CELLS,
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
    ),
    impK_
    (
        solvePressure()
      ? 2.0*mechanical().shearModulus()
      : mechanical().impK()
    ),
    impKf_(fvc::interpolate(impK_)),
    rImpK_(1.0/impK_),
    rKappa_(1.0/mechanical().bulkModulus()),
    A_
    (
        IOobject
        (
            "A",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength/pow(dimTime, 2), vector::zero)
    ),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    ),
    ds_
    (
        IOobject
        (
            "ds",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("ds", (dimForce/dimVolume)/dimVelocity, 1.0)
    )
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // For consistent restarts, we will calculate the gradient field
    D().correctBoundaryConditions();
    D().storePrevIter();
    mechanical().grad(D(), gradD());

    Info<< "solvePressure = " << solvePressure() << endl;

    if (solvePressure())
    {
        if (solutionAlg() != solutionAlgorithm::PETSC_SNES)
        {
            FatalErrorInFunction
                << "The solution algorithm must be "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " when solvePressure is enabled" << abort(FatalError);
        }

        // Ensure p is created
        p();
    }

    if (predictor_)
    {
        // Construct U while the run time still points at the restart directory
        U();

        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // Check the gradScheme
    const word gradDScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh().gradScheme("grad(" + D().name() +')')
#else
        mesh().schemesDict().gradScheme("grad(" + D().name() +')')
#endif
    );

    if
    (
        solutionAlg() == solutionAlgorithm::PETSC_SNES
     || solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED
    )
    {
        if (gradDScheme != "leastSquaresS4f")
        {
            WarningIn(type() + "::" + type())
                << "Typically, the `leastSquaresS4f` gradScheme should be used "
                << "for `grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " and "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
                   ]
                << " solution algorithms" << endl;
        }

        // Set extrapolateValue to true for solidTraction boundaries
        forAll(D().boundaryField(), patchI)
        {
            if
            (
                isA<solidTractionFvPatchVectorField>
                (
                    D().boundaryField()[patchI]
                )
            )
            {
                Info<< "    Setting `extrapolateValue` to `true` on the "
                    << mesh().boundary()[patchI].name() << " patch of the D "
                    << "field" << endl;

                solidTractionFvPatchVectorField& tracPatch =
                    refCast<solidTractionFvPatchVectorField>
                    (
#ifdef OPENFOAM_NOT_EXTEND
                        D().boundaryFieldRef()[patchI]
#else
                        D().boundaryField()[patchI]
#endif
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void linGeomTotalDispSolid::setDeltaT(Time& runTime)
{
    if (solutionAlg() == solutionAlgorithm::EXPLICIT)
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
#ifdef OPENFOAM_NOT_EXTEND
                DimensionedField<scalar, Foam::surfaceMesh>
#endif
                (
                    mesh().surfaceInterpolation::
                    deltaCoeffs().internalField()
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
}


bool linGeomTotalDispSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    // else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    // {
    //     // Not yet implmented, although coupledUnsLinGeomLinearElasticSolid
    //     // could be combined with PETSc to achieve this.. todo!
    //     return evolveImplicitCoupled();
    // }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            // << solutionAlgorithmNames_.names() << endl;
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::EXPLICIT
               ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


#ifdef USE_PETSC

label linGeomTotalDispSolid::initialiseJacobian(Mat& jac)
{
#ifndef FOAMEXTEND
    if (highOrderJacobian())
    {
        return hofvm::initialiseJacobian
        (
            jac,
            *this,
            displacementLeastSquares(),
            D(),
            blockSize_
        );
    }
#endif

    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label linGeomTotalDispSolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label linGeomTotalDispSolid::formResidual
(
    Vec f,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();

    // Copy x into the D field
    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    D.correctBoundaryConditions();

    if (solvePressure() && highOrderResidual())
    {
        FatalErrorInFunction
                << "solvePressure must be disbled when using high order "
                << "residual calculation. Mixed approach not yet supported!"
                << abort(FatalError);
    }

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Traction vectors at the faces
    surfaceVectorField traction
    (
        IOobject
        (
            "traction",
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimPressure, vector::zero)
    );

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Update cell-centre gradient of displacement
        // Consider switching to mechanical().grad() interface
        gradD() = displacementLeastSquares().grad(D);

        // Update gradient of displacement at face quadrature points
        mechanical().grad(D, gradDQuad());

        // Calculate sigma at quadrature points
        mechanical().correct(gradDQuad(), sigmaQuad());

        // Integration over face quadrature points to get face traction
        traction = hofvc::surfaceIntegrate(sigmaQuad(), mesh);
#endif
    }
    else
    {
        // Update gradient of displacement
        mechanical().grad(D, gradD());

        // Enforce the boundary conditions again for any conditions that
        // use gradD
        //D.correctBoundaryConditions();

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }

    // Update velocity
    U() = fvc::ddt(D);

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        sigma() = dev(sigma()) - p*I;

        // Calculate the pressure gradient (we should store this!)
        const volVectorField gradp(fvc::grad(p));

        // Re-calculate the pressure stabilisation parameter
        pressureStabilisation().updateScalar(p, &gradp);

        // Dimensional consistency factor
        const dimensionedScalar one
        (
            "one", dimensionSet(-2, 4, 4, 0, 0, 0, 0), 1.0
        );

        // Compute the positive face-interpolated reciprocal of the approximate
        // momentum equation diagonal. This is the solid analogue of rAUf in
        // pressure-velocity coupling and has units of [Pa].
        {
            fvVectorMatrix approxMomJ
            (
                fvm::laplacian(impKf_, D, "laplacian(DD,D)")
              - rho()*fvm::d2dt2(D)
            );
            approxMomJ.relax();
            rAUf() = -1.0/(fvc::interpolate(approxMomJ.A())*one);
        }

        // Calculate pressure equation residual
        scalarField pressureResidual
        (
          - p*rKappa_
          + pressureStabilisation().cellScalar(&rAUf(), true)*one
          - tr(gradD())
        );

        // Make residual extensive
        pressureResidual *= mesh.V();

        // Copy the pressureResidual into the f field as the 4th equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1
        );
    }

    // Calculate the traction at the faces. This must be placed after the
    // pressure solve so that sigma() already carries the pressure component
    // when solvePressure() is active.
    if (!highOrderResidual())
    {
        traction = (n & fvc::interpolate(sigma()));
    }

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the
    // stabilisation is set to zero on traction boundaries
    momentumStabilisation().updateVector(D, &gradD());
    traction += impKf_*momentumStabilisation().faceVector();

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    vectorField residual
    (
        fvc::div(mesh.magSf()*traction)
      + rho()
       *(
            g() - dampingCoeff()*fvc::ddt(D)
        )
    );

    if (highOrderResidual())
    {
#ifndef FOAMEXTEND
        residual -= rho()*hofvc::d2dt2(D);
#endif
    }
    else
    {
        residual -= rho()*fvc::d2dt2(D);
    }

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

#ifdef OPENFOAM_COM
    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();
#endif

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    return 0;
}


label linGeomTotalDispSolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Copy x into the D field
    volVectorField& D = const_cast<volVectorField&>(this->D());
    vectorField& DI = D;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        DI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    D.correctBoundaryConditions();

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();

        {
            // Dimensional consistency factor
            const dimensionedScalar one
            (
                "one", dimensionSet(-2, 4, 4, 0, 0, 0, 0), 1.0
            );

            // Compute the positive face-interpolated reciprocal of the approximate
            // momentum equation diagonal (solid analogue of rAUf), [Pa]
            {
                fvVectorMatrix approxMomJ
                (
                    fvm::laplacian(impKf_, D, "laplacian(DD,D)")
                  - rho()*fvm::d2dt2(D)
                );
                approxMomJ.relax();
                rAUf() = -1.0/(fvc::interpolate(approxMomJ.A())*one);
            }

            fvScalarMatrix approxPressureJ
            (
              - fvm::Sp(rKappa_, p)
              + one*pressureStabilisation().scalarJacobian(p, &rAUf())
            );

            // Insert the pressure equation
            foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
            (
                approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
            );
        }

        // Insert D-in-p equation coeffs coming from tr(grad(D)) == div(D)
        foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
        (
            p,
            D,
            jac,
            blockSize_ - 1,            // row offset
            0,                         // column offset
            solidModel::twoD() ? 2 : 3 // number of scalar components of D
        );

        // Insert p-in-D term
        // Insert "-grad(p)" (equivalent to "-div(p*I)") into the D equation
        foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
        (
            p,
            jac,
            0,                         // row offset
            blockSize_ - 1,            // column offset
            solidModel::twoD() ? 2 : 3 // number of scalar equations to insert
        );
    }

    if (highOrderJacobian())
    {
#ifndef FOAMEXTEND
        // Note: unlike the fallback fvVectorMatrix approxJ path below, we do
        // not currently apply matrix under-relaxation to the high-order
        // Jacobian assembled directly into PETSc. If this becomes important
        // for robustness, an equivalent relaxation step may need to be added.
        tmp<volScalarField> tK = mechanical().bulkModulus();
        const volScalarField& K = tK();

        tmp<volScalarField> tMu = (impK_ - K)*(3.0/4.0);
        const volScalarField& mu = tMu();

        tmp<volScalarField> tLambda = impK_ - 2.0*mu;
        const volScalarField& lambda = tLambda();

        const leastSquaresScheme& reconstruction =
            displacementLeastSquares();

        hofvm::laplacianIntoPETScMatrix
        (
            jac,
            *this,
            reconstruction,
            D,
            mu
        );

        hofvm::laplacianTransposeIntoPETScMatrix
        (
            jac,
            *this,
            reconstruction,
            D,
            mu
        );

        hofvm::laplacianTraceIntoPETScMatrix
        (
            jac,
            *this,
            reconstruction,
            D,
            lambda
        );

        if (momentumStabilisation().supportsHighOrderResidual())
        {
            hofvm::insertAlphaStabIntoPETScMatrix
            (
                jac,
                *this,
                reconstruction,
                D,
                impKf_,
                momentumStabilisation().scaleFactor()
            );
        }

        fvVectorMatrix transientJ
        (
          - rho()*hofvm::d2dt2(D)
        );

        if (dampingCoeff().value() > SMALL)
        {
            transientJ -= dampingCoeff()*rho()*fvmDdtVectorCompat(D);
        }

        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
        (
            transientJ, jac, 0, 0, solidModel::twoD() ? 2 : 3
        );
#endif
    }
    else
    {
        // Calculate a segregated approximation of the Jacobian
        fvVectorMatrix approxJ
        (
            momentumStabilisation().vectorJacobian(D, &impKf_)
          - rho()*fvm::d2dt2(D)
        );

        if (dampingCoeff().value() > SMALL)
        {
            approxJ -= dampingCoeff()*rho()*fvm::ddt(D);
        }

        // Optional: under-relaxation of the linear system
        approxJ.relax();

        // Convert fvMatrix matrix to PETSc matrix
        foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix
        (
            approxJ, jac, 0, 0, solidModel::twoD() ? 2 : 3
        );
    }

    return 0;
}


label linGeomTotalDispSolid::precondition
(
    Vec y,         // Output: y = M^{-1} x
    const Vec x    // Input: vector to be preconditioned
)
{
#ifdef OPENFOAM_COM
    // Take references
    const fvMesh& mesh = this->mesh();
    volVectorField& D = this->D();

    // Extract algebraic RHS (per-cell vector)
    vectorField rhs(mesh.nCells(), vector::zero);
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        rhs,
        0, // first component location
        solidModel::twoD() ? makeList<label>({0,1})
                           : makeList<label>({0,1,2})
    );

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    // Build scalar Laplacian for this component
    fvVectorMatrix DEqn(fvm::laplacian(impKf_, D, "preconditionD"));

    // Overwrite the source
    DEqn.source() = rhs;

    // Solve
    DEqn.solve("preconditionD");

    // Write back to PETSc y
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        D,
        y,
        0, // first component location
        solidModel::twoD() ? makeList<label>({0,1})
                           : makeList<label>({0,1,2})
    );

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

#else // Not OPENFOAM_COM
    FatalErrorInFunction
        << "precondition(...) not implemented for this version of OpenFOAM"
        << exit(FatalError);
#endif

    return 0;
}

#endif // USE_PETSC


tmp<vectorField> linGeomTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
