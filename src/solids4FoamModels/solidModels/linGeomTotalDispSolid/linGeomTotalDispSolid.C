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

            traction.boundaryFieldRef()[patchI] =
                tracPatch.traction() - nPatch*tracPatch.pressure();
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
            traction.boundaryFieldRef()[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
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
            // We add this before enforcing the traction condition as the stabilisation
            // is set to zero on traction boundaries
            // To-do: add a stabilisation traction function to momentumStabilisation
            const scalar scaleFactor =
                readScalar(stabilisation().dict().lookup("scaleFactor"));
            const surfaceTensorField gradDf(fvc::interpolate(gradD()));
            traction += scaleFactor*impKf_*(fvc::snGrad(D()) - (n & gradDf));

            // Enforce traction boundary conditions
            enforceTractionBoundaries(traction, D(), n);

            // Linear momentum equation total displacement form
            fvVectorMatrix DEqn
            (
                rho()*fvm::d2dt2(D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(mesh().magSf()*traction)
              + rho()*g()
              + fvOptions()(ds_, D())
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
                            D().primitiveField()
                          - D().prevIter().primitiveField()
                        )
                    )
                );

            // Norm of the solution
            xNorm = sqrt(gSum(magSqr(D().primitiveField())));

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
            D().primitiveFieldRef(),
            foamPetscSnesHelper::solution(),
            0, // Location of first component
            solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
        );
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the D field
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        D().primitiveFieldRef(),
        0, // Location of first component
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    D().correctBoundaryConditions();

    if (solvePressure())
    {
        // Map the PETSc solution to the p field
        // p is located in the 4th component
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            foamPetscSnesHelper::solution(),
            p().primitiveFieldRef(),
            blockSize_ - 1 // Location of p component
        );

        p().correctBoundaryConditions();
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

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
                D.primitiveFieldRef().replace(dirI, 0.0);
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
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // Solve the momentum equation for acceleration
    A_ = fvc::div(mesh.magSf()*traction)/rho
       + g()
       - dampingCoeff()*fvc::ddt(D);

    return true;
}

void linGeomTotalDispSolid::makePDiffusivity() const
{
    if (pDiffusivityPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!" << abort(FatalError);
    }

    const scalar pressureSmoothingCoeff
    (
        readScalar(solidModelDict().lookup("pressureSmoothingCoeff"))
    );

    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D())
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(D());
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    pDiffusivityPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "pDiffusivity",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -pressureSmoothingCoeff*impKf_/fvc::interpolate(approxJ.A())
        )
    );
}


const surfaceScalarField& linGeomTotalDispSolid::pDiffusivity() const
{
    if (pDiffusivityPtr_.empty())
    {
        makePDiffusivity();
    }

    return pDiffusivityPtr_.ref();
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
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        mesh().nCells(),
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3),
        solvePressure()
      ? (
            solidModel::twoD()
          ? labelListList({ {0, 1}, {2} })    // 0-1: momentum. 2: pressure
          : labelListList({ {0, 1, 2}, {3} }) // 0-2: momentum. 3: pressure
        )
      : labelListList({}), // No need to assign fields
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
    pDiffusivityPtr_(),
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
        mesh().gradScheme("grad(" + D().name() +')')
    );

    if
    (
        solutionAlg() == solutionAlgorithm::PETSC_SNES
     || solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED
    )
    {
        if (gradDScheme != "leastSquaresS4f")
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " and "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithms" << abort(FatalError);
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
                        D().boundaryFieldRef()[patchI]
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
                DimensionedField<scalar, Foam::surfaceMesh>
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


label linGeomTotalDispSolid::initialiseJacobian(Mat jac)
{
    // Initialise based on compact stencil fvMesh
    return Foam::initialiseJacobian(jac, mesh(), blockSize_);
}


label linGeomTotalDispSolid::formResidual
(
    PetscScalar *f,
    const PetscScalar *x
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
        blockSize_, // Block size of x
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    // Enforce the boundary conditions
    D.correctBoundaryConditions();

    // Update gradient of displacement
    mechanical().grad(D, gradD());

    // Enforce the boundary conditions again for any conditions that use gradD
    D.correctBoundaryConditions();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Copy x into the p field
        volScalarField& p = const_cast<volScalarField&>(this->p());
        scalarField& pI = p;
        foamPetscSnesHelper::ExtractFieldComponents<scalar>
        (
            x, pI, blockSize_ - 1, blockSize_
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();

        // Replace the pressure component of stress
        sigma() = dev(sigma()) - p*I;
    }

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma()));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD()));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

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
            g() - fvc::d2dt2(D) - dampingCoeff()*fvc::ddt(D)
        )
    );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();

    // Copy the residual into the f field
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
        residual,
        f,
        0, // Location of first component
        blockSize_, // Block size of x
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
    );

    if (solvePressure())
    {
        volScalarField& p = const_cast<volScalarField&>(this->p());

        // Calculate pressure equation residual
        // Res = p/k + div(D) - gamma*laplacian(p) + gamma*div(grad(p))
        // where
        //   - k: bulk modulus
        //   - gamma: "pDiffusivity" controls the amount of smoothing

        // scalarField pressureResidual
        // (
        //   - p
        //   + fvc::laplacian(pDiffusivity(), p, "laplacian(Dp,p)")
        //   - fvc::div(pDiffusivity()*mesh.Sf() & fvc::interpolate(fvc::grad(p)))
        //   - mechanical().bulkModulus()*tr(gradD())
        //   //- mechanical().bulkModulus()*fvc::div(D)
        // );

        // Divided by bulkModulus form
        const volScalarField kappa("kappa", mechanical().bulkModulus());
        const surfaceScalarField kappaf(fvc::interpolate(kappa));
        scalarField pressureResidual
        (
          - p/kappa
          + fvc::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
          - fvc::div((pDiffusivity()/kappaf)*mesh.Sf() & fvc::interpolate(fvc::grad(p)))
          - tr(gradD())
        );

        // Make residual extensive
        pressureResidual *= mesh.V();

        // Copy the pressureResidual into the f field as the 4th equation
        foamPetscSnesHelper::InsertFieldComponents<scalar>
        (
            pressureResidual, f, blockSize_ - 1, blockSize_
        );
    }

    return 0;
}


label linGeomTotalDispSolid::formJacobian
(
    Mat jac,
    const PetscScalar *x
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
        blockSize_, // Block size of x
        solidModel::twoD() ? labelList({0,1}) : labelList({0,1,2})
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
            x, pI, blockSize_ - 1, blockSize_
        );

        // Enforce the boundary conditions
        p.correctBoundaryConditions();
    }

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D, "laplacian(DD,D)")
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

    if (solvePressure())
    {
        const volScalarField& p = this->p();

        const volScalarField kappa("kappa", mechanical().bulkModulus());
        //const volScalarField rKappa(1.0/mechanical().bulkModulus());
        const volScalarField rKappa(1.0/kappa);
        const surfaceScalarField kappaf(fvc::interpolate(kappa));
        {
            // Calculate pressure equation matrix
            const dimensionedScalar one("one", dimless, 1);
            // fvScalarMatrix approxPressureJ
            // (
            //   - fvm::Sp(one, p)
            //   + fvm::laplacian(pDiffusivity(), p, "laplacian(Dp,p)")
            // );
            fvScalarMatrix approxPressureJ
            (
              - fvm::Sp(rKappa, p)
              + fvm::laplacian(pDiffusivity()/kappaf, p, "laplacian(Dp,p)")
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

    return 0;
}


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
