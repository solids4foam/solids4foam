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

#include "nonLinGeomTotalLagVelocitySolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "compatibilityFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagVelocitySolid, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagVelocitySolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void nonLinGeomTotalLagVelocitySolid::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    // D() = D().oldTime() + U()*runTime().deltaT();
    D() = D().oldTime() + U()*runTime().deltaT()
        + 0.5*sqr(runTime().deltaT())*A_;

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        // Predict p using the dp/dt field
        p() = p().oldTime() + autoPtrRef(dpdtPtr_)*runTime().deltaT();
        // p() = p().oldTime() + dpdt*runTime().deltaT()
        //     + 0.5*sqr(runTime().deltaT())*d2pdt2;

        sigma() = dev(sigma()) - p()*I;
    }
}


void nonLinGeomTotalLagVelocitySolid::enforceTractionBoundaries
(
    surfaceVectorField& force,
    const volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
#ifdef OPENFOAM_NOT_EXTEND
        vectorField& forceP = force.boundaryFieldRef()[patchI];
#else
        vectorField& forceP = force.boundaryField()[patchI];
#endif

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

            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // traction.boundaryFieldRef()[patchI] =
            //     tracPatch.traction() - nPatch*tracPatch.pressure();
            if (tracPatch.useUndeformedArea())
            {
                const scalarField& magSfPatch =
                    D.mesh().boundary()[patchI].magSf();

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfPatch;
            }
            else
            {
                const scalarField& magSfCurrentPatch =
                    magSfCurrent.boundaryField()[patchI];

                forceP =
                    (
                        tracPatch.traction() - nPatch*tracPatch.pressure()
                    )*magSfCurrentPatch;
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
         || isA<slipFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const vectorField& nPatch = nCurrent.boundaryField()[patchI];

            // Set shear traction to zero
            // traction.boundaryFieldRef()[patchI] =
                // sqr(nPatch) & traction.boundaryField()[patchI];
            forceP = sqr(nPatch) & force.boundaryField()[patchI];
        }
    }
}


bool nonLinGeomTotalLagVelocitySolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for U using PETSc SNES" << endl;

    // If needed, update the current time fields from the future time
    // predictions
    if (curTimeIndex_ != runTime().timeIndex())
    {
        curTimeIndex_ = runTime().timeIndex();

        // Note that forceCurrentTime_.oldTime() is already correctly updated by
        // the OpenFOAM time object

        // Extrapolated forces at the faces
        forceCurrentTime_ =
            surfaceVectorField(forceCurrentTime_.name(), forceFutureTime_);

        // Extrapolated cell displacments
        D() = volVectorField(D().name(), DFutureTime_);
    }

    // Update D and U boundary conditions
    D().correctBoundaryConditions();
    U().correctBoundaryConditions();

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the U field
    vectorField& UI = U();
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        foamPetscSnesHelper::solution(),
        UI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    U().correctBoundaryConditions();

    if (solvePressure())
    {
        notImplemented("Not implemented for solvePressure");

        // // Map the PETSc solution to the p field
        // // p is located in the last ("blockSize - 1") component
        // scalarField& pI = p();
        // foamPetscSnesHelper::ExtractFieldComponents<scalar>
        // (
        //     foamPetscSnesHelper::solution(),
        //     pI,
        //     blockSize_ - 1 // Location of p component
        // );

        // p().correctBoundaryConditions();

        // // Update dpdt
        // autoPtrRef(dpdtPtr_) = fvc::ddt(p());
    }

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();

    // Acceleration
    A_ = fvc::ddt(U());

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


void nonLinGeomTotalLagVelocitySolid::makePDiffusivity() const
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


const surfaceScalarField& nonLinGeomTotalLagVelocitySolid::pDiffusivity() const
{
    if (pDiffusivityPtr_.empty())
    {
        makePDiffusivity();
    }

    return autoPtrRef(pDiffusivityPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagVelocitySolid::nonLinGeomTotalLagVelocitySolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        "U",
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
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    A_
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::d2dt2(D())
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
    dpdtPtr_(),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    ),
    DFutureTime_
    (
        IOobject
        (
            "DFutureTime",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    forceCurrentTime_
    (
        IOobject
        (
            "forceCurrentTime",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce, vector::zero)
    ),
    forceFutureTime_
    (
        IOobject
        (
            "forceFutureTime",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce, vector::zero)
    ),
    curTimeIndex_(-1)
{
    DisRequired();

    // Ensure U is created (and possibly read)
    U();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // Store the old time stress
    sigma().storeOldTime();
    forceCurrentTime_.storeOldTime();

    // It is important to call the stress calculation procedure during the
    // constructor to allow it to correctly initialise fields
    mechanical().correct(sigma());

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

        // Ensure p is created and the oldTime is stored
        p().oldTime();

        // Initialise dpdt field
        dpdtPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                fvc::ddt(p())
            )
        );
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

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        DD() = D() - D().oldTime();
        mechanical().grad(D(), gradD());
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
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

    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
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
                << " solution algorithm" << abort(FatalError);
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


bool nonLinGeomTotalLagVelocitySolid::evolve()
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
    // else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    // {
    //     return evolveImplicitSegregated();
    // }
    // else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    // {
    //     return evolveExplicit();
    // }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            // << solutionAlgorithmNames_.names() << endl;
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            // << solidModel::solutionAlgorithmNames_
            //    [
            //        solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
            //    ]
            // << solidModel::solutionAlgorithmNames_
            //    [
            //        solidModel::solutionAlgorithm::EXPLICIT
            //    ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


#ifdef USE_PETSC

label nonLinGeomTotalLagVelocitySolid::initialiseJacobian(Mat& jac)
{
    // Initialise based on compact stencil fvMesh
    return foamPetscSnesHelper::initialiseJacobian(jac, mesh(), blockSize_);
}


label nonLinGeomTotalLagVelocitySolid::initialiseSolution(Vec& x)
{
    // Initialise based on mesh.nCells()
    return foamPetscSnesHelper::initialiseSolution(x, mesh(), blockSize_);
}


label nonLinGeomTotalLagVelocitySolid::formResidual
(
    Vec f,
    const Vec x
)
{
    const fvMesh& mesh = this->mesh();
    const dimensionedScalar& deltaT = runTime().deltaT();

    // Copy x into the U field
    volVectorField& U = const_cast<volVectorField&>(this->U());
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    U.correctBoundaryConditions();

    // In the momentum equation, we will approximate div(sigma) as
    //     div(sigma) = 0.5*div(sigma[n-1] + sigma[n+1])
    // where n is the current time, n-1 is the old time, and n+1 is the future
    // time, so sigma[n+1] is calculated using D[n+1], i.e. D at the future time
    // step (DFutureTime_)
    // We approximate D[n+1] using a second-order Adams–Bashforth extrapolation
    // from the current time using the velocity. This is the approach given in
    // Jaiman and Joshi, 2022, in Eq. 6.16.

    // Approximate the displacement at the future time step
    DFutureTime_ = D() + 0.5*deltaT*(3*U - U.oldTime());

    // Calculate future-time gradient of displacement
    mechanical().grad(DFutureTime_, gradD());

    // Total deformation gradient at the future time
    F_ = I + gradD().T();

    // Inverse of the deformation gradient at the future time
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient at the future time
    J_ = det(F_);

    // Calculate the stress at the future time using run-time selectable
    // mechanical law
    mechanical().correct(sigma());

    if (solvePressure())
    {
        notImplemented("Not checked yet for solvePressure");

        // // Copy x into the p field
        // volScalarField& p = const_cast<volScalarField&>(this->p());
        // scalarField& pI = p;
        // foamPetscSnesHelper::ExtractFieldComponents<scalar>
        // (
        //     x, pI, blockSize_ - 1
        // );

        // // Enforce the boundary conditions
        // p.correctBoundaryConditions();

        // // Replace the pressure component of stress
        // sigma() = dev(sigma()) - p*I;
    }

    // Calculate the deformed normals at the future step
    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfNew
    (
        fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
    );
    const surfaceScalarField magSfNew(mag(SfNew));
    const surfaceVectorField nNew(SfNew/magSfNew);

    // Future-time traction vectors at the faces
    surfaceVectorField tractionNew
    (
        nNew & fvc::interpolate(sigma())
    );

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // Todo: consider stabilisation in the deformed configuration (or formulate
    // the traction in the reference configuration)
    // CAREFUL: D may not have corrected BCs with non-ortho corrections; maybe
    // formulate in terms of alpha scheme instead
    {
        const scalar scaleFactor =
            readScalar(stabilisation().dict().lookup("scaleFactor"));
        const surfaceTensorField gradDf(fvc::interpolate(gradD()));
        tractionNew +=
            scaleFactor*impKf_*(fvc::snGrad(DFutureTime_) - (n & gradDf));
    }

    // Update the future-time force at the faces
    forceFutureTime_ = magSfNew*tractionNew;

    // Enforce traction boundary conditions at the new time
    // TODO: we need to access traction at the new time (if available) ...
    // but how do we extrapolate this consistently ...
    // Warning
    //     << "I need to consistently extrapolate the tractions" << endl;
    enforceTractionBoundaries(forceFutureTime_, D(), nNew, magSfNew);

    // The residual vector is defined as
    // F = 0.5*div(sigma[n-1] + sigma[n+1]) + rho*g
    //     - rho*ddt(U) - dampingCoeff*rho*U + stabilisationTerm
    // where the stabilisationTerm has been rolled into the div(sigma) terms
    vectorField residual
    (
        0.5*fvc::div(forceCurrentTime_.oldTime() + forceFutureTime_)
      + rho()
       *(
            g() - fvc::ddt(U) - dampingCoeff()*U
        )
    );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh.V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    //residual -= fvOptions()(ds_, const_cast<volVectorField&>(U))().source();

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

    if (solvePressure())
    {
        notImplemented("Not checked yet for solvePressure");

        // volScalarField& p = const_cast<volScalarField&>(this->p());

        // // Divided by bulkModulus form
        // const volScalarField kappa("kappa", mechanical().bulkModulus());
        // const surfaceScalarField kappaf(fvc::interpolate(kappa));
        // const dimensionedScalar omega(solidModelDict().lookup("omega"));
        // scalarField pressureResidual
        // (
        //   - p/kappa
        //   + fvc::laplacian
        //     (
        //         omega/sqr(mesh.deltaCoeffs()/impKf_), p, "laplacian(Dp,p)"
        //     )
        //   - 0.5*(pow(J_, 2.0) - 1.0)/J_
        // );

        // // Make residual extensive
        // pressureResidual *= mesh.V();

        // // Copy the pressureResidual into the f field as the final equation
        // foamPetscSnesHelper::InsertFieldComponents<scalar>
        // (
        //     pressureResidual, f, blockSize_ - 1
        // );
    }

    return 0;
}


label nonLinGeomTotalLagVelocitySolid::formJacobian
(
    Mat jac,
    const Vec x
)
{
    // Copy x into the U field
    volVectorField& U = const_cast<volVectorField&>(this->U());
    vectorField& UI = U;
    foamPetscSnesHelper::ExtractFieldComponents<vector>
    (
        x,
        UI,
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

    // Enforce the boundary conditions
    U.correctBoundaryConditions();

    if (solvePressure())
    {
        notImplemented("Not checked yet for solvePressure");

        // // Copy x into the p field
        // volScalarField& p = const_cast<volScalarField&>(this->p());
        // scalarField& pI = p;
        // foamPetscSnesHelper::ExtractFieldComponents<scalar>
        // (
        //     x, pI, blockSize_ - 1
        // );

        // // Enforce the boundary conditions
        // p.correctBoundaryConditions();
    }

    // Is "impKf_*runTime().deltaT()" the best Laplacian coefficient?
    // We calculate div(sigma) as 0.5*div(sigma[n-1] + sigma[n+1]) where
    // sigma[n+1] is a function of
    //     DFutureTime = D + 0.5*deltaT*(3*U - U.oldTime())
    // Since we use the Laplacian(impK,D) as an approximation of div(sigma),
    // then the implicit part (in terms of U) would be 0.5*sigma[n+1] as a
    // function of 0.5*deltaT*3*U, so the implicit Laplacian becomes
    // laplacian(1.5*impK*deltaT, U)

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian(1.5*impKf_*runTime().deltaT(), U, "laplacian(DU,U)")
      - rho()*fvm::ddt(U)
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*U;
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
        notImplemented("Not checked yet for solvePressure");

        // const volScalarField& p = this->p();

        // const volScalarField kappa("kappa", mechanical().bulkModulus());
        // //const volScalarField rKappa(1.0/mechanical().bulkModulus());
        // const volScalarField rKappa(1.0/kappa);
        // const surfaceScalarField kappaf(fvc::interpolate(kappa));
        // const dimensionedScalar omega(solidModelDict().lookup("omega"));
        // {
        //     // Calculate pressure equation matrix
        //     //const dimensionedScalar one("one", dimless, 1);
        //     fvScalarMatrix approxPressureJ
        //     (
        //       - fvm::Sp(rKappa, p)
        //       + fvm::laplacian
        //         (
        //             omega/sqr(mesh().deltaCoeffs())/impKf_,
        //             p,
        //             "jacobian-laplacian(rAU,p)"
        //         )
        //     );

        //     // Insert the pressure equation
        //     foamPetscSnesHelper::InsertFvMatrixIntoPETScMatrix<scalar>
        //     (
        //         approxPressureJ, jac, blockSize_ - 1, blockSize_ - 1, 1
        //     );
        // }

        // // Insert D-in-p equation coeffs coming from tr(grad(D)) == div(D)
        // foamPetscSnesHelper::InsertFvmDivUIntoPETScMatrix
        // (
        //     p,
        //     D,
        //     jac,
        //     blockSize_ - 1,            // row offset
        //     0,                         // column offset
        //     solidModel::twoD() ? 2 : 3 // number of scalar components of D
        // );

        // // Insert p-in-D term
        // // Insert "-grad(p)" (equivalent to "-div(p*I)") into the D equation
        // foamPetscSnesHelper::InsertFvmGradIntoPETScMatrix
        // (
        //     p,
        //     jac,
        //     0,                         // row offset
        //     blockSize_ - 1,            // column offset
        //     solidModel::twoD() ? 2 : 3 // number of scalar equations to insert
        // );
    }

    return 0;
}

#endif // USE_PETSC

tmp<vectorField> nonLinGeomTotalLagVelocitySolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
