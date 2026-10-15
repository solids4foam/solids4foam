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
#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "compatibilityFunctions.H"
#include "mechanicalConstitutiveLawManager.H"


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
    gradD() = fvc::grad(D());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    correctStress();

    if (solvePressure())
    {
        // Predict p using the dp/dt field
        p() = p().oldTime() + autoPtrRef(dpdtPtr_)*runTime().deltaT();
        // p() = p().oldTime() + dpdt*runTime().deltaT()
        //     + 0.5*sqr(runTime().deltaT())*d2pdt2;

        sigma() = dev(sigma()) - p()*I;
    }
}


void nonLinGeomTotalLagVelocitySolid::correctStress()
{
    if (linearGeometryFormulation_)
    {
        mechanicalManager().updateStressSmallStrain
        (
            gradD(),
            gradD().oldTime(),
            mesh().time().deltaTValue(),
            sigma()
        );
    }
    else
    {
        mechanicalManager().updateStressFiniteStrain
        (
            F_,
            F_.oldTime(),
            Finv_,
            Finv_.oldTime(),
            J_,
            J_.oldTime(),
            mesh().time().deltaTValue(),
            sigma()
        );
    }
}


tmp<volScalarField> nonLinGeomTotalLagVelocitySolid::makeImpK() const
{
    // With the mixed displacement-pressure formulation the implicit stiffness
    // is the deviatoric scalar tangent, as in nonLinGeomTotalLagTotalDispSolid
    const tangentRequest req =
        solvePressure()
      ? tangentRequest::scalarDeviatoric
      : tangentRequest::scalar;

    if (linearGeometryFormulation_)
    {
        return lawImpK(mechanicalManager(), req);
    }

    tmp<volScalarField> tImpK
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimForce/dimArea, 0),
            calculatedFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    volScalarField& impK = tImpK.ref();
#else
    volScalarField& impK = tImpK();
#endif

    // A tangent query at the current deformation gradient, which neither
    // writes a stress nor disturbs history; F_, Finv_ and J_ are declared
    // before impK_, so they are already built
    mechanicalManager().updateScalarTangentFiniteStrain
    (
        F_,
        F_.oldTime(),
        Finv_,
        Finv_.oldTime(),
        J_,
        J_.oldTime(),
        mesh().time().deltaTValue(),
        impK,
        req
    );

    return tImpK;
}


void nonLinGeomTotalLagVelocitySolid::makeVelocityMomentumStabilisationCoeff()
{
    if (!velocityMomentumStabilisationCoeffPtr_.empty())
    {
        FatalErrorInFunction
            << "pointer already set!" << abort(FatalError);
    }

    surfaceScalarField cSqrtRho(sqrt(impKf_*fvc::interpolate(rho())));
#ifdef OPENFOAM_COM
    cSqrtRho.setOriented(false);
#endif

    velocityMomentumStabilisationCoeffPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "velocityMomentumStabilisationCoeff",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            cSqrtRho/mesh().deltaCoeffs()
        )
    );
#ifdef OPENFOAM_COM
    velocityMomentumStabilisationCoeffPtr_().setOriented(false);
#endif
}


const surfaceScalarField&
nonLinGeomTotalLagVelocitySolid::velocityMomentumStabilisationCoeff()
{
    if (velocityMomentumStabilisationCoeffPtr_.empty())
    {
        makeVelocityMomentumStabilisationCoeff();
    }

    return velocityMomentumStabilisationCoeffPtr_();
}


void nonLinGeomTotalLagVelocitySolid::enforceTractionBoundaries
(
    surfaceVectorField& force,
    volVectorField& D,
    const surfaceVectorField& nCurrent,
    const surfaceScalarField& magSfCurrent
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
        vectorField& forceP = boundaryFieldRef(force)[patchI];

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

            // For now, disallow use of the "pressure" field
            if (gMax(mag(tracPatch.pressure())) > SMALL)
            {
                FatalErrorInFunction
                    << "The 'pressure' on traction patches is not yet "
                    << "implemented; it must be zero" << exit(FatalError);
            }

            // The current traction (set by the boundary condition or the FSI
            // coupling) is applied at the new time, as is the stress
            const vectorField& trac = tracPatch.traction();

            if (tracPatch.useUndeformedArea())
            {
                forceP = trac*D.mesh().boundary()[patchI].magSf();
            }
            else
            {
                forceP = trac*magSfCurrent.boundaryField()[patchI];
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


            /// Hmnnnn.. do we need to extrapolate the traction on all
            // boundaries!?
        }
    }
}


scalar nonLinGeomTotalLagVelocitySolid::backwardDeltaT0() const
{
    // Follow the backward ddt scheme of each fork in deciding the first
    // time-step
#ifdef OPENFOAM_COM
    if (runTime().timeIndex() < 2)
#else
    if (U().oldTime().timeIndex() == U().oldTime().oldTime().timeIndex())
#endif
    {
        return GREAT;
    }

    return runTime().deltaT0Value();
}


dimensionedScalar
nonLinGeomTotalLagVelocitySolid::DImplicitVelocityCoeff() const
{
    // Coefficient of the new time in the backward scheme
    const scalar deltaT = runTime().deltaTValue();
    const scalar deltaT0 = backwardDeltaT0();
    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);

    return runTime().deltaT()/coefft;
}


void nonLinGeomTotalLagVelocitySolid::updateDImplicit
(
    const volVectorField& U
)
{
    // The backward scheme applied to the displacement,
    //     coefft*D - coefft0*D[n-1] + coefft00*D[n-2] = deltaT*U
    // so that DImplicit_ and ddt(U) are consistent
    const scalar deltaT = runTime().deltaTValue();
    const scalar deltaT0 = backwardDeltaT0();
    const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0 = coefft + coefft00;

    DImplicit_ =
        (
            runTime().deltaT()*U
          + coefft0*D().oldTime()
          - coefft00*D().oldTime().oldTime()
        )/coefft;

    // Impose the prescribed displacements, as D holds them at the new time
    forAll(useBoundaryFaceValuesDImplicit_, patchI)
    {
        if (useBoundaryFaceValuesDImplicit_[patchI])
        {
            boundaryFieldRef(DImplicit_)[patchI] =
                D().boundaryField()[patchI];
        }
    }
}


bool nonLinGeomTotalLagVelocitySolid::evolveSnes()
{
#ifdef USE_PETSC
    Info<< "Solving the momentum equation for U using PETSc SNES" << endl;

    // Update D and U boundary conditions
    D().correctBoundaryConditions();
    U().correctBoundaryConditions();

    // Seed the SNES solution vector from the current velocity on every time
    // step, so that the initial velocity is not overwritten by a stale
    // solution vector
    foamPetscSnesHelper::InsertFieldComponents<vector>
    (
#ifdef OPENFOAM_NOT_EXTEND
        U().primitiveFieldRef(),
#else
        U().internalField(),
#endif
        foamPetscSnesHelper::solution(),
        0, // Location of first component
        solidModel::twoD()
      ? makeList<label>({0,1})
      : makeList<label>({0,1,2})
    );

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
    else
    {
        // Update p field, which we use
    }

    // Update the displacement from the converged velocity
    updateDisplacement();

#else

    FatalErrorInFunction
        << "To use PETSc with solids4foam, set the PETSC_DIR to point to your "
        << "PETSC installation directory and re-build solids4foam"
        << exit(FatalError);

#endif

    return true;
}


word nonLinGeomTotalLagVelocitySolid::modelTypeName()
{
    return typeName;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagVelocitySolid::nonLinGeomTotalLagVelocitySolid
(
    Time& runTime,
    const word& region
)
:
    nonLinGeomTotalLagVelocitySolid(modelTypeName(), runTime, region, false)
{
}


nonLinGeomTotalLagVelocitySolid::nonLinGeomTotalLagVelocitySolid
(
    const word& modelType,
    Time& runTime,
    const word& region,
    const bool linearGeometryFormulation
)
:
    solidModel(modelType, runTime, region),
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
    linearGeometryFormulation_(linearGeometryFormulation),
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
    impK_(makeImpK()),
    impKf_(fvc::interpolate(impK_)),
    rImpK_(1.0/impK_),
    dpdtPtr_(),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    velocityBasedMomentumStabilisation_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "velocityBasedMomentumStabilisation",
            false
        )
    ),
    blockSize_
    (
        solvePressure()
      ? label(solidModel::twoD() ? 3 : 4)
      : label(solidModel::twoD() ? 2 : 3)
    ),
    velocityMomentumStabilisationCoeffPtr_(),
    DImplicit_
    (
        IOobject
        (
            "DImplicit",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    useBoundaryFaceValuesDImplicit_
    (
        IOobject
        (
            "useBoundaryFaceValues_DImplicit",
            runTime.constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        boolList(mesh().boundary().size(), false)
    ),
    useBoundaryFaceValuesU_
    (
        IOobject
        (
            "useBoundaryFaceValues_U",
            runTime.constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        boolList(mesh().boundary().size(), false)
    )
{
    Info<< "linearGeometryFormulation = " << linearGeometryFormulation_ << nl
        << "predictor = " << predictor_ << nl
        << "velocityBasedMomentumStabilisation = "
        << velocityBasedMomentumStabilisation_ << nl
        << "solvePressure = " << solvePressure() << endl;

    // The momentum stabilisation model dimensions depend on which field is
    // used to update it: DImplicit_ has dimensions [L], while U has [L/T].
    if (velocityBasedMomentumStabilisation_)
    {
        resetMomentumStabilisation(dimless/dimTime);
    }
    else
    {
        resetMomentumStabilisation(dimless);
    }

    DisRequired();

    // Ensure U is created (and possibly read)
    U();

    // Use boundary face values only on fixed-displacement patches
    forAll(useBoundaryFaceValuesDImplicit_, patchI)
    {
        if
        (
            isA<fixedDisplacementFvPatchVectorField>
            (
                D().boundaryField()[patchI]
            )
        )
        {
            useBoundaryFaceValuesDImplicit_[patchI] = true;
            useBoundaryFaceValuesU_[patchI] = true;
        }
    }

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // Store the old time stress
    sigma().storeOldTime();

    // It is important to call the stress calculation procedure during the
    // constructor to allow it to correctly initialise fields
    correctStress();

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
        gradD() = fvc::grad(D());
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);

        gradD().storeOldTime();

        // The framework restores its own constitutive state through the
        // registered state fields, so there is nothing to tell it here
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
                        boundaryFieldRef(D())[patchI]
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
    else
    {
        FatalErrorIn("bool nonLinGeomTotalLagVelocitySolid::evolve()")
            << "The " << type() << " solid model only supports the "
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << " solution algorithm; the selected algorithm is "
            << solidModel::solutionAlgorithmNames_[solutionAlg()] << nl
            << "Add 'solutionAlgorithm "
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << ";' to the " << type() << "Coeffs dictionary."
            << abort(FatalError);
    }

    // Keep compiler happy
    return true;
}


void nonLinGeomTotalLagVelocitySolid::updateDisplacement()
{
    // Displacement at the new time from the converged velocity
    updateDImplicit(U());
    D() = volVectorField(D().name(), DImplicit_);
    D().correctBoundaryConditions();

    // Update the gradient of displacement and the stress
    gradD() = fvc::grad(D());

    if (!linearGeometryFormulation_)
    {
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);
    }

    correctStress();

    // Interpolate the cell displacements to the vertices
    interpolatePointDisplacement(D(), gradD(), pointD());
    pointD().correctBoundaryConditions();
    pointDD() = pointD() - pointD().oldTime();

    // Acceleration
    A_ = fvc::ddt(U());
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
    if (solvePressure())
    {
        notImplemented("Not checked yet for solvePressure");
    }

    const fvMesh& mesh = this->mesh();

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
    D().correctBoundaryConditions();

    // The stress is evaluated at the new time, with the displacement given by
    // the backward scheme applied to the velocity, consistent with ddt(U)
    updateDImplicit(U);

    // Calculate the gradient of displacement at the new time
    gradD() = fvc::grad(DImplicit_);

    if (linearGeometryFormulation_)
    {
        // Small strain stress at the new time
        correctStress();
    }
    else
    {
        // Total deformation gradient at the new time
        F_ = I + gradD().T();

        // Inverse of the deformation gradient at the new time
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient at the new time
        J_ = det(F_);

        // Calculate the stress at the new time using run-time selectable
        // mechanical law
        correctStress();

        if (solvePressure())
        {
            notImplemented("Not checked yet for solvePressure");
        }
    }

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Force at the faces at the new time
    autoPtr<surfaceVectorField> forcePtr;

    if (linearGeometryFormulation_)
    {
        // Traction vectors at the faces
        surfaceVectorField traction(n & fvc::interpolate(sigma()));

        // Add stabilisation to the traction
        // We add this before enforcing the traction condition as the stabilisation
        // is set to zero on traction boundaries
        if (velocityBasedMomentumStabilisation_)
        {
            const volTensorField gradU("gradU", fvc::grad(U));
            momentumStabilisation().updateVector(U, &gradU);
            traction +=
                velocityMomentumStabilisationCoeff()
               *momentumStabilisation().faceVector();
        }
        else
        {
            momentumStabilisation().updateVector(DImplicit_, &gradD());
            traction += impKf_*momentumStabilisation().faceVector();
        }

        forcePtr.reset(new surfaceVectorField("force", mesh.magSf()*traction));

        // Enforce traction boundary conditions
        enforceTractionBoundaries(forcePtr(), D(), n, mesh.magSf());
    }
    else
    {
        // Calculate the deformed normals at the new time
        const surfaceVectorField SfNew
        (
            fvc::interpolate(J_*Finv_.T()) & mesh.Sf()
        );
        const surfaceScalarField magSfNew(mag(SfNew));
        const surfaceVectorField nNew(SfNew/magSfNew);

        // Traction vectors at the faces
        surfaceVectorField tractionNew(nNew & fvc::interpolate(sigma()));

        // Add stabilisation to the traction
        // We add this before enforcing the traction condition as the stabilisation
        // is set to zero on traction boundaries
        // Todo: consider stabilisation in the deformed configuration (or formulate
        // the traction in the reference configuration)
        // CAREFUL: D may not have corrected BCs with non-ortho corrections; maybe
        // formulate in terms of alpha scheme instead
        if (velocityBasedMomentumStabilisation_)
        {
            const volTensorField gradU("gradU", fvc::grad(U));
            momentumStabilisation().updateVector(U, &gradU);
            tractionNew +=
                velocityMomentumStabilisationCoeff()
               *momentumStabilisation().faceVector();
        }
        else
        {
            momentumStabilisation().updateVector(DImplicit_, &gradD());
            tractionNew += impKf_*momentumStabilisation().faceVector();
        }

        forcePtr.reset(new surfaceVectorField("force", magSfNew*tractionNew));

        // Enforce traction boundary conditions
        enforceTractionBoundaries(forcePtr(), D(), nNew, magSfNew);
    }

    // The residual vector is defined as
    // F = div(sigma) + rho*g - rho*ddt(U) - dampingCoeff*rho*U
    // where the stabilisation term has been rolled into the div(sigma) term
    vectorField residual
    (
        fvc::div(forcePtr())
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

    // The stress is a function of DImplicit_, whose derivative with respect
    // to U is DImplicitVelocityCoeff() (2*deltaT/3 for a constant
    // time-step). Since we use laplacian(impK, D) as an approximation of
    // div(sigma), the implicit Laplacian becomes
    // laplacian(impK*DImplicitVelocityCoeff(), U)

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian
        (
            impKf_*DImplicitVelocityCoeff(), U, "laplacian(DU,U)"
        )
      - rho()*fvm::ddt(U)
    );

    if (dampingCoeff().value() > SMALL)
    {
        // Implicit, so that the damping enters the Jacobian coefficients
        approxJ -= fvm::Sp(dampingCoeff()*rho(), U);
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

        // const volScalarField kappa("kappa", mechanicalManager().kappa());
        // //const volScalarField rKappa(1.0/mechanicalManager().kappa());
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
