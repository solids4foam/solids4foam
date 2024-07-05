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

#include "linGeomVelocitySolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "momentumStabilisation.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomVelocitySolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomVelocitySolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void linGeomVelocitySolid::updateTimeIntergrationFields()
{
    tmp<fvVectorMatrix> ddtD(fvm::ddt(D()));
#ifdef OPENFOAM_NOT_EXTEND
    dtUCoeff_.internalFieldRef() 
#else
    dtUCoeff_.internalField() 
#endif                    
                        = ddtD().A();
    dtUCoeff_.correctBoundaryConditions();
    ddtD.clear();

    rdtUCoeff_ == 1.0/dtUCoeff_;

    volVectorField U_Star(fvc::ddt(D()));

#ifdef OPENFOAM_NOT_EXTEND
    dtUSource_.internalFieldRef() 
#else
    dtUSource_.internalField()
#endif 
        = dtUCoeff_.internalField() * D().internalField() 
                - U_Star.internalField();

    forAll(D().boundaryField(),patchI)
    {
#ifdef OPENFOAM_NOT_EXTEND
        dtUSource_.boundaryFieldRef()[patchI] 
#else
        dtUSource_.boundaryField()[patchI] 
#endif  
            = dtUCoeff_.boundaryField()[patchI] * D().boundaryField()[patchI] 
                - U_Star.boundaryField()[patchI];
    }
}

void linGeomVelocitySolid::timeIntegrateU()
{
    D() == rdtUCoeff_ * (U() + dtUSource_);
}

void linGeomVelocitySolid::predict()
{
    Info<< "Applying linear predictor to D" << endl;

    // Predict D using previous time steps
    timeIntegrateU();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomVelocitySolid::linGeomVelocitySolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    dtUCoeff_
    (
        IOobject
        (
            "timeIntegrateUCoeffs",
            U().instance(),
            U().db()
        ),
        U().mesh(),
        dimless/dimTime,
        extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    rdtUCoeff_
    (
        IOobject
        (
            "rTimeIntegrateUCoeffs",
            U().instance(),
            U().db()
        ),
        U().mesh(),
        dimTime,
        extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    dtUSource_
    (
        IOobject
        (
            "timeIntegrateUSource",
            U().instance(),
            U().db()
        ),
        U().mesh(),
        dimVelocity,
        extrapolatedCalculatedFvPatchScalarField::typeName
    )
{
    if(!U().headerOk())
    {
        FatalError("U must be defined in 0 directory!");
    }

    // Force all required old-time fields to be created
    fvm::ddt(U());

    // For consistent restarts, we will calculate the gradient field
    U().correctBoundaryConditions();
    U().storePrevIter();
    D().correctBoundaryConditions();
    D().storePrevIter();

    mechanical().grad(D(), gradD());

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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomVelocitySolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Updating the integration fields at start of timestep
    updateTimeIntergrationFields();

    // Store for efficiency
    volScalarField impKByA_ = impK_*rdtUCoeff_;
    surfaceScalarField impKbyAf_ = impKf_ * fvc::interpolate(rdtUCoeff_);


    if (predictor_)
    {
        predict();
    }

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        SolverPerformance<vector> solverPerfD;
        SolverPerformance<vector>::debug = 0;
#else
        lduSolverPerformance solverPerfD;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the momentum equation for D" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            U().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix UEqn
            (
                rho()*fvm::ddt(U())
             == fvm::laplacian(impKbyAf_, U(), "laplacian(DD,U)")
              - fvc::laplacian(impKbyAf_, U(), "laplacian(DD,U)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(D(), gradD(), impKByA_)
            );

            // Add damping
            if (dampingCoeff().value() > SMALL)
            {
                UEqn += fvm::Sp(dampingCoeff()*rho(),U());
            }

            // Under-relaxation the linear system
            UEqn.relax();

            // Enforce any cell displacements
            //! Not sure about this!
            solidModel::setCellDisps(UEqn);

            // Solve the linear system
            solverPerfD = UEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = rdtUCoeff_ * U();
            
            // Time-integrate U to D
            timeIntegrateU();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Update the momentum equation inverse diagonal field
            // This may be used by the mechanical law when calculating the
            // hydrostatic pressure
            // ! Still applicable?
            const volScalarField DEqnA("DEqnA", UEqn.A());

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());

            // Update impKf to improve convergence
            // Note: impK and rImpK are not updated as they are used for
            // traction boundaries
            //if (iCorr % 10 == 0)
            //{
            //    impKf_ = mechanical().impKf();
            //}
        }
        while
        (
            !converged
            (
                iCorr,
#ifdef OPENFOAM_NOT_EXTEND
                mag(solverPerfD.initialResidual()),
                cmptMax(solverPerfD.nIterations()),
#else
                solverPerfD.initialResidual(),
                solverPerfD.nIterations(),
#endif
                D()
            )
         && ++iCorr < nCorr()
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
    while (mesh().update());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> linGeomVelocitySolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch time integration coeffs
    const scalarField& dtUCoeff = dtUCoeff_.boundaryField()[patchID];

    // Patch reciprocal time integration coeffs
    const scalarField& rdtUCoeff = rdtUCoeff_.boundaryField()[patchID];

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
              - (n & (pSigma - rdtUCoeff*impK*pGradD))
            )*dtUCoeff*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
