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


void linGeomVelocitySolid::predict()
{
    Info<< "Linear predictor using DD" << endl;

    // Predict D using the increment of displacement field from the previous
    // time-step
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement increment
    mechanical().grad(DD(), gradDD());

    // Update gradient of total displacement
    gradD() = gradD().oldTime() + gradDD();

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
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false))
{
    UisRequired();

    // Force all required old-time fields to be created
    fvm::ddt(U());

    // For consistent restarts, we will calculate the gradient field
    U().storePrevIter();
    D().storePrevIter();
    mechanical().grad(D(), gradD());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomVelocitySolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    if (predictor_)
    {
        predict();
    }

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
        SolverPerformance<vector> solverPerfU;
        SolverPerformance<vector>::debug = 0;
#else
        lduSolverPerformance solverPerfU;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the momentum equation for U" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            U().storePrevIter();

            // Linear momentum equation velocity form
            fvVectorMatrix UEqn
            (
                rho()*fvm::ddt(U())
             == fvm::laplacian(impKf_*runTime().deltaT(), U(), "laplacian(DD,U)")
              - fvc::laplacian(impKf_*runTime().deltaT(), U(), "laplacian(DD,U)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(U(), gradU(), impK_*runTime().deltaT())
            );

            // Under-relaxation the linear system
            UEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(UEqn);

            // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
            const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif

            // Solve the linear system
            solverPerfU = UEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(U(), iCorr);
            
            // Update the increment of displacement
            DD() = U()*runTime().deltaT();
            // Update the total displacement
            D() = D().oldTime() + DD();

            // Update gradient of velocity
            mechanical().grad(U(), gradU());

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Calculate the stress using run-time selectable mechanical law
            const volScalarField UEqnA("DEqnA", UEqn.A());
            mechanical().correct(sigma());
        }
        while
        (
            !converged
            (
                iCorr,
#ifdef OPENFOAMESIORFOUNDATION
                mag(solverPerfU.initialResidual()),
                cmptMax(solverPerfU.nIterations()),
#else
                solverPerfU.initialResidual(),
                solverPerfU.nIterations(),
#endif
                U()
            ) && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

    }
    while (mesh().update());

#ifdef OPENFOAMESIORFOUNDATION
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

    // Patch mechanical property
    const scalarField& pImpK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& pRImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradU = gradU().boundaryField()[patchID];

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
              - (n & (pSigma - pImpK*runTime().deltaTValue()*pGradU))
            )*pRImpK/runTime().deltaTValue()
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
