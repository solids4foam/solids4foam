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

#include "linGeomSolid.H"
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

defineTypeNameAndDebug(linGeomSolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void linGeomSolid::predict()
{
    Info<< "Linear predictor using DD" << endl;

    // Predict DD using the previous time steps
    DD() = 2.0*DD().oldTime() - DD().oldTime().oldTime();

    // Update the total displacement
    D() = D().oldTime() + DD();

    // Update gradient of displacement increment
    mechanical().grad(DD(), gradDD());

    // Update gradient of total displacement
    gradD() = gradD().oldTime() + gradDD();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomSolid::linGeomSolid
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
    DDisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(DD());
    fvc::d2dt2(D().oldTime());
    rhoD2dt2D().oldTime();

    // For consistent restarts, we will calculate the gradient field
    DD().correctBoundaryConditions();
    DD().storePrevIter();
    mechanical().grad(DD(), gradDD());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomSolid::evolve()
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
        SolverPerformance<vector> solverPerfDD;
        SolverPerformance<vector>::debug = 0;
#else
        lduSolverPerformance solverPerfDD;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the momentum equation for DD" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            DD().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DDEqn
            (
                rho()*fvm::d2dt2(DD())
              + rhoD2dt2D().oldTime()
             == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
              - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(DD(), gradDD(), impK_)
            );

            // Under-relaxation the linear system
            DDEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DDEqn);

            // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
            const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif

            // Solve the linear system
            solverPerfDD = DDEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(DD(), iCorr);

            // Update the total displacement
            D() = D().oldTime() + DD();

            // Update gradient of displacement increment
            mechanical().grad(DD(), gradDD());

            // Update gradient of total displacement
            gradD() = gradD().oldTime() + gradDD();

            // Calculate the stress using run-time selectable mechanical law
            const volScalarField DDEqnA("DEqnA", DDEqn.A());
            mechanical().correct(sigma());
        }
        while
        (
            !converged
            (
                iCorr,
#ifdef OPENFOAMESIORFOUNDATION
                mag(solverPerfDD.initialResidual()),
                cmptMax(solverPerfDD.nIterations()),
#else
                solverPerfDD.initialResidual(),
                solverPerfDD.nIterations(),
#endif
                DD()
            ) && ++iCorr < nCorr()
        );

        // Update point displacement increment
        mechanical().interpolate(DD(), pointDD());

        // Update point displacement
        pointD() = pointD().oldTime() + pointDD();

        // Update velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    // Store d2dt2
    rhoD2dt2D() = rho()*fvc::d2dt2(D());

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> linGeomSolid::tractionBoundarySnGrad
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
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

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
              - (n & (pSigma - pImpK*pGradDD))
            )*pRImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
