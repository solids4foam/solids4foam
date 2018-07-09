/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "linGeomTotalDispSolid.H"
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

defineTypeNameAndDebug(linGeomTotalDispSolid, 0);
addToRunTimeSelectionTable(physicsModel, linGeomTotalDispSolid, solid);
addToRunTimeSelectionTable(solidModel, linGeomTotalDispSolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomTotalDispSolid::linGeomTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool linGeomTotalDispSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Mesh update loop
    do
    {
        int iCorr = 0;
        lduSolverPerformance solverPerfD;
        blockLduMatrix::debug = 0;

        Info<< "Solving the momentum equation for D" << endl;

        // Stablisation viscosity
        // const dimensionedScalar eta_ =
        //    solidModelDict().lookupOrDefault<dimensionedScalar>
        //    (
        //        "numericalViscosity",
        //        dimensionedScalar("eta", dimless/dimTime, 0.0)
        //    );

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DEqn
            (
                rho()*fvm::d2dt2(D())
                //+ eta_*rho_*fvm::ddt(D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + mechanical().RhieChowCorrection(D(), gradD())
            );

            // Under-relaxation the linear system
            DEqn.relax();

            // Solve the linear system
            solverPerfD = DEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

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
                solverPerfD.initialResidual(), solverPerfD.nIterations(), D()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    return true;
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
    const vectorField n = patch.nf();

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
