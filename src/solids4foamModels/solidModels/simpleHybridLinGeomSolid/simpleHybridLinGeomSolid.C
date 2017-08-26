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

#include "simpleHybridLinGeomSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleHybridLinGeomSolid, 0);
addToRunTimeSelectionTable(physicsModel, simpleHybridLinGeomSolid, solid);
addToRunTimeSelectionTable(solidModel, simpleHybridLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool simpleHybridLinGeomSolid::converged
(
    const int iCorr,
    const lduSolverPerformance& solverPerfD,
    const lduSolverPerformance& solverPerfp,
    const volVectorField& D,
    const volScalarField& p
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate displacement residual
    const scalar residualD =
        gMax
        (
            mag(D.internalField() - D.prevIter().internalField())
           /max
            (
                gMax(mag(D.internalField() - D.oldTime().internalField())),
                SMALL
            )
        );

    // Calculate pressure residual
    const scalar residualp =
        gMax
        (
            mag(p.internalField() - p.prevIter().internalField())
           /max
            (
                gMax(mag(p.internalField() - p.oldTime().internalField())),
                SMALL
            )
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        if
        (
            solverPerfD.initialResidual() < solutionTol()
         && residualD < solutionTol()
         && solverPerfp.initialResidual() < solutionTol()
         && residualp < solutionTol()
        )
        {
            Info<< "    All residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            residualD < alternativeTol()
         && residualp < alternativeTol()
        )
        {
            Info<< "    The relative residuals have converged" << endl;
            converged = true;
        }
        else if
        (
            solverPerfD.initialResidual() < alternativeTol()
         && solverPerfp.initialResidual() < alternativeTol()
        )
        {
            Info<< "    The solver residuals have converged" << endl;
            converged = true;
        }
        else
        {
            converged = false;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, resD, relResD, resP, relResP, matRes, "
            << "itersD, itersP" << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfD.initialResidual()
            << ", " << residualD
            << ", " << solverPerfp.initialResidual()
            << ", " << residualp
            << ", " << materialResidual
            << ", " << solverPerfD.nIterations()
            << ", " << solverPerfp.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within momentum loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

simpleHybridLinGeomSolid::simpleHybridLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    gradp_
    (
        IOobject
        (
            "grad(" + p_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimPressure/dimLength, vector::zero)
    ),
    rK_(1.0/mechanical().K()),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool simpleHybridLinGeomSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfD;
    lduSolverPerformance solverPerfp;
    blockLduMatrix::debug = 0;

    Info<< "Solving the momentum equation for D and p using a transient "
        << "SIMPLE approach" << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();
        p_.storePrevIter();

        // Transient SIMPLE segregated solution methodlogy

        // Linear momentum equation total displacement form with pressure
        // term removed i.e. we add it back on
        tmp<fvVectorMatrix> HDEqn
        (
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(sigma(), "div(sigma)")
          + gradp_
          + rho()*g()
          + mechanical().RhieChowCorrection(D(), gradD())
        );

        // Under-relaxation the linear system
        HDEqn().relax();

        // Solve the linear system
        solverPerfD =
            solve
            (
                HDEqn()
                ==
                -gradp_
            );

        // Under-relax the field
        relaxField(D(), iCorr);

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update p boundaries in case they depend on D
        p_.boundaryField().updateCoeffs();

        // Prepare clean 1/Ap without contribution from under-relaxation
        const volScalarField rDA("(1|A(D))", 1/HDEqn().A());

        // Clear equation to reduce memory overhead
        HDEqn.clear();

        // Pressure equation with Rhie-Chow corrections to avoid
        // oscillations
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rDA, p_, "laplacian(Dp,p)")
          - fvc::div(rDA*gradp_, "skewCorrectedDiv(D)")
         == fvc::div(D(), "skewCorrectedDiv(D)")
          + fvm::Sp(rK_, p_)
        );

        // Under-relax the equation
        pEqn.relax();

        // Solve the linear system
        solverPerfp = pEqn.solve();

        // Explicitly relax pressure for momentum corrector
        p_.relax();

        // Update gradient of pressure
        gradp_ = fvc::grad(p_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    while
    (
        !converged(iCorr, solverPerfD, solverPerfp, D(), p_)
     && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Velocity
    U() = fvc::ddt(D());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    return true;
}


tmp<vectorField> simpleHybridLinGeomSolid::tractionBoundarySnGrad
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
