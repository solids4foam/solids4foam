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

#include "poroLinGeomSolid.H"
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

defineTypeNameAndDebug(poroLinGeomSolid, 0);
addToRunTimeSelectionTable(solidModel, poroLinGeomSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool poroLinGeomSolid::converged
(
    const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
    const SolverPerformance<vector>& solverPerfD,
    const SolverPerformance<scalar>& solverPerfp,
#else
    const lduSolverPerformance& solverPerfD,
    const lduSolverPerformance& solverPerfp,
#endif
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
#ifdef OPENFOAM_NOT_EXTEND
            DimensionedField<double, volMesh>
#endif
            (
                mag(D.internalField() - D.prevIter().internalField())
               /max
                (
                    gMax
                    (
#ifdef OPENFOAM_NOT_EXTEND
                        DimensionedField<double, volMesh>
#endif
                        (
                            mag(D.internalField() - D.oldTime().internalField())
                        )
                    ),
                    SMALL
                )
            )
        );

    // Calculate pressure residual
    const scalar residualp =
        gMax
        (
#ifdef OPENFOAM_NOT_EXTEND
            DimensionedField<double, volMesh>
#endif
            (
                mag(p.internalField() - p.prevIter().internalField())
               /max
                (
                    gMax
                    (
#ifdef OPENFOAM_NOT_EXTEND
                        DimensionedField<double, volMesh>
#endif
                        (
                            mag(p.internalField() - p.oldTime().internalField())
                        )
                    ),
                    SMALL
                )
            )
        );

    // Calculate material residual
    const scalar materialResidual = this->materialResidual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        if
        (
            mag(solverPerfD.initialResidual()) < solutionTol()
         && solverPerfp.initialResidual() < solutionTol()
         && residualD < solutionTol()
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
            mag(solverPerfD.initialResidual()) < alternativeTol()
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
        Info<< "    Corr, resD, resP, relResD, relResP, matRes, iters" << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << mag(solverPerfD.initialResidual())
            << ", " << solverPerfp.initialResidual()
            << ", " << residualD
            << ", " << residualp
            << ", " << materialResidual
            << ", " << solverPerfD.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within momentum-pressure loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

poroLinGeomSolid::poroLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    useMechanicalConstitutiveLawManager_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "useMechanicalConstitutiveLawManager", false
        )
    ),
    mechanicalManagerPtr_(),
    impK_(makeImpK()),
    impKf_(makeImpKf()),
    rImpK_(1.0/impK_),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    gradp_
    (
        IOobject
        (
            "grad(p)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimPressure/dimLength, vector::zero)
    ),
    hydraulicConductivity_(solidModelDict().lookup("hydraulicConductivity")),
    gammaWater_(solidModelDict().lookup("waterSpecificWeight")),
    porosity_(solidModelDict().lookup("porosity")),
    saturation_(solidModelDict().lookup("degreeOfSaturation")),
    KWater_(solidModelDict().lookup("waterBulkModulus")),
    rKprime_
    (
        (saturation_/KWater_)
      + (1.0 - saturation_)
        /dimensionedScalar("atmosphericPressure", dimPressure, 1e+05)
    )
{
    // Store old time of p
    p_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::mechanicalConstitutiveLawManager&
poroLinGeomSolid::mechanicalManager() const
{
    if (mechanicalManagerPtr_.empty())
    {
        // mechanicalModel is itself the mechanicalProperties IOdictionary, so
        // both frameworks are built from exactly the same entries
        mechanicalManagerPtr_.set
        (
            new mechanicalConstitutiveLawManager(mesh(), mechanical())
        );
    }

    return mechanicalManagerPtr_();
}


void poroLinGeomSolid::correctStress()
{
    if (!useMechanicalConstitutiveLawManager_)
    {
        mechanical().correct(sigma());
        return;
    }

    // The framework is a pure function of the displacement gradient and the
    // old-time state, so the gradient is passed explicitly rather than looked
    // up from the registry, and the old-time state is rolled over by the
    // manager rather than by a separate call
    mechanicalManager().updateStressSmallStrain
    (
        gradD(),
        gradD().oldTime(),
        mesh().time().deltaTValue(),
        sigma()
    );
}


Foam::tmp<Foam::volScalarField> poroLinGeomSolid::makeImpK() const
{
    if (!useMechanicalConstitutiveLawManager_)
    {
        return mechanical().impK();
    }

    return frameworkImpK(mechanicalManager(), tangentRequest::scalar);
}


Foam::tmp<Foam::surfaceScalarField> poroLinGeomSolid::makeImpKf() const
{
    if (!useMechanicalConstitutiveLawManager_)
    {
        return mechanical().impKf();
    }

    // The framework has no separate face tangent: the face value is the
    // interpolate of the cell one
    return fvc::interpolate(makeImpK()());
}


Foam::scalar poroLinGeomSolid::materialResidual()
{
    if (!useMechanicalConstitutiveLawManager_)
    {
        return mechanical().residual();
    }

    // The framework keeps its own state and rolls it over itself, so it has no
    // residual of its own to report and contributes nothing to convergence
    return 0.0;
}


void poroLinGeomSolid::updateTotalFields()
{
    solidModel::updateTotalFields();

    // The framework keeps its own state, and its laws may have end-of-step
    // work or diagnostics
    if (useMechanicalConstitutiveLawManager_)
    {
        mechanicalManager().endTimeStep();
    }
}


bool poroLinGeomSolid::evolve()
{
    Info << "Evolving poro solid solver" << endl;

    int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<scalar> solverPerfp;
    SolverPerformance<scalar>::debug = 0;
    SolverPerformance<vector>::debug = 0;
#else
    lduSolverPerformance solverPerfD;
    lduSolverPerformance solverPerfp;
    blockLduMatrix::debug = 0;
#endif

    Info<< "Solving the pressure equation for p and momentum equation for D"
        << endl;

    // Pressure-displacement coupling outer loop
    do
    {
        // Pressure equation

        // Store fields for under-relaxation and residual calculation
        p_.storePrevIter();

        // Pressure equation
        fvScalarMatrix pEqn
        (
            (porosity_*rKprime_)*fvm::ddt(p_)
          + fvc::div(U())
         == (hydraulicConductivity_/gammaWater_)*fvm::laplacian(p_)
        );

        // Under-relaxation the linear system
        pEqn.relax();

        // Solve the linear system
        solverPerfp = pEqn.solve();

        // Under-relax the field
        p_.relax();

        // Update gradient of pressure
        gradp_ = fvc::grad(p_);


        // Momentum equation

        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        momentumStabilisation().updateVector(D(), &gradD());

        // Linear momentum equation total displacement form
        // Assemble the RHS in stages.
        tmp<fvVectorMatrix> tRhsEqn
        (
            fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
        );
        tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, D(), "laplacian(DD,D)");
        tmpRef(tRhsEqn) += fvc::div(sigma(), "div(sigma)");
        tmpRef(tRhsEqn) += rho()*g();
        tmpRef(tRhsEqn) += impK_*momentumStabilisation().cellVector(nullptr, true);

        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D()) == tRhsEqn
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Under-relax the field
        relaxField(D(), iCorr);

        // Update increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Update velocity as it is used in the pEqn
        U() = fvc::ddt(D());

        // Calculate the stress using run-time selectable mechanical law
        correctStress();

        // Update impKf to improve convergence
        // Note: impK and rImpK are not updated as they are used for traction
        // boundaries
        if (iCorr % 10 == 0)
        {
            impKf_ = makeImpKf();
        }
    }
    while
    (
        !converged(iCorr, solverPerfD, solverPerfp, D(), p_)
     && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<scalar>::debug = 1;
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> poroLinGeomSolid::tractionBoundarySnGrad
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
