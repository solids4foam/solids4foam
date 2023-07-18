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

#include "thermalSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "fixedGradientFvPatchFields.H"
#include "thermalRobinFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalSolid, 0);
addToRunTimeSelectionTable(solidModel, thermalSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool thermalSolid::converged
(
    const int iCorr,
#ifdef OPENFOAMESIORFOUNDATION
    const SolverPerformance<scalar>& solverPerfT,
#else
    const lduSolverPerformance& solverPerfT,
#endif
    const volScalarField& T
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

     // Calculate relative residuals
    const scalar absResidualT =
        gMax
        (
#ifdef OPENFOAMESIORFOUNDATION
            DimensionedField<double, volMesh>
            (
                mag(T.internalField() - T.prevIter().internalField())
            )
#else
            mag(T.internalField() - T.prevIter().internalField())
#endif
        );
    const scalar residualT =
        absResidualT
       /max
        (
            gMax
            (
#ifdef OPENFOAMESIORFOUNDATION
                DimensionedField<double, volMesh>
                (
                    mag(T.internalField() - T.oldTime().internalField())
                )
#else
                    mag(T.internalField() - T.oldTime().internalField())
#endif
            ),
            SMALL
        );

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1)
    {
        bool convergedT = false;

        if
        (
            (
                solverPerfT.initialResidual() < solutionTol()
             && residualT < solutionTol()
            )
         || solverPerfT.initialResidual() < alternativeTol()
         || residualT < alternativeTol()
         || absResidualT < absTTol_
        )
        {
            convergedT = true;
        }


        if (convergedT)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (T), relRes (T), iters (T)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfT.initialResidual()
            << ", " << residualT
            << ", " << solverPerfT.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the enery loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalSolid::thermalSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    thermal_(mesh()),
    rhoC_
    (
        IOobject
        (
            "rhoC",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermal_.C()*mechanical().rho()
    ),
    k_(thermal_.k()),
    T_
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    gradT_
    (
        IOobject
        (
            "grad(T)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimTemperature/dimLength, vector::zero)
    ),
    absTTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "absoluteTemperatureTolerance",
            1e-06
        )
    )
{
    // Store T old time
    T_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool thermalSolid::evolve()
{
    Info<< "Evolving thermal solid solver" << endl;

    int iCorr = 0;
#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<scalar> solverPerfT;
    SolverPerformance<scalar>::debug = 0;
#else
    lduSolverPerformance solverPerfT;
    blockLduMatrix::debug = 0;
#endif

    Info<< "Solving energy equation for T"
        << endl;

    // Energy eq loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        T().storePrevIter();

        // Heat equation
        fvScalarMatrix TEqn
        (
            rhoC_*fvm::ddt(T_)
          - fvm::laplacian(k_, T_, "laplacian(k,T)")
        );

        // Under-relaxation the linear system
        TEqn.relax();

        // Solve the linear system
        solverPerfT = TEqn.solve();

        // Under-relax the field
        T_.relax();

        // Update gradient of temperature
        gradT_ = fvc::grad(T_);

        // Hack to avoid expensive copy of residuals
#ifdef OPENFOAMESI
        const_cast<dictionary&>(mesh().solverPerformanceDict()).clear();
#endif
    }
    while
    (
        !converged(iCorr, solverPerfT, T_)
     && ++iCorr < nCorr()
    );

#ifdef OPENFOAMESIORFOUNDATION
    SolverPerformance<scalar>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> thermalSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField(patch.size(), vector::zero)
    );
}


void thermalSolid::writeFields(const Time& runTime)
{
    Info<< "Max T = " << max(T_).value() << nl
        << "Min T = " << min(T_).value() << endl;

    // Heat flux
    volVectorField heatFlux
    (
        IOobject
        (
            "heatFlux",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       -k_*gradT_
    );

    Info<< "Max magnitude of heat flux = " << max(mag(heatFlux)).value()
        << endl;

    solidModel::writeFields(runTime);
}


void thermalSolid::setTemperatureAndHeatFlux
(
    fvPatchScalarField& temperaturePatch,
    const scalarField& temperature,
    const scalarField& heatFlux
)
{
    if (temperaturePatch.type() == thermalRobinFvPatchScalarField::typeName)
    {
        thermalRobinFvPatchScalarField& patchT =
            refCast<thermalRobinFvPatchScalarField>(temperaturePatch);

        patchT.temperature() = temperature;
        patchT.heatFlux() = heatFlux;
    }
    else
    {
        FatalErrorIn
        (
            "void thermalSolid::setTemperatureAndHeatFlux\n"
            "(\n"
            "    fvPatchScalarField& temperaturePatch,\n"
            "    const scalarField& temperature,\n"
            "    const scalarField& heatFlux\n"
            ")"
        )   << "Boundary condition "
            << temperaturePatch.type()
            << " for patch " << temperaturePatch.patch().name()
            << " should instead be type "
            << thermalRobinFvPatchScalarField::typeName
            << abort(FatalError);
    }
}


void thermalSolid::setTemperatureAndHeatFlux
(
    const label interfaceI,
    const label patchID,
    const scalarField& faceZoneTemperature,
    const scalarField& faceZoneHeatFlux
)
{
    const scalarField patchTemperature
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTemperature)
    );

    const scalarField patchHeatFlux
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneHeatFlux)
    );
    
#ifdef OPENFOAMESIORFOUNDATION
    setTemperatureAndHeatFlux
    (
        T_.boundaryFieldRef()[patchID],
        patchTemperature,
        patchHeatFlux
    );
#else
    setTemperatureAndHeatFlux
    (
        T_.boundaryField()[patchID],
        patchTemperature,
        patchHeatFlux
    );
#endif
}


void thermalSolid::setEqInterHeatTransferCoeff
(
    fvPatchScalarField& temperaturePatch,
    const scalarField& HTC
)
{
    if (temperaturePatch.type() == thermalRobinFvPatchScalarField::typeName)
    {
        thermalRobinFvPatchScalarField& patchT =
            refCast<thermalRobinFvPatchScalarField>(temperaturePatch);

        patchT.eqInterHeatTransferCoeff() = HTC;
    }
    else
    {
        FatalErrorIn
        (
            "void thermalSolid::setEqInterHeatTransferCoeff\n"
            "(\n"
            "    fvPatchScalarField& temperaturePatch,\n"
            "    const scalarField& HTc\n"
            ")"
        )   << "Boundary condition "
            << temperaturePatch.type()
            << " for patch " << temperaturePatch.patch().name()
            << " should instead be type "
            << thermalRobinFvPatchScalarField::typeName
            << abort(FatalError);
    }
}


void thermalSolid::setEqInterHeatTransferCoeff
(
    const label interfaceI,
    const label patchID,
    const scalarField& faceZoneHTC
)
{
    const scalarField patchHTC
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneHTC)
    );

#ifdef OPENFOAMESIORFOUNDATION
    setEqInterHeatTransferCoeff
    (
        T_.boundaryFieldRef()[patchID],
        patchHTC
    );
#else
    setEqInterHeatTransferCoeff
    (
        T_.boundaryField()[patchID],
        patchHTC
    );
#endif
}


tmp<scalarField> thermalSolid::faceZoneTemperature
(
    const label interfaceI
) const
{
    return globalPatches()[interfaceI].patchFaceToGlobal
    (
        T_.boundaryField()[globalPatches()[interfaceI].patch().index()]
    );
}

    
tmp<scalarField> thermalSolid::faceZoneHeatFlux
(
    const label interfaceI
) const
{
    scalarField patchHeatFlux =
        k_.boundaryField()[globalPatches()[interfaceI].patch().index()]*
        T_.boundaryField()[globalPatches()[interfaceI].patch().index()]
       .snGrad();
        
    return globalPatches()[interfaceI].patchFaceToGlobal(patchHeatFlux);
}


tmp<scalarField> thermalSolid::faceZoneHeatTransferCoeff
(
    const label interfaceI
) const
{
    const scalarField& patchDeltaCoeff =
        mesh().deltaCoeffs().boundaryField()
        [globalPatches()[interfaceI].patch().index()];
    
    scalarField patchHeatTransferCoeff =
        (1.0/patchDeltaCoeff)/
        k_.boundaryField()[globalPatches()[interfaceI].patch().index()];
        
    return globalPatches()[interfaceI].patchFaceToGlobal
    (
        patchHeatTransferCoeff
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
