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

#include "adaptiveRhieChowStab.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvm.H"
#include "compatibilityFunctions.H"

namespace Foam
{
    defineTypeNameAndDebug(adaptiveRhieChowStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, adaptiveRhieChowStab, stabModel
    );
}


Foam::adaptiveRhieChowStab::adaptiveRhieChowStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    pressureSmoothingCoeff_
    (
        dict.lookupOrDefault<scalar>
        (
            "pressureSmoothingCoeff",
            dict.lookupOrDefault<scalar>("scaleFactor", 1.0)
        )
    ),
    adaptiveGammaPtr_(),
    adaptiveGammaTimeIndex_(-1)
{}


Foam::adaptiveRhieChowStab::~adaptiveRhieChowStab()
{}


void Foam::adaptiveRhieChowStab::updateScalar
(
    const volScalarField& p,
    const volVectorField* gradPtr
) const
{
    clearCellScalarCache();

    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << p.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    const volVectorField& gradP = *gradPtr;

    if (faceScalarPtr().empty())
    {
        faceScalarPtr().set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "faceStabilisation(" + p.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("0", dims(), 0.0)
            )
        );
    }

    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    autoPtrRef(faceScalarPtr()) =
        fvc::snGrad(p) - (n & fvc::interpolate(gradP));
}


const Foam::fvScalarMatrix& Foam::adaptiveRhieChowStab::scalarJacobian
(
    const volScalarField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    if (scalarJacobianPtr().empty() || rebuild)
    {
        if (gammaPtr == nullptr)
        {
            FatalErrorInFunction
                << "adaptiveRhieChow requires an adaptive gamma field for the "
                << "scalar Jacobian" << exit(FatalError);
        }

        const word schemeName
        (
            "laplacian(" + gammaPtr->name() + "," + field.name() + ")"
        );

        scalarJacobianPtr().reset
        (
            new fvScalarMatrix
            (
                scaleFactorJacobian()
               *fvm::laplacian(*gammaPtr, field, schemeName)
            )
        );
    }

    return scalarJacobianPtr()();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::adaptiveRhieChowStab::adaptiveScalarGamma
(
    const volVectorField& D,
    const surfaceScalarField& impKf,
    const volScalarField& rho,
    const dimensionedScalar& dampingCoeff
) const
{
    if
    (
        adaptiveGammaPtr_.empty()
     || adaptiveGammaTimeIndex_ != mesh().time().timeIndex()
    )
    {
        fvVectorMatrix approxJ
        (
            fvm::laplacian(impKf, D, "laplacian(DD,D)")
          - rho*fvm::d2dt2(D)
        );

        if (dampingCoeff.value() > SMALL)
        {
            approxJ -= dampingCoeff*rho*fvm::ddt(D);
        }

        // Match the legacy implementation by caching the adaptive
        // coefficient instead of rebuilding it every nonlinear evaluation.
        approxJ.relax();

        adaptiveGammaPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "adaptiveScalarGamma(" + D.name() + ")",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                -pressureSmoothingCoeff_/fvc::interpolate(approxJ.A())
            )
        );

        adaptiveGammaTimeIndex_ = mesh().time().timeIndex();
    }

    return tmp<surfaceScalarField>(new surfaceScalarField(adaptiveGammaPtr_()));
}

// ************************************************************************* //
