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

#include "alphaStab.H"
#include "addToRunTimeSelectionTable.H"
#include "compatibilityFunctions.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(alphaStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, alphaStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::alphaStab::alphaStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::alphaStab::~alphaStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaStab::updateScalar
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

    // If required, initialise the face stabilisation field
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

        // This is an oriented field
#ifdef OPENFOAM_COM
        autoPtrRef(faceScalarPtr()).setOriented(true);
#endif
    }

    const solidModel& solMod = lookupSolidModel(p.mesh());

    if (solMod.highOrderResidual())
    {
#ifndef FOAMEXTEND
        computeDiffStencilHighOrderScalar
        (
            p,
            gradP,
            autoPtrRef(faceScalarPtr()),
            scaleFactor_
        );
#else
        computeDiffStencil(p, gradP, autoPtrRef(faceScalarPtr()), scaleFactor_);
#endif
    }
    else
    {
        // Update the stabilisation
        computeDiffStencil(p, gradP, autoPtrRef(faceScalarPtr()), scaleFactor_);
    }
}


void Foam::alphaStab::updateVector
(
    const volVectorField& p,
    const volTensorField* gradPtr
) const
{
    clearCellVectorCache();

    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << p.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    const volTensorField& gradP = *gradPtr;

    // If required, initialise the face stabilisation field
    if (faceVectorPtr().empty())
    {
        faceVectorPtr().set
        (
            new surfaceVectorField
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
                dimensionedVector("0", dims(), vector::zero)
            )
        );

        // This is an oriented field
#ifdef OPENFOAM_COM
        autoPtrRef(faceVectorPtr()).setOriented(true);
#endif
    }

    const solidModel& solMod = lookupSolidModel(p.mesh());

    if (solMod.highOrderResidual())
    {
#ifndef FOAMEXTEND
        // Update the high-order stabilisation
        computeDiffStencilHighOrderVector
        (
            p,
            gradP,
            autoPtrRef(faceVectorPtr()),
            scaleFactor_
        );
#else
        computeDiffStencil
        (
            p,
            gradP,
            autoPtrRef(faceVectorPtr()),
            scaleFactor_
        );
#endif
    }
    else
    {
        // Update the second-order stabilisation
        computeDiffStencil
        (
            p,
            gradP,
            autoPtrRef(faceVectorPtr()),
            scaleFactor_
        );
    }
}


const Foam::fvVectorMatrix& Foam::alphaStab::vectorJacobian
(
    const volVectorField& field,
    const surfaceScalarField* gammaPtr,
    const bool rebuild
) const
{
    checkGamma(gammaPtr);

    // Vector Jacobian taken from diffStencilLaplacianStab
    // Residual can use alpha stabilisation but Jacobian is build using
    // implicit Laplacian. Alpha implicit Laplacian will be similar so this
    // is ok workaround
    if (vectorJacobianPtr().empty() || rebuild)
    {
        if (gammaPtr == nullptr)
        {
            vectorJacobianPtr().reset
            (
                new fvVectorMatrix
                (
                    scaleFactorJacobian()
                   *fvm::laplacian(field, "laplacian(" + field.name() + ")")
                )
            );
        }
        else
        {
            const word schemeName
            (
                "laplacian(" + gammaPtr->name() + "," + field.name() + ")"
            );

            vectorJacobianPtr().reset
            (
                new fvVectorMatrix
                (
                    scaleFactorJacobian()
                   *fvm::laplacian(*gammaPtr, field, schemeName)
                )
            );
        }
    }

    return vectorJacobianPtr()();
}

// ************************************************************************* //
