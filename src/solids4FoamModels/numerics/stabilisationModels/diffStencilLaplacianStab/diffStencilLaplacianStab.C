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

#include "diffStencilLaplacianStab.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffStencilLaplacianStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, diffStencilLaplacianStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::diffStencilLaplacianStab::diffStencilLaplacianStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
    scaleFactorJacobian_
    (
        dict.lookupOrDefault<scalar>("scaleFactorJacobian", scaleFactor_)
    ),
    faceScalarPtr_(),
    faceVectorPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffStencilLaplacianStab::~diffStencilLaplacianStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diffStencilLaplacianStab::updateScalar
(
    const volScalarField& p,
    const volVectorField* gradPtr
) const
{
    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << p.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    const volVectorField& gradP = *gradPtr;

    // If required, initialise the face stabilisation field
    if (faceScalarPtr_.empty())
    {
        faceScalarPtr_.set
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

    // Update the stabilisation
    computeDiffStencil(p, gradP, faceScalarPtr_.ref(), scaleFactor_);
}


void Foam::diffStencilLaplacianStab::updateVector
(
    const volVectorField& p,
    const volTensorField* gradPtr
) const
{
    if (gradPtr == nullptr)
    {
        FatalErrorInFunction
            << "grad(" << p.name() << ") must be provided with this "
            << "stabilisation method" << exit(FatalError);
    }

    const volTensorField& gradP = *gradPtr;

    // If required, initialise the face stabilisation field
    if (faceVectorPtr_.empty())
    {
        faceVectorPtr_.set
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
    }

    // Update the stabilisation
    computeDiffStencil(p, gradP, faceVectorPtr_.ref(), scaleFactor_);
}

const Foam::fvScalarMatrix& Foam::diffStencilLaplacianStab::scalarJacobian
(
    const volScalarField& field,
    const surfaceScalarField* gammaPtr
) const
{
    // If required, initialise the face stabilisation field
    if (scalarJacobianPtr_.empty())
    {
        if (gammaPtr == nullptr)
        {
            scalarJacobianPtr_.set
            (
                new fvScalarMatrix
                (
                    fvm::laplacian(field, "laplacian(" + field.name() + ")")
                )
            );
        }
        else
        {
            const word schemeName
            (
                "laplacian(" + gammaPtr->name() + "," + field.name() + ")"
            );

            scalarJacobianPtr_.set
            (
                new fvScalarMatrix
                (
                    fvm::laplacian(*gammaPtr, field, schemeName)
                )
            );
        }
    }

    return scalarJacobianPtr_();
}


const Foam::fvVectorMatrix& Foam::diffStencilLaplacianStab::vectorJacobian
(
    const volVectorField& field,
    const surfaceScalarField* gammaPtr
) const
{
    // If required, initialise the face stabilisation field
    if (vectorJacobianPtr_.empty())
    {
        if (gammaPtr == nullptr)
        {
            vectorJacobianPtr_.set
            (
                new fvVectorMatrix
                (
                    scaleFactorJacobian_
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

            vectorJacobianPtr_.set
            (
                new fvVectorMatrix
                (
                    scaleFactorJacobian_
                   *fvm::laplacian(*gammaPtr, field, schemeName)
                )
            );
        }
    }

    return vectorJacobianPtr_();
}


// ************************************************************************* //
