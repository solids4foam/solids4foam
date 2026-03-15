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

#include "generalisedEvenOrderLaplacianStab.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(generalisedEvenOrderLaplacianStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, generalisedEvenOrderLaplacianStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::generalisedEvenOrderLaplacianStab::generalisedEvenOrderLaplacianStab
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dimensionSet& dims
)
:
    stabilisationModel(mesh, dict, dims),
    scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
    laplacianPower_(readScalar(dict.lookup("laplacianPower")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::generalisedEvenOrderLaplacianStab::~generalisedEvenOrderLaplacianStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::generalisedEvenOrderLaplacianStab::updateScalar
(
    const volScalarField& p,
    const volVectorField* gradPtr
) const
{
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
    }

    // Update the stabilisation
    computeDiffStencil(p, faceScalarPtr().ref(), scaleFactor_, laplacianPower_);
}


void Foam::generalisedEvenOrderLaplacianStab::updateVector
(
    const volVectorField& p,
    const volTensorField* gradPtr
) const
{
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
    }

    // Update the stabilisation
    computeDiffStencil(p, faceVectorPtr().ref(), scaleFactor_, laplacianPower_);
}

// ************************************************************************* //
