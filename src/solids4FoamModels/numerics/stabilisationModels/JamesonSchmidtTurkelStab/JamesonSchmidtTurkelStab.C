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

#include "JamesonSchmidtTurkelStab.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(JamesonSchmidtTurkelStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, JamesonSchmidtTurkelStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::JamesonSchmidtTurkelStab::JamesonSchmidtTurkelStab
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

Foam::JamesonSchmidtTurkelStab::~JamesonSchmidtTurkelStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JamesonSchmidtTurkelStab::updateScalar
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
    computeDiffStencil(p, faceScalarPtr().ref(), scaleFactor_);
}


void Foam::JamesonSchmidtTurkelStab::updateVector
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
    computeDiffStencil(p, faceVectorPtr().ref(), scaleFactor_);
}

// ************************************************************************* //
