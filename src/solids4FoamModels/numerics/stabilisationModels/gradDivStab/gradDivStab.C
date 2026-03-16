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

#include "gradDivStab.H"
#include "addToRunTimeSelectionTable.H"
#include "orthogonalSnGrad.H"
#include "fvc.H"
#include "compatibilityFunctions.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gradDivStab, 0);
    addToRunTimeSelectionTable
    (
        stabilisationModel, gradDivStab, stabModel
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::gradDivStab::gradDivStab
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

Foam::gradDivStab::~gradDivStab()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradDivStab::updateVector
(
    const volVectorField& p,
    const volTensorField* gradPtr
) const
{
    clearCellVectorCache();

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

    // Create the orthogonal normal gradient operator
    fv::orthogonalSnGrad<scalar> op(mesh());

    // Inner pass: divergence of field
    const volScalarField divP(fvc::div(p));

    // Face-based h^2 scale via deltaCoeffs
    const surfaceScalarField h2f(1.0/sqr(mesh().deltaCoeffs()));

    // Face unit normals
    surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Outer pass: orthogonal diffusion of L1, scaled by h^2
    // Note that we ignore the tangential component as we only care about
    // stabilisation (not the physical operator) and the tangential components
    // could inject noise
    tmp<surfaceVectorField> trhs
    (
        -scaleFactor_*Foam::pow(h2f, 0.5)*op.snGrad(divP)*n
    );
#ifdef OPENFOAM_COM
    tmpRef(trhs).setOriented(true);
#endif

    autoPtrRef(faceVectorPtr()) = tmpRef(trhs);
}

// ************************************************************************* //
