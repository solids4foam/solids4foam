/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "linearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IOdictionary.H"
#include "mechanicalModel.H"

// WIP
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElastic::linearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    lambda_
    (
        planeStress()
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    mu_(E_/(2.0*(1.0 + nu_)))
{
    // PC: some mechanical laws are only appropriate for some solidModels
    // what is a nice way to check this?
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElastic::~linearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElastic::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::linearElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            2.0*mu_ + lambda_
        )
    );
}


const Foam::dimensionedScalar& Foam::linearElastic::mu() const
{
    return mu_;
}


const Foam::dimensionedScalar& Foam::linearElastic::lambda() const
{
    return lambda_;
}


void Foam::linearElastic::correct(volSymmTensorField& sigma)
{
    // Lookup the strain tensor from the solver
    // const volSymmTensorField epsilon =
    //     mesh().db().lookupObject<fvMesh>
    //     (
    //         baseMeshRegionName()
    //     ).lookupObject<mechanicalModel>
    //     (
    //         "mechanicalProperties"
    //     ).lookupBaseMeshVolField<symmTensor>("epsilon", mesh());

    // WIP
    // What is a nice way to do this?
    // We cannot use lookupBaseMeshField because it linearly interpolates fields
    // to bi-material interfaces, which is bad.
    // This method is OK but we need to re-evaluate gradD; in the case of finite
    // strain approaches, we would have to re-evaluate F, relF, etc.
    // Hmnn...
    const volVectorField& D = mesh().lookupObject<volVectorField>("D");
    const volSymmTensorField epsilon = symm(fvc::grad(D));

    // Calculate stress based on Hooke's law
    sigma = 2.0*mu_*epsilon + lambda_*tr(epsilon)*I;
}


void Foam::linearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Lookup the strain tensor from the solver
    const surfaceSymmTensorField epsilon =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshSurfaceField<symmTensor>("epsilonf", mesh());

    // Calculate stress based on Hooke's law
    sigma = 2.0*mu_*epsilon + lambda_*tr(epsilon)*I;
}


// ************************************************************************* //
