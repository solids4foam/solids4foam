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

#include "thermoLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicalModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, thermoLinearElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::thermoLinearElastic::thermoLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    linearElastic(name, mesh, dict),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    lambda_
    (
        planeStress()
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_(lambda_ + (2.0/3.0)*mu_),
    alpha_(dict.lookup("alpha"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoLinearElastic::~thermoLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::thermoLinearElastic::correct(volSymmTensorField& sigma)
{
    // Calculate linear elastic stress
    linearElastic::correct(sigma);

    // Lookup the temperature field from the solver
    const volScalarField T =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T", mesh());

    // Lookup the stress-free reference temperature field from the solver
    const volScalarField T0 =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T0", mesh());

    // Add thermal stress component
    sigma -= 3.0*K_*alpha_*(T - T0)*symmTensor(I);
}


void Foam::thermoLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate linear elastic stress
    linearElastic::correct(sigma);

    // Lookup the temperature field from the solver
    const volScalarField T =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T", mesh());

    const surfaceScalarField& Tf = fvc::interpolate(T);

    // Lookup the stress-free reference temperature field from the solver
    const volScalarField T0 =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("T0", mesh());

    const surfaceScalarField& T0f = fvc::interpolate(T0);

    // Add thermal stress component
    sigma -= 3.0*K_*alpha_*(Tf - T0f)*symmTensor(I);
}


// ************************************************************************* //
