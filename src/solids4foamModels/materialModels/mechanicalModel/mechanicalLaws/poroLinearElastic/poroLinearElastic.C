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

#include "poroLinearElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroLinearElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, poroLinearElastic, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::poroLinearElastic::poroLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    linearElastic(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::poroLinearElastic::~poroLinearElastic()
{}


void Foam::poroLinearElastic::correct(volSymmTensorField& sigma)
{
    // Calculate effective stress
    linearElastic::correct(sigma);

    // Lookup the pressure field from the solver
    const volScalarField p =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("p", mesh());

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma -= p*symmTensor(I);
}


void Foam::poroLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate effective stress
    linearElastic::correct(sigma);

    // Lookup the pressure field from the solver
    const volScalarField p =
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        ).lookupBaseMeshVolField<scalar>("p", mesh());

    const surfaceScalarField pf = fvc::interpolate(p);

    // Calculate the total stress as the sum of the effective stress and the
    // pore-pressure
    sigma -= pf*symmTensor(I);
}


// ************************************************************************* //
