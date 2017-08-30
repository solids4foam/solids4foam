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
        mechanicalLaw, thermoLinearElastic, linGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::thermoLinearElastic::thermoLinearElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    linearElastic(name, mesh, dict, nonLinGeom),
    alpha_(dict.lookup("alpha")),
    T0_(dict.lookup("T0")),
    TPtr_(NULL)
{
    // Check if the T field exists in memory
    // THIS IS NOT CORRECT: of the mechanical model is created first then T will
    // not exist yet
    // One option would be to check again inside the correct function
    // Anther problem is the T thermalConvection gives out looking for 'k' if we
    // read T here
    // if (!mesh.foundObject<volScalarField>("T"))
    // {
    //     Info<< type() << ": reading T field from " << mesh.time().timeName()
    //         << endl;

    //     // Read the temperature field from disk
    //     TPtr_.set
    //     (
    //         new volScalarField
    //         (
    //             IOobject
    //             (
    //                 "T",
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::MUST_READ,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh
    //         )
    //     );
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoLinearElastic::~thermoLinearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::thermoLinearElastic::correct(volSymmTensorField& sigma)
{
    // Calculate linear elastic stress
    linearElastic::correct(sigma);

    if (TPtr_.valid())
    {
        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(TPtr_() - T0_)*symmTensor(I);
    }
    else
    {
        // Lookup the temperature field from the solver
        const volScalarField& T = mesh().lookupObject<volScalarField>("T");

        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(T - T0_)*symmTensor(I);
    }
}


void Foam::thermoLinearElastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate linear elastic stress
    linearElastic::correct(sigma);

    if (TPtr_.valid())
    {
        // Interpolate the volField temperature to the faces
        const surfaceScalarField Tf = fvc::interpolate(TPtr_());

        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(Tf - T0_)*symmTensor(I);
    }
    else
    {
        // Lookup the temperature field from the solver
        const volScalarField& T = mesh().lookupObject<volScalarField>("T");

        // Interpolate the volField temperature to the faces
        const surfaceScalarField Tf = fvc::interpolate(T);

        // Add thermal stress component
        sigma -= 3.0*K()*alpha_*(Tf - T0_)*symmTensor(I);
    }
}


// ************************************************************************* //
