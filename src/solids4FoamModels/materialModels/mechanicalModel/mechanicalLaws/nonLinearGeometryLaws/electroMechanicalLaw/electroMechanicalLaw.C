/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "electroMechanicalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroMechanicalLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, electroMechanicalLaw, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

// void Foam::electroMechanicalLaw::makeP0f() const
// {
//     if (p0fPtr_)
//     {
//         FatalErrorIn("void Foam::electroMechanicalLaw::makeP0f() const")
//             << "pointer already set" << abort(FatalError);
//     }

//     p0fPtr_ =
//         new surfaceScalarField
//         (
//             "p0f",
//             fvc::interpolate(p0_)
//         );
// }


// const Foam::surfaceScalarField& Foam::electroMechanicalLaw::p0f() const
// {
//     if (!p0fPtr_)
//     {
//         makeP0f();
//     }

//     return *p0fPtr_;
// }


// const Foam::volScalarField& Foam::electroMechanicalLaw::lookupPressureField() const
// {
//     if (mesh().thisDb().parent().foundObject<objectRegistry>(pRegion_))
//     {
//         return mesh().thisDb().parent().subRegistry
//         (
//             pRegion_
//         ).lookupObject<volScalarField>(pName_);
//     }
//     else if
//     (
//         mesh().thisDb().parent().foundObject<objectRegistry>("solid")
//     )
//     {
//         return mesh().thisDb().parent().subRegistry
//         (
//             "solid"
//         ).lookupObject<volScalarField>(pName_);
//     }
//     else
//     {
//         FatalErrorIn("Foam::electroMechanicalLaw::lookupPressureField()")
//             << "Cannot find " << pName_ << " field in " << pRegion_
//             << " or in 'solid'" << abort(FatalError);
//     }

//     // Keep compiler happy
//     return mesh().lookupObject<volScalarField>("null");
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::electroMechanicalLaw::electroMechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    passiveMechLawPtr_
    (
        mechanicalLaw::NewNonLinGeomMechLaw
        (
            word(dict.subDict("passiveMechanicalLaw").lookup("type")),
            mesh,
            dict.subDict("passiveMechanicalLaw"),
            nonLinGeom
        )
    ),
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime")))
{
    if (rampTime_ < 0.0)
    {
        FatalErrorIn("electroMechanicalLaw::electroMechanicalLaw(...)")
            << "rampTime should be greater than or equal to zero"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroMechanicalLaw::~electroMechanicalLaw()
{}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::impK() const
{
    return passiveMechLawPtr_->impK();
}


void Foam::electroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Lookup the fibre directions
    // How best should we do this to avoid duplicating the fibre field?
    // For now, let's hard-code in the field name
    const volSymmTensorField f0f0 =
        mesh().lookupObject<volSymmTensorField>("f0f0");

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // For now, we will assume a constant active stress
    // The next step will be to include an active-stress model to convert
    // muscle activation to fibre tension

    // Calculate current value of Ta
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    // Add active stress to the passive stress
    // Note that the active stress is converted from a 2nd Piola-Kirchhoff
    // stress to a Cauchy stress
    sigma += J*symm(F & (currentTa*f0f0) & F.T());
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    notImplemented
    (
        "void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)"
    );
}


// ************************************************************************* //
