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
    Ta_(dict.lookup("activeTension"))
    // b_
    // (
    //     mechanicalLaw::dict().lookupOrAddDefault<dimensionedScalar>
    //     (
    //         "biotCoeff", dimensionedScalar("0", dimless, 1.0)
    //     )
    // ),
    // pName_(mechanicalLaw::dict().lookupOrAddDefault<word>("pressureFieldName", "p")),
    // pRegion_(mechanicalLaw::dict().lookupOrAddDefault<word>("pressureFieldRegion", "region0")),
    // p0_
    // (
    //     IOobject
    //     (
    //         "p0",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     mechanicalLaw::dict().lookupOrAddDefault<dimensionedScalar>
    //     (
    //         "p0",
    //         dimensionedScalar("zero", dimPressure, 0.0)
    //     )
    // ),
    // p0fPtr_(NULL)
{
    // if (gMax(mag(p0_)()) > SMALL)
    // {
    //     Info<< "Reading p0 initial/residual pore-pressure field" << endl;
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroMechanicalLaw::~electroMechanicalLaw()
{
    // deleteDemandDrivenData(p0fPtr_);
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::impK() const
{
    return passiveMechLawPtr_->impK();
}


Foam::tmp<Foam::Field<Foam::RectangularMatrix<Foam::scalar>>>
Foam::electroMechanicalLaw::materialTangentField() const
{
    return passiveMechLawPtr_->materialTangentField();
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::bulkModulus() const
{
    return passiveMechLawPtr_->bulkModulus();
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

    // Add active stress to the passive stress
    // Note that the active stress is converted from a 2nd Piola-Kirchhoff
    // stress to a Cauchy stress
    sigma += J*symm(F & (Ta_*f0f0) & F.T());
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Lookup the fibre directions
    // How best should we do this to avoid duplicating the fibre field?
    // For now, let's hard-code in the field name
    const surfaceSymmTensorField f0f0 =
        mesh().lookupObject<surfaceSymmTensorField>("f0f0f");

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const surfaceTensorField& F = this->Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    // For now, we will assume a constant active stress
    // The next step will be to include an active-stress model to convert
    // muscle activation to fibre tension

    // Add active stress to the passive stress
    // Note that the active stress is converted from a 2nd Piola-Kirchhoff
    // stress to a Cauchy stress
    sigma += J*symm(F & (Ta_*f0f0) & F.T());
}


// ************************************************************************* //
