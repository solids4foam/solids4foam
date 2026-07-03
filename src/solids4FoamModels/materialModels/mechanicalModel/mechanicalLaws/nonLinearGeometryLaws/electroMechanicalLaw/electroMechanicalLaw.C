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
    f0_
    (
        IOobject
        (
            "f0",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    f0f0_("f0f0", sqr(f0_)),
    f0f_
    (
        IOobject
        (
            "f0f",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    f0f0f_("f0f0f", sqr(f0f_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroMechanicalLaw::~electroMechanicalLaw()
{}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::impK() const
{
    return passiveMechLawPtr_->impK();
}


void Foam::electroMechanicalLaw::materialTangentField(List<mat66>& matTan) const
{
    passiveMechLawPtr_->materialTangentField(matTan);
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::bulkModulus() const
{
    return passiveMechLawPtr_->bulkModulus();
}


Foam::tmp<Foam::volScalarField> Foam::electroMechanicalLaw::shearModulus() const
{
    return passiveMechLawPtr_->shearModulus();
}


void Foam::electroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    if (!mesh().foundObject<volScalarField>("Ta"))
    {
        FatalErrorInFunction
            << "electroMechanicalLaw requires a volScalarField named Ta in "
            << "the solid mesh objectRegistry. Ta must be supplied by the "
            << "electromechanical coupling, e.g. LandNiederer active tension."
            << abort(FatalError);
    }

    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Take a reference to the deformation gradient
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    const volScalarField& Ta =
        mesh().lookupObject<volScalarField>("Ta");

    // Add active stress: convert 2nd Piola-Kirchhoff to Cauchy
    sigma += symm(F & (Ta*f0f0_) & F.T())/J;
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    if (!mesh().foundObject<volScalarField>("Ta"))
    {
        FatalErrorInFunction
            << "electroMechanicalLaw requires a volScalarField named Ta in "
            << "the solid mesh objectRegistry. Ta must be supplied by the "
            << "electromechanical coupling, e.g. LandNiederer active tension."
            << abort(FatalError);
    }

    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Take a reference to the deformation gradient
    const surfaceTensorField& F = this->Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    const volScalarField& Ta =
        mesh().lookupObject<volScalarField>("Ta");

    const surfaceScalarField Taf(fvc::interpolate(Ta));

    // Add active stress: convert 2nd Piola-Kirchhoff to Cauchy
    sigma += symm(F & (Taf*f0f0f_) & F.T())/J;
}


// ************************************************************************* //
