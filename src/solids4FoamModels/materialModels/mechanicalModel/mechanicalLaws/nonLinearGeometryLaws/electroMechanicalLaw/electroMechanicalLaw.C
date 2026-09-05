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
    f0f0f_("f0f0f", sqr(f0f_)),
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime"))),
    useFieldTa_(false),
    fieldTaChecked_(false)
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
    // Lazy check for field-based active tension
    if (!fieldTaChecked_)
    {
        fieldTaChecked_ = true;
        useFieldTa_ = mesh().foundObject<volScalarField>("Ta");

        if (useFieldTa_)
        {
            Info<< "    electroMechanicalLaw: using field-based active tension"
                << " (Ta volScalarField from objectRegistry)" << endl;
        }
        else
        {
            Info<< "    electroMechanicalLaw: using constant active tension"
                << " Ta = " << Ta_.value()
                << " with rampTime = " << rampTime_ << endl;
        }
    }

    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Take a reference to the deformation gradient maintained by the passive law
    const mechanicalLaw& passiveLaw = passiveMechLawPtr_();
    const volTensorField& F = passiveLaw.F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    if (useFieldTa_)
    {
        // Field-based active tension from the coupling model
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        // Add active stress: convert 2nd Piola-Kirchhoff to Cauchy
        sigma += symm(F & (Ta*f0f0_) & F.T())/J;
    }
    else
    {
        // Constant active tension with optional ramp
        dimensionedScalar currentTa = Ta_;
        if (mesh().time().value() < rampTime_)
        {
            currentTa = (mesh().time().value()/rampTime_)*Ta_;
        }

        sigma += symm(F & (currentTa*f0f0_) & F.T())/J;
    }
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Lazy check for field-based active tension (same as vol variant)
    if (!fieldTaChecked_)
    {
        fieldTaChecked_ = true;
        useFieldTa_ = mesh().foundObject<volScalarField>("Ta");

        if (useFieldTa_)
        {
            Info<< "    electroMechanicalLaw: using field-based active tension"
                << " (Ta volScalarField from objectRegistry)" << endl;
        }
        else
        {
            Info<< "    electroMechanicalLaw: using constant active tension"
                << " Ta = " << Ta_.value()
                << " with rampTime = " << rampTime_ << endl;
        }
    }

    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    // Take a reference to the face deformation gradient maintained by the
    // passive law
    const mechanicalLaw& passiveLaw = passiveMechLawPtr_();
    const surfaceTensorField& F = passiveLaw.Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    if (useFieldTa_)
    {
        // Interpolate field-based active tension to faces
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        const surfaceScalarField Taf(fvc::interpolate(Ta));

        // Add active stress: convert 2nd Piola-Kirchhoff to Cauchy
        sigma += J*symm(F & (Taf*f0f0f_) & F.T());
    }
    else
    {
        // Constant active tension with optional ramp
        dimensionedScalar currentTa = Ta_;
        if (mesh().time().value() < rampTime_)
        {
            currentTa = (mesh().time().value()/rampTime_)*Ta_;
        }

        sigma += J*symm(F & (currentTa*f0f0f_) & F.T());
    }
}


// ************************************************************************* //
