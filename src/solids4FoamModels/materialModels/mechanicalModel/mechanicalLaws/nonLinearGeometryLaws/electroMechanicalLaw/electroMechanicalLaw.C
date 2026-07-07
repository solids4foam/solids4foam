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


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::electroMechanicalLaw::checkFieldTa() const
{
    if (!fieldTaChecked_)
    {
        fieldTaChecked_ = true;
        useFieldTa_ = mesh().foundObject<volScalarField>("Ta");

        if (useFieldTa_)
        {
            Info<< "    electroMechanicalLaw: using field-based active tension"
                << " (Ta volScalarField from objectRegistry)" << endl;
        }
        else if (mag(Ta_.value()) > SMALL)
        {
            Info<< "    electroMechanicalLaw: using constant active tension"
                << " Ta = " << Ta_.value()
                << " with rampTime = " << rampTime_ << endl;
        }
    }
}


Foam::dimensionedScalar Foam::electroMechanicalLaw::activeTension() const
{
    dimensionedScalar currentTa = Ta_;
    if (mesh().time().value() < rampTime_)
    {
        currentTa = (mesh().time().value()/rampTime_)*Ta_;
    }

    return currentTa;
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


bool Foam::electroMechanicalLaw::hasActiveStress() const
{
    checkFieldTa();

    if (useFieldTa_)
    {
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        return max(mag(Ta)).value() > SMALL;
    }

    return mag(activeTension().value()) > SMALL;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::electroMechanicalLaw::activeCauchyStress() const
{
    checkFieldTa();

    // Take a reference to the deformation gradient
    const volTensorField& F =
        const_cast<electroMechanicalLaw&>(*this).F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    if (useFieldTa_)
    {
        // Field-based active tension from the coupling model
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
        return tmp<volSymmTensorField>
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "sigmaActive",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                symm(F & (Ta*f0f0_) & F.T())/J
            )
        );
    }

    // Constant active tension with optional ramp
    const dimensionedScalar currentTa = activeTension();

    // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "sigmaActive",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            symm(F & (currentTa*f0f0_) & F.T())/J
        )
    );
}


Foam::tmp<Foam::surfaceSymmTensorField>
Foam::electroMechanicalLaw::activeCauchyStressf() const
{
    checkFieldTa();

    // Take a reference to the deformation gradient
    const surfaceTensorField& F =
        const_cast<electroMechanicalLaw&>(*this).Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    if (useFieldTa_)
    {
        // Interpolate field-based active tension to faces
        const volScalarField& Ta =
            mesh().lookupObject<volScalarField>("Ta");

        const surfaceScalarField Taf(fvc::interpolate(Ta));

        // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
        return tmp<surfaceSymmTensorField>
        (
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "sigmaActivef",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                (1.0/J)*symm(F & (Taf*f0f0f_) & F.T())
            )
        );
    }

    // Constant active tension with optional ramp
    const dimensionedScalar currentTa = activeTension();

    // Convert active 2nd Piola-Kirchhoff stress to Cauchy stress
    return tmp<surfaceSymmTensorField>
    (
        new surfaceSymmTensorField
        (
            IOobject
            (
                "sigmaActivef",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (1.0/J)*symm(F & (currentTa*f0f0f_) & F.T())
        )
    );
}


void Foam::electroMechanicalLaw::correct(volSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    if (hasActiveStress())
    {
        const tmp<volSymmTensorField> tSigmaActive = activeCauchyStress();
        sigma += tSigmaActive();
    }
}


void Foam::electroMechanicalLaw::correct(surfaceSymmTensorField& sigma)
{
    // Calculate passive stress
    passiveMechLawPtr_->correct(sigma);

    if (hasActiveStress())
    {
        const tmp<surfaceSymmTensorField> tSigmaActive = activeCauchyStressf();
        sigma += tSigmaActive();
    }
}


// ************************************************************************* //
