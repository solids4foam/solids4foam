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

#include "poroMechanicalLawMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(poroMechanicalLawMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        poroMechanicalLawMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poroMechanicalLawMechanicalConstitutiveLaw::
poroMechanicalLawMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    subLawPtr_
    (
        mechanicalConstitutiveLaw::New
        (
            subLawDict(dict, "effectiveStressMechanicalLaw")
        )
    ),
    subStateName_("effectiveStressMechanicalLaw"),
    b_
    (
        dict.lookupOrDefault<dimensionedScalar>
        (
            "biotCoeff", dimensionedScalar("b", dimless, 1.0)
        )
    ),
    p0_
    (
        dict.lookupOrDefault<dimensionedScalar>
        (
            "p0", dimensionedScalar("p0", dimPressure, 0.0)
        )
    ),
    pName_(dict.lookupOrDefault<word>("pressureFieldName", "porePressure"))
{
    if (b_.dimensions() != dimless)
    {
        FatalIOErrorInFunction(dict)
            << "The Biot coefficient biotCoeff must be dimensionless. Got "
            << b_.dimensions()
            << exit(FatalIOError);
    }

    if (p0_.dimensions() != dimPressure)
    {
        FatalIOErrorInFunction(dict)
            << "The reference pressure p0 must have dimensions " << dimPressure
            << ". Got " << p0_.dimensions()
            << exit(FatalIOError);
    }

    Info<< "    Pore pressure '" << pName_ << "' on top of "
        << subLawPtr_->type() << ", biotCoeff = " << b_.value()
        << ", p0 = " << p0_.value() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::mechanicalConstitutiveLaw&
Foam::poroMechanicalLawMechanicalConstitutiveLaw::childLaw
(
    const word& name
) const
{
    if (name != subStateName_)
    {
        FatalErrorInFunction
            << "Unknown child law '" << name << "'. This law owns only '"
            << subStateName_ << "'."
            << exit(FatalError);
    }

    return subLawPtr_();
}


void Foam::poroMechanicalLawMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // The effective stress the sub-law works in. Persistent, because it is
    // seeded once and then carried
    spec.addSymmTensor
    (
        "sigmaEff",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        symmTensor::zero
    );

    // Whether the seeding has happened. Kept per point rather than as a flag
    // on the law, because a law has one state per material and one per
    // boundary patch, and each is seeded when it is first evaluated
    spec.addScalar
    (
        "sigmaEffSeeded",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        0.0
    );
}


void Foam::poroMechanicalLawMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const scalarField& p = inputs.getScalar(pName_);

    Field<symmTensor>& sigmaEff = state.symmTensorField("sigmaEff");
    Field<scalar>& seeded = state.scalarField("sigmaEffSeeded");

    const scalar b = b_.value();
    const scalar p0 = p0_.value();

    // Seed the effective stress from the total stress standing at these
    // points, once per point. That stress arrives as a declared input rather
    // than being read back out of the response: a law that reads its own
    // output buffer depends on whatever the caller last left there, which is
    // not something that can be restarted, perturbed for a tangent, or
    // reasoned about. This is what makes the sub-law's first evaluation see an
    // effective stress rather than whatever happened to be in the field
    forAll(sigma, i)
    {
        if (seeded[i] < 0.5)
        {
            sigmaEff[i] =
                inputs.incomingStress()[i] + (b*(p[i] + p0))*symmTensor(I);
            seeded[i] = 1.0;
        }
    }

    // Hand the sub-law the effective stress to work in. It is passed through
    // the caller's storage rather than through a list of its own: the response
    // wraps an indirect list, and copying in and out costs two passes where
    // building an identity-indexed view would cost an allocation. The sub-law
    // sees the effective stress it left behind, which is what the legacy law
    // gives it
    forAll(sigma, i)
    {
        sigma[i] = sigmaEff[i];
    }

    // What the sub-law finds already standing at these points is the
    // effective stress, not the total. It used to learn that by reading the
    // caller's storage, which this law had just written sigmaEff into; now it
    // is said explicitly, because a law reading its own output buffer is the
    // thing being removed
    inputs.setIncomingStress(sigmaEff);

    subLawPtr_->evaluate(kin, inputs, state.child(subStateName_), response);

    // What the sub-law wrote is the new effective stress. Keep it, and hand
    // the caller the total stress
    forAll(sigma, i)
    {
        sigmaEff[i] = sigma[i];
        sigma[i] -= (b*(p[i] + p0))*symmTensor(I);
    }
}


// ************************************************************************* //
