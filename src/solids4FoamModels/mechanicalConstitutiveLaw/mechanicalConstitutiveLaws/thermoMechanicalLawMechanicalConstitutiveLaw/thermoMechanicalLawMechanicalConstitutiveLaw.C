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

#include "thermoMechanicalLawMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoMechanicalLawMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        thermoMechanicalLawMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoMechanicalLawMechanicalConstitutiveLaw::
thermoMechanicalLawMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    subLawPtr_
    (
        mechanicalConstitutiveLaw::New(subLawDict(dict, "mechanicalLaw"))
    ),
    subStateName_("mechanicalLaw"),
    alpha_(dict.lookup("alpha")),
    T0_(dict.lookup("T0")),
    TName_(dict.lookupOrDefault<word>("TName", "T"))
{
    if (alpha_.dimensions() != dimless/dimTemperature)
    {
        FatalIOErrorInFunction(dict)
            << "The thermal expansion coefficient alpha must have dimensions "
            << dimless/dimTemperature << ". Got " << alpha_.dimensions()
            << exit(FatalIOError);
    }

    if (T0_.dimensions() != dimTemperature)
    {
        FatalIOErrorInFunction(dict)
            << "The reference temperature T0 must have dimensions "
            << dimTemperature << ". Got " << T0_.dimensions()
            << exit(FatalIOError);
    }

    Info<< "    Thermal stress on top of " << subLawPtr_->type()
        << ", alpha = " << alpha_.value()
        << ", T0 = " << T0_.value() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::mechanicalConstitutiveLaw&
Foam::thermoMechanicalLawMechanicalConstitutiveLaw::childLaw
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


void Foam::thermoMechanicalLawMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    // The mechanical response, written straight into the caller's storage.
    // The sub-law gets its own child state, so its history is its own, and the
    // response is a view rather than a copy, so there is nothing to write back
    subLawPtr_->evaluate(kin, inputs, state.child(subStateName_), response);

    UIndirectList<symmTensor>& sigma = response.stress();

    const scalarField& T = inputs.getScalar(TName_);

    // The thermal term uses the sub-law's bulk modulus, which is what the
    // legacy law does. Note that it is taken as a single value: a sub-law
    // whose bulk modulus varies from point to point cannot be expressed
    // through this interface, which is recorded in DESIGN-state-io.md
    // section 8a.2
    const scalar threeKAlpha =
        3.0*subLawPtr_->kappa().value()*alpha_.value();

    const scalar T0 = T0_.value();

    forAll(sigma, i)
    {
        sigma[i] -= (threeKAlpha*(T[i] - T0))*symmTensor(I);
    }
}


// ************************************************************************* //
