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

#include "electroMechanicalLawMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroMechanicalLawMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        electroMechanicalLawMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroMechanicalLawMechanicalConstitutiveLaw::
electroMechanicalLawMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    subLawPtr_
    (
        mechanicalConstitutiveLaw::New
        (
            subLawDict(dict, "passiveMechanicalLaw")
        )
    ),
    subStateName_("passiveMechanicalLaw"),
    Ta_(dict.lookup("activeTension")),
    rampTime_(readScalar(dict.lookup("rampTime"))),
    TaName_(dict.lookupOrDefault<word>("activeTensionFieldName", "Ta")),
    f0Default_
    (
        dict.lookupOrDefault<Switch>("uniformFibreField", false)
      ? vector(dict.lookup("f0"))
      : vector::zero
    )
{
    if (rampTime_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "rampTime must not be negative." << exit(FatalIOError);
    }

    Info<< "    Active tension law over " << subLawPtr_->type()
        << ", Ta = " << Ta_.value() << ", rampTime = " << rampTime_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mechanicalConstitutiveLaw&
Foam::electroMechanicalLawMechanicalConstitutiveLaw::childLaw
(
    const word& name
) const
{
    if (name != subStateName_)
    {
        FatalErrorInFunction
            << "Mechanical constitutive law '" << type()
            << "' has no child law '" << name << "'."
            << exit(FatalError);
    }

    return subLawPtr_();
}


void Foam::electroMechanicalLawMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // The direction the tension pulls along. Declared here as well as by the
    // passive law below, because each reads it from its own state and the two
    // are separate namespaces - they are given the same field
    spec.addVector
    (
        "f0",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        f0Default_
    );
}


void Foam::electroMechanicalLawMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    // The passive stress first, written into the caller's storage
    subLawPtr_->evaluate(kin, inputs, state.child(subStateName_), response);

    UIndirectList<symmTensor>& sigma = response.stress();

    const UIndirectList<tensor>& F = kin.F();
    const UIndirectList<scalar>& J = kin.J();

    // Read at old time through a const reference: a prescribed field holds the
    // same value at both times, and the old one is what a shadow state aliases
    const mechanicalConstitutiveLawState& cState = state;
    const Field<vector>& f0 = cState.vectorField0("f0");

    // A field of active tension if the run has one, otherwise the constant,
    // ramped. The field wins, which is how a coupled electrophysiology model
    // overrides the standalone ramp without either knowing about the other
    const bool haveTaField = inputs.foundScalar(TaName_);

    scalar TaUniform = Ta_.value();

    if (!haveTaField && rampTime_ > 0 && inputs.time() < rampTime_)
    {
        TaUniform *= inputs.time()/rampTime_;
    }

    forAll(sigma, i)
    {
        const scalar magF0 = mag(f0[i]);

        if (magF0 < SMALL)
        {
            FatalErrorInFunction
                << "The fibre direction has zero length at index " << i
                << '.' << nl
                << "Active tension acts along the fibre, so this law needs "
                << "one. Either supply an f0 field - setFibreField writes one "
                << "- or set 'uniformFibreField yes' and give 'f0' in this "
                << "material's dictionary."
                << exit(FatalError);
        }

        const symmTensor f0f0(sqr(f0[i]/magF0));

        const scalar Ta =
            haveTaField ? inputs.getScalar(TaName_)[i] : TaUniform;

        // The active second Piola-Kirchhoff stress pushed forward to Cauchy
        sigma[i] += symm(F[i] & (Ta*f0f0) & F[i].T())/J[i];
    }
}


// ************************************************************************* //
