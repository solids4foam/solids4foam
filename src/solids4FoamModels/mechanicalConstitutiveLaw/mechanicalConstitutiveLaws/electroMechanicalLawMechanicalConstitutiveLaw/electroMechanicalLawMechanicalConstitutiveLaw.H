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


Class
    Foam::electroMechanicalLawMechanicalConstitutiveLaw

Description
    Active contraction along a fibre direction, added to the stress of a
    passive law.

        sigma = sigma_passive + (1/J)*F & (Ta*f0 (x) f0) & F^T

    The active tension Ta is either a field, supplied by an electrophysiology
    model or read from the case, or a constant ramped over rampTime. The field
    is read when activeTensionFromField is set, under the name given by
    activeTensionFieldName (default Ta), and the constant is used otherwise.
    Unlike the removed legacy law, a registered field is not picked up
    unasked: a coupled case sets activeTensionFromField, and one that does not
    is warned that it is running on the constant.

    The passive law is given in a passiveMechanicalLaw sub-dictionary and owns
    a child state of its own, exactly as poroMechanicalLaw's effective stress
    law does.

    Tangents. The active stress is added after the passive law has run, so a
    tangent that law wrote describes the passive response only. A
    fourthOrderFiniteDifference request is therefore recomputed here over the
    whole stress, and the analytical fourth-order request is refused rather
    than answered with the passive part of the answer.

    The scalar tangent is the passive law's, and the active tension is not in
    it. That one is a stiffness estimate for the segregated solver rather than
    a quantity the answer depends on, so it costs iterations rather than
    accuracy, and there is no finite-difference scalar tangent in the base
    class to recompute it with. Said here rather than silently.

SourceFiles
    electroMechanicalLawMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef electroMechanicalLawMechanicalConstitutiveLaw_H
#define electroMechanicalLawMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class electroMechanicalLawMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class electroMechanicalLawMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- The passive law this one adds tension to
        autoPtr<mechanicalConstitutiveLaw> subLawPtr_;

        //- Name of the child state the passive law keeps its history in
        const word subStateName_;

        //- Active tension, when it is not supplied as a field
        const dimensionedScalar Ta_;

        //- Time over which the constant tension is ramped from zero
        const scalar rampTime_;

        //- Name of the active tension field, when one is supplied
        const word TaName_;

        //- Whether the active tension is to come from a field.
        //  Asked for explicitly rather than detected, because a law that
        //  quietly looks for a field it was never given the chance to receive
        //  reports nothing and uses the constant instead - which is a wrong
        //  answer that runs to completion
        const Switch TaFromField_;

        //- The fibre direction to use where a case supplies no f0 field
        const vector f0Default_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        electroMechanicalLawMechanicalConstitutiveLaw
        (
            const electroMechanicalLawMechanicalConstitutiveLaw&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=
        (
            const electroMechanicalLawMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("electroMechanicalLaw");

    // Constructors

        //- Construct from dictionary
        explicit electroMechanicalLawMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~electroMechanicalLawMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- Return the density, which is the passive law's
        virtual dimensionedScalar rho() const override
        {
            return subLawPtr_->rho();
        }

        //- Return the bulk modulus, which is the passive law's
        virtual dimensionedScalar kappa() const override
        {
            return subLawPtr_->kappa();
        }

        //- The coupling inputs this law reads.
        //
        //  Declared only when the case asks for a field-based tension, since
        //  the manager stops a run whose declared input is neither registered
        //  nor given in the case, and a standalone run has neither
        virtual wordList requiredScalarInputs() const override
        {
            if (TaFromField_)
            {
                return wordList(1, TaName_);
            }

            return wordList();
        }

        //- The active tension field, when the constant is used in its place
        virtual wordList optionalScalarInputs() const override
        {
            if (!TaFromField_)
            {
                return wordList(1, TaName_);
            }

            return wordList();
        }

        //- Whether the two parts can be separated, which is the passive
        //  law's answer.
        //
        //  This composite adds an active tension whose mean stress is not
        //  zero, and that mean is physical: a fibre carrying tension along
        //  itself and nothing across it has one. It is also not a function of
        //  J, so it is not dU/dJ and does not belong in the volumetric
        //  channel. It stays with the stress, which is the part a solved
        //  pressure must not overwrite, so this composite adds nothing to the
        //  volumetric response and its answer is the passive law's
        //- A composite needs whatever its sub-law needs. Without this the
        //  manager never attaches the standing stress for, say, a poro law
        //  nested inside this one, and that law fails at the point of use
        virtual bool requiresIncomingStress() const override
        {
            return subLawPtr_->requiresIncomingStress();
        }

        virtual bool providesVolumetricSplit() const override
        {
            return subLawPtr_->providesVolumetricSplit();
        }

        //- The active tension is added as a Cauchy stress
        //  symm(F & (Ta*f0f0) & F.T())/J, which is not derived from a
        //  potential and scales with a superposed dilation. So the split this
        //  law returns is honest but not dilation invariant, whatever the
        //  passive law underneath does. It happens to be invariant while Ta
        //  is zero, which is not a property worth claiming
        virtual bool isochoricStressIsDilationInvariant() const override
        {
            return false;
        }

        //- The passive law's convergence scale, over its own state, so that a
        //  plastic passive law is judged by the scale the manager reduces
        virtual scalar localConvergenceScale
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawState& state
        ) const override
        {
            return subLawPtr_->localConvergenceScale
            (
                kin, state.child(subStateName_)
            );
        }

        //- The child state the passive law owns
        virtual wordList childStateNames() const override
        {
            return wordList(1, subStateName_);
        }

        //- The passive law that uses the named child state
        virtual const mechanicalConstitutiveLaw& childLaw
        (
            const word& name
        ) const override;

        //- Declare the fibre direction the active tension acts along
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

        //- Finite-strain constitutive evaluation
        virtual void evaluate
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
