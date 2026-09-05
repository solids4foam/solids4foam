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

    The active tension Ta is either a field, when an electrophysiology model
    supplies one, or a constant ramped over rampTime. The field wins if it is
    there, which is how a coupled run overrides the standalone one without
    either knowing about the other.

    The passive law is given in a passiveMechanicalLaw sub-dictionary and owns
    a child state of its own, exactly as poroMechanicalLaw's effective stress
    law does.

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
        //  Declared only when the case asks for a field-based tension. The
        //  manager gathers these on the paths that support them and refuses on
        //  the paths that do not, which is why this is conditional: a law that
        //  declared it unconditionally would stop every finite-strain case
        //  from running, and one that never declared it would silently ignore
        //  the field it was given
        virtual wordList requiredScalarInputs() const override
        {
            if (TaFromField_)
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
        virtual bool providesVolumetricSplit() const override
        {
            return subLawPtr_->providesVolumetricSplit();
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
