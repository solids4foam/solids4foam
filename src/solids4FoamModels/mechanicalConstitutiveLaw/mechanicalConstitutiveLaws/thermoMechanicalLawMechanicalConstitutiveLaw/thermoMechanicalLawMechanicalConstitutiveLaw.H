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
    Foam::thermoMechanicalLawMechanicalConstitutiveLaw

Description
    Thermal stress added to the stress of another law.

    This is a composite: it owns a sub-law, delegates the mechanical response
    to it, and subtracts the thermal term

        sigma -= 3*K*alpha*(T - T0)*I

    where K is the sub-law's bulk modulus, alpha the linear expansion
    coefficient and T0 the stress-free reference temperature.

    The temperature is a live coupling input, declared through
    requiredScalarInputs() and gathered by the manager on every evaluation. It
    is deliberately not state: it is solved for by the solid model and changes
    every iteration, so storing it, rolling it over, or letting a file of the
    same name stand in for it would all be wrong.

    The sub-law's history lives in a child state, so the sub-law looks it up by
    its own unqualified names and nothing collides with this law's.

SourceFiles
    thermoMechanicalLawMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef thermoMechanicalLawMechanicalConstitutiveLaw_H
#define thermoMechanicalLawMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class thermoMechanicalLawMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class thermoMechanicalLawMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- The law the mechanical response is delegated to
        autoPtr<mechanicalConstitutiveLaw> subLawPtr_;

        //- Name of the child state the sub-law's history lives in
        const word subStateName_;

        //- Linear thermal expansion coefficient
        dimensionedScalar alpha_;

        //- Stress-free reference temperature
        dimensionedScalar T0_;

        //- Name of the temperature field
        const word TName_;


    // Private Member Functions


        //- No copy construct
        thermoMechanicalLawMechanicalConstitutiveLaw
        (
            const thermoMechanicalLawMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const thermoMechanicalLawMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("thermoMechanicalLaw");

    // Constructors

        //- Construct from dictionary
        explicit thermoMechanicalLawMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermoMechanicalLawMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- The density is the sub-law's: this law adds a stress, not a
        //  material
        virtual dimensionedScalar rho() const override
        {
            return subLawPtr_->rho();
        }

        //- As is the bulk modulus
        virtual dimensionedScalar kappa() const override
        {
            return subLawPtr_->kappa();
        }

        //- The temperature, gathered by the manager on every evaluation
        virtual wordList requiredScalarInputs() const override
        {
            return wordList(1, TName_);
        }

        //- The sub-law's history lives in a child state
        virtual wordList childStateNames() const override
        {
            return wordList(1, subStateName_);
        }

        //- The sub-law that uses the named child state
        virtual const mechanicalConstitutiveLaw& childLaw
        (
            const word& name
        ) const override;

    // Edit

        //- Small-strain constitutive evaluation
        virtual void evaluate
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
