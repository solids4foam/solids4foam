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
    Foam::poroMechanicalLawMechanicalConstitutiveLaw

Description
    Total stress of a porous medium, from the effective stress of a sub-law and
    a pore pressure.

    This is a composite. It owns a sub-law which supplies the effective stress,
    and returns

        sigma = sigmaEff - b*(p + p0)*I

    where b is the Biot coefficient, p the pore pressure and p0 a reference
    pressure.

    It is not simply "the sub-law's stress minus the pressure". The effective
    stress is a field of its own, seeded once from the incoming total stress as
    sigmaEff = sigma + b*(p + p0)*I, and the sub-law is given that field to
    work in. For a sub-law that computes its stress purely from the kinematics
    and its own history this makes no difference to the answer. For one whose
    strength depends on the stress state - Mohr-Coulomb is the case the legacy
    law calls out - it is the difference between the sub-law seeing effective
    stress and seeing total stress, so it is reproduced rather than simplified
    away.

    The pore pressure is a live coupling input, gathered by the manager on
    every evaluation. It is not state: it is solved for elsewhere.

SourceFiles
    poroMechanicalLawMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef poroMechanicalLawMechanicalConstitutiveLaw_H
#define poroMechanicalLawMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class poroMechanicalLawMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class poroMechanicalLawMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- The law that supplies the effective stress
        autoPtr<mechanicalConstitutiveLaw> subLawPtr_;

        //- Name of the child state the sub-law's history lives in
        const word subStateName_;

        //- Biot coefficient
        dimensionedScalar b_;

        //- Reference pore pressure
        dimensionedScalar p0_;

        //- Name of the pore pressure field
        const word pName_;


    // Private Member Functions


        //- No copy construct
        poroMechanicalLawMechanicalConstitutiveLaw
        (
            const poroMechanicalLawMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const poroMechanicalLawMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("poroMechanicalLaw");

    // Constructors

        //- Construct from dictionary
        explicit poroMechanicalLawMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~poroMechanicalLawMechanicalConstitutiveLaw() = default;


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

        //- The pore pressure, gathered on every evaluation
        virtual wordList requiredScalarInputs() const override
        {
            return wordList(1, pName_);
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

        //- The effective stress the sub-law works in, and whether it has been
        //  seeded yet
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
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
