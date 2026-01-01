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
    Foam::mechanicalConstitutiveLaw

Description
    Abstract base class for mechanical constitutive laws in solids4foam.

    A mechanicalConstitutiveLaw is a stateless constitutive operator that
    computes stress (and optionally tangent-related quantities) from the
    supplied kinematic fields and constitutive state.

    All history-dependent data required by a constitutive law is stored
    externally in one or more mechanicalConstitutiveLawState objects,
    owned and managed by the mechanicalConstitutiveLawManager. The law
    itself does not own or persist state variables between evaluations.

    Constitutive evaluation is performed via the evaluate() interface,
    which accepts a kinematics container, a state object, and a response
    container. The response container provides non-owning views onto the
    output fields to be populated by the law.

    This base class defines the interface only. Concrete constitutive
    laws (e.g. linear elasticity, plasticity, hyperelasticity) implement
    one or more evaluate() overloads corresponding to the kinematic
    formulations they support.

    Notes:
    - The law interface is independent of OpenFOAM GeometricField types.
      Field extraction/slicing and discretisation-specific concerns are
      handled exclusively by the mechanicalConstitutiveLawManager.
    - Tangent-related outputs are optional and may be requested by the
      calling solver via the response container.

SourceFiles
    mechanicalConstitutiveLaw.C
    mechanicalConstitutiveLawNew.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLaw_H
#define mechanicalConstitutiveLaw_H

#include "dictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "smallStrainMechanicalConstitutiveLawKinematics.H"
#include "finiteStrainMechanicalConstitutiveLawKinematics.H"
#include "mechanicalConstitutiveLawState.H"
#include "mechanicalConstitutiveLawResponse.H"
#include "dimensionedScalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLaw
{
    // Private Member Functions

        //- Disallow copy construct
        mechanicalConstitutiveLaw(const mechanicalConstitutiveLaw&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalConstitutiveLaw&);


public:

    //- Runtime type information
    TypeName("mechanicalConstitutiveLaw");


    // Declare run-time constructor selection table

        //- Declare run-time selection table for the mechanical law
        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalConstitutiveLaw,
            mechanicalConstitutiveLaw,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Create a mechanical law
        static autoPtr<mechanicalConstitutiveLaw> New
        (
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        explicit mechanicalConstitutiveLaw(const dictionary& dict);


    // Destructor

        virtual ~mechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- Return the density for this constitutive law
        //  Default: not specified
        virtual dimensionedScalar rho() const
        {
            FatalErrorInFunction
                << "Density (rho) not defined for mechanical constitutive "
                << "law of type " << type() << exit(FatalError);

            // Keep the compiler happy
            return dimensionedScalar();
        }

        //- Return the bulk modulus for this constitutive law
        //  Default: not specified
        virtual dimensionedScalar kappa() const
        {
            FatalErrorInFunction
                << "Bulk modulus (kapp) not defined for mechanical constitutive "
                << "law of type " << type() << exit(FatalError);

            // Keep the compiler happy
            return dimensionedScalar();
        }

    // Edit

        //- Initialise constitutive state variables
        //  Called once after state allocation.
        //  Default implementation does nothing.
        virtual void initialiseState
        (
            mechanicalConstitutiveLawState& state
        ) const
        {
            // Do nothing by default
        }

        //- Small-strain constitutive evaluation
        //  Default implementation indicates unsupported formulation.
        virtual void evaluate
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const
        {
            FatalErrorInFunction
                << "Small-strain evaluation is not implemented for this "
                << "mechanicalConstitutiveLaw." << exit(FatalError);
        }

        //- Finite-strain constitutive evaluation (optional for a given law)
        //  Default implementation indicates unsupported formulation.
        virtual void evaluate
        (
            const finiteStrainMechanicalConstitutiveLawKinematics&,
            mechanicalConstitutiveLawState&,
            mechanicalConstitutiveLawResponse&
        ) const
        {
            FatalErrorInFunction
                << "Finite-strain evaluation is not implemented for this "
                << "mechanicalConstitutiveLaw." << exit(FatalError);
        }

        //- Optional hook called once at the end of each time step.
        //  Default implementation does nothing.
        //  Can be used for diagnostics, logging, or state bookkeeping.
        virtual void endTimeStep
        (
            const mechanicalConstitutiveLawState& state,
            const scalar time,
            const label timeIndex
        ) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
