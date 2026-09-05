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
#include "mechanicalConstitutiveLawInputs.H"
#include "mechanicalConstitutiveLawStateSpec.H"
#include "mechanicalConstitutiveLawResponse.H"
#include "dimensionedScalar.H"
#include "dimensionSets.H"


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


protected:

    //- Finite-difference fourth-order tangent (small strain)
    //  Finite-difference tangent evaluation is guaranteed to leave
    //  the supplied constitutive state unmodified.
    void finiteDifferenceFourthOrder
    (
        const smallStrainMechanicalConstitutiveLawKinematics& kin,
        const mechanicalConstitutiveLawInputs& inputs,
        mechanicalConstitutiveLawState& state,
        mechanicalConstitutiveLawResponse& response
    ) const;

        //- Finite-strain counterpart of the above.
        //  Perturbs the deformation gradient rather than the displacement
        //  gradient, recomputing its inverse and determinant for each
        //  perturbed state, and differences the resulting Cauchy stress.
        //  The perturbation is defined on the same Voigt strain vector, so a
        //  tangent taken here is directly comparable with a small-strain one:
        //  at F close to the identity a hyperelastic law should return its
        //  linearised isotropic tangent
        void finiteDifferenceFourthOrder
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const;


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

            // Keep the compiler happy.
            // Note: dimensioned<Type> has no default constructor in
            // foam-extend, so the dummy must be constructed explicitly
            return dimensionedScalar("rho", dimDensity, 0.0);
        }

        //- Return the bulk modulus for this constitutive law
        //  Default: not specified
        virtual dimensionedScalar kappa() const
        {
            FatalErrorInFunction
                << "Bulk modulus (kapp) not defined for mechanical constitutive "
                << "law of type " << type() << exit(FatalError);

            // Keep the compiler happy.
            // Note: dimensioned<Type> has no default constructor in
            // foam-extend, so the dummy must be constructed explicitly
            return dimensionedScalar("kappa", dimPressure, 0.0);
        }

    // Edit

        //- Declare the state fields this law uses.
        //  Called once, before any evaluation. The manager allocates them at
        //  the size of the topology, applies the declared defaults, and for a
        //  prescribed field looks for a volField of the same name on the mesh
        //  and broadcasts it to this law's integration points.
        //
        //  A law that declares nothing keeps working exactly as before: this
        //  is additive, and a prescribed field that no case supplies takes its
        //  declared default
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const
        {}

        //- Names of the live coupling inputs this law reads, such as "T" for a
        //  temperature. The manager gathers each from the field of that name
        //  and hands it over in the inputs object, as a view of this law's own
        //  integration points.
        //
        //  These are not state. A coupling input is solved for by something
        //  else and changes every iteration, so it is gathered afresh on every
        //  evaluation: never stored, never rolled over to old time, and never
        //  read from a file that could shadow the live field
        virtual wordList requiredScalarInputs() const
        {
            return wordList();
        }

        //- A sub-law's dictionary, carrying the entries it needs from the
        //  composite's own.
        //
        //  Two kinds of entry live on the parent and have to reach the child.
        //  A case gives rho once, on the composite, because it describes one
        //  material; the legacy laws get away with not passing it down because
        //  they read it only when asked and nothing asks a sub-law, where a
        //  framework law reads its material properties when it is constructed.
        //  And the manager injects planeStress and solutionD into each
        //  top-level law's dictionary, which is not where a sub-law looks.
        //
        //  Shared so that the composites cannot drift apart on which entries
        //  those are
        static dictionary subLawDict
        (
            const dictionary& dict,
            const word& subDictName
        );

        //- Names of the child states this law owns, one per sub-law.
        //
        //  A composite law delegates to laws of its own, and each of those has
        //  history it looks up by its own unqualified name. Naming the
        //  children here is what lets the manager reach them: it applies each
        //  sub-law's declared state to the matching child, recursively, so a
        //  sub-law's prescribed fields and defaults work exactly as they would
        //  if it were the top-level law
        virtual wordList childStateNames() const
        {
            return wordList();
        }

        //- The sub-law that uses the named child state
        virtual const mechanicalConstitutiveLaw& childLaw
        (
            const word& name
        ) const
        {
            FatalErrorInFunction
                << "Mechanical constitutive law '" << type()
                << "' has no child law '" << name << "'."
                << exit(FatalError);

            return *this;
        }

        //- Can this law separate its isochoric and volumetric responses?
        //
        //  A mixed displacement-pressure formulation replaces the volumetric
        //  part of the stress with a solved pressure, and to do that correctly
        //  it needs the law's isochoric stress rather than the trace-free part
        //  of its total. The two coincide only when the law's energy is
        //  written on an isochoric measure; for a law written on the full
        //  deformation they differ, and the difference is a different material
        //  rather than a visible error.
        //
        //  False by default, which is the safe direction: a law that has not
        //  considered the question does not silently qualify. A solid model
        //  running a mixed formulation refuses to start with a law that
        //  answers false, naming it, rather than taking a trace and hoping
        virtual bool providesVolumetricSplit() const
        {
            return false;
        }

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
            const mechanicalConstitutiveLawInputs& inputs,
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
            const mechanicalConstitutiveLawInputs&,
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
