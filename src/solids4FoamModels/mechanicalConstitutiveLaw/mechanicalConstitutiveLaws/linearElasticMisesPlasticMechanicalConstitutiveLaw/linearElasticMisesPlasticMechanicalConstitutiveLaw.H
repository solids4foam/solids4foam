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
    Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw

Description
    Small-strain elasto-plastic constitutive law with linear (Hookean)
    elasticity and isotropic J2 (von Mises) plasticity.

    The formulation follows Simo & Hughes (1998), Box 3.2, using a radial
    return mapping algorithm. Both perfect plasticity and isotropic
    hardening are supported via a user-specified yield-stress versus
    equivalent plastic strain curve.

    This class implements only the constitutive response at a collection
    of integration points. Storage and management of history variables
    (plastic strain, yield stress, etc.) are handled via
    mechanicalConstitutiveLawState objects owned by the
    mechanicalConstitutiveLawManager.

    The law is agnostic to spatial discretisation (cell-centred,
    vertex-centred, higher-order integration points, etc.) and operates
    exclusively on indirect views provided by the manager.

    Elastic parameters may be specified as:
      - Young's modulus (E) and Poisson's ratio (nu), or
      - Shear modulus (mu) and bulk modulus (kappa).

    References:
      - J.C. Simo & T.J.R. Hughes, Computational Inelasticity, Springer, 1998.
      - P. Cardiff et al., Int. J. Numer. Meth. Eng., doi:10.1002/nme.5345.

SourceFiles
    linearElasticMisesPlasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMisesPlasticMechanicalConstitutiveLaw_H
#define linearElasticMisesPlasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class linearElasticMisesPlasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMisesPlasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private Data

        //- Density
        const dimensionedScalar rho_;

        //- Young's modulus
        //  Set in the constructor body, since the material may be given as
        //  mu and K instead
        dimensionedScalar E_;

        //- Poisson's ration
        dimensionedScalar nu_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;

        //- Table of plastic strain versus post-yield stress
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- True if plasticity law is nonlinear (table > 2 points)
        const bool nonLinearPlasticity_;

        //- Linear hardening modulus (used only for linear plasticity)
        scalar Hp_;

        //- Newton loop tolerance
        const scalar loopTol_;

        //- Maximum Newton iterations
        const label maxNewtonIter_;

        //- Finite-difference step for numerical derivatives
        const scalar finiteDiff_;

    // Private Member Functions

        //- Current yield stress as a function of equivalent plastic strain
        inline scalar yieldStress
        (
            const scalar epsilonPEq
        ) const;

        //- Yield function value
        inline scalar yieldFunction
        (
            const scalar epsilonPEqOld,
            const scalar magSTrial,
            const scalar dLambda,
            const scalar muBar
        ) const;

        //- Newton iteration for plastic multiplier
        inline void newtonLoop
        (
            scalar& dLambda,
            scalar& curSigmaY,
            const scalar epsilonPEqOld,
            const scalar magSTrial,
            const scalar muBar,
            const scalar maxMagDEpsilon
        ) const;

        //- No copy construct
        linearElasticMisesPlasticMechanicalConstitutiveLaw
        (
            const linearElasticMisesPlasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticMisesPlasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElasticMisesPlastic");

    // Constructors

        //- Construct from dictionary
        explicit linearElasticMisesPlasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearElasticMisesPlasticMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- Return the density
        virtual dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Return the bulk modulus
        virtual dimensionedScalar kappa() const override
        {
            return kappa_;
        }

    // Edit

        //- Declare this law's history, so that the manager can carry it
        //  across a restart
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

        //- The stress is a deviatoric part plus K*tr(epsilon), and the
        //  return mapping only touches the deviatoric part, so the two can be
        //  handed back separately
        virtual bool providesVolumetricSplit() const override
        {
            return true;
        }

        //- Small-strain constitutive evaluation
        virtual void evaluate
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;

        //- Optional hook called once at the end of each time step.
        virtual void endTimeStep
        (
            const mechanicalConstitutiveLawState& state,
            const scalar time,
            const label timeIndex,
            DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics
        ) const override;

        //- Say what was found, once the manager has reduced it
        virtual void reportDiagnostics
        (
            const UList<mechanicalConstitutiveLawDiagnostic>& diagnostics
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
