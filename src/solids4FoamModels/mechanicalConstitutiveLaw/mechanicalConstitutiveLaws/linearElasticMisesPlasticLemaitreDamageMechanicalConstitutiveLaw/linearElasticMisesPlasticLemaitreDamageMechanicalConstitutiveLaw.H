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
    Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw

Description
    Small-strain J2 plasticity with isotropic hardening coupled to Lemaitre
    ductile damage, local form.

    With D the damage, r the accumulated plastic multiplier (the hardening
    variable, stored as epsilonPEq) and sigma the stress, the yield function
    and flow rule are written on the effective stress sigma/(1 - D):

        Phi       = sigmaEq/(1 - D) - sigmaY(r)
        dEpsilonP = (3/2) dr/(1 - D) dev(sigma)/sigmaEq
        sigma     = (1 - D) C : epsilonE

    and the damage evolves with the energy release rate Y:

        dD = dr/(1 - D) (-Y/S)^b,   if r > epsilonD and eta > -1/3

    where eta = tr(sigma)/(3 sigmaEq) is the stress triaxiality. The
    triaxiality cut-off, on by default, stops damage under compressive
    triaxiality. S is the damage denominator, selected by lemaitreType:

        classic      S = S0
        Malcher      S = S0Tension/(3|eta| + (S0Tension/S0Shear)(1 - xi^2))
        triaxiality  S = 2 S0/(1 + 3 eta)

    with xi the normalised third invariant of the deviatoric stress. The
    last is the form proposed by Whelan et al. (2025) for wire drawing, and
    needs the cut-off, being singular at eta = -1/3.

    Y is the classic

        Y = -[(1 + nu) sigma:sigma - nu tr(sigma)^2]/(2 E (1 - D)^2)

    or, with a crack closure parameter h < 1, the form of Ladeveze and
    Lemaitre that weights the compressive principal stresses by h:

        Y = -[(1 + nu) sigma+:sigma+ - nu <tr(sigma)>^2]/(2 E (1 - D)^2)
          - h [(1 + nu) sigma-:sigma- - nu <-tr(sigma)>^2]/(2 E (1 - h D)^2)

    When the damage would pass the critical damage Dc, the point is failed:
    its damage is set to 0.99 and stays there.

    The return mapping and the damage update are solved together at each
    point, implicitly in time: for a trial damage the plastic multiplier is
    found by the radial return, and the damage is then the root of
    D - D0 - dD(D) = 0 in [D0, Dc]. This is the converged state of the
    staggered scheme of Whelan et al., which lags the damage by one outer
    iteration, but it is a function of the kinematics and the old-time state
    alone, as a law here must be.

    The law implements the small-strain update and the elastic predictor
    interface, so it is used at finite strain as the sub-law of logStrain,
    which is the formulation of Whelan et al.

    References:
      - J. Lemaitre, J. Eng. Mater. Technol. 107 (1985) 83-89.
      - E.A. de Souza Neto, D. Peric & D.R.J. Owen, Computational Methods
        for Plasticity, Wiley, 2008, chapter 12.
      - L. Malcher, F.M. Andrade Pires & J.M.A. Cesar de Sa, Int. J.
        Plasticity 30-31 (2012) 81-115.
      - A. Whelan, T. Tang, V. Pakrashi & P. Cardiff, Int. J. Numer. Meth.
        Eng. (2025), doi:10.1002/nme.7640.

SourceFiles
    linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.
    Based on the Lemaitre damage laws of Andrew Whelan, UCD.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw_H
#define linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw
                                  Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
public:

    // Public Types

        //- The form of the damage denominator S
        enum class lemaitreType
        {
            classic,
            Malcher,
            triaxiality
        };


    // Public Classes

        //- The result of the coupled return mapping at one point
        struct pointResult
        {
            //- Elastic strain at the end of the step
            symmTensor epsilonE;

            //- Plastic strain increment
            symmTensor DEpsilonP;

            //- Stress
            symmTensor sigma;

            //- Plastic multiplier increment, the hardening variable's
            scalar DLambda;

            //- Damage at the end of the step
            scalar D;

            //- Yield stress at the end of the step
            scalar sigmaY;

            //- Factor the return scaled the trial deviatoric stress by
            scalar theta;
        };


private:

    // Private Data

        //- Density
        const dimensionedScalar rho_;

        //- Young's modulus
        const dimensionedScalar E_;

        //- Poisson's ratio
        const dimensionedScalar nu_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;

        //- Table of plastic strain versus yield stress
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- True if the hardening is nonlinear (table of more than 2 points)
        const bool nonLinearPlasticity_;

        //- Linear hardening modulus (used only for linear hardening)
        scalar Hp_;

        //- Form of the damage denominator
        lemaitreType lemaitreType_;

        //- Damage denominator for classic and triaxiality
        scalar S0_;

        //- Malcher damage denominator in uniaxial tension, S_0.33
        scalar S0Tension_;

        //- Malcher damage denominator in pure shear, S_0.0
        scalar S0Shear_;

        //- Damage exponent
        const scalar b_;

        //- Accumulated plastic strain below which there is no damage
        const scalar epsilonD_;

        //- Critical damage, past which a point is failed
        const scalar Dc_;

        //- Crack closure parameter: 1 for none
        const scalar h_;

        //- Whether damage stops for triaxialities at or below -1/3
        const bool triaxialityCutOff_;

        //- Newton loop tolerance
        const scalar loopTol_;

        //- Maximum Newton iterations
        const label maxNewtonIter_;

        //- Finite-difference step for the hardening slope
        const scalar finiteDiff_;

        //- Absolute tolerance on the damage of the coupled solution
        const scalar damageTol_;

        //- Maximum iterations of the coupled solution
        const label maxDamageIter_;


    // Private Member Functions

        //- Yield stress as a function of the accumulated plastic multiplier
        inline scalar yieldStress(const scalar r) const;

        //- Radial return from a trial elastic strain, at a given damage
        void returnMap
        (
            const symmTensor& epsilonETrial,
            const scalar epsilonPEq0,
            const scalar D,
            const scalar maxMagDEpsilon,
            pointResult& result
        ) const;

        //- The energy release rate Y of a stress at a damage
        scalar energyReleaseRate
        (
            const symmTensor& sigma,
            const scalar D
        ) const;

        //- The damage increment of a returned step, at the damage the
        //  step was returned at
        scalar damageIncrement
        (
            const pointResult& result,
            const scalar epsilonPEq0
        ) const;

        //- The coupled return mapping and local damage update at one point
        void localUpdate
        (
            const symmTensor& epsilonETrial,
            const scalar epsilonPEq0,
            const scalar D0,
            const scalar maxMagDEpsilon,
            pointResult& result
        ) const;

        //- No copy construct
        linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw
        (
            const linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElasticMisesPlasticLemaitreDamage");


    // Static Data Members

        //- Damage of a failed point
        static const scalar Dfailed;


    // Constructors

        //- Construct from dictionary
        explicit linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw()
        = default;


    // Member Functions

    // Access

        //- Return the density
        virtual dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Return the bulk modulus of the undamaged material
        virtual dimensionedScalar kappa() const override
        {
            return kappa_;
        }

    // Edit

        //- Declare the plastic and damage history
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

        //- The largest strain magnitude over these points, which the Newton
        //  residual is normalised by
        virtual scalar smallStrainConvergenceScale
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawState& state
        ) const override;

        //- The convergence scale above is used
        virtual bool usesConvergenceScale() const override
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

        //- Isotropic, with scalar history, so this law can be driven by a
        //  trial elastic strain, as logStrain drives it
        virtual bool providesElasticPredictor() const override
        {
            return true;
        }

        //- Evaluation from a trial elastic strain. epsilonP is not used
        //  here and stays at its initial value
        virtual void evaluateElasticPredictor
        (
            const UList<symmTensor>& epsilonETrial,
            const mechanicalConstitutiveLawInputs& inputs,
            mechanicalConstitutiveLawState& state,
            UList<symmTensor>& epsilonE,
            mechanicalConstitutiveLawResponse& response
        ) const override;

        //- Count the yielding, damaging and failed points of a step
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
