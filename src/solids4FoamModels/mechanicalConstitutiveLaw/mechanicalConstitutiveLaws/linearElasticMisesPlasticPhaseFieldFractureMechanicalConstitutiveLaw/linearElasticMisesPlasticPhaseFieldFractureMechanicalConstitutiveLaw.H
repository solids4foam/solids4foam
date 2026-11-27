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
    Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw

Description
    Small-strain J2 plasticity with isotropic hardening coupled to the
    ductile phase field fracture model of Borden et al., as implemented by
    Whelan et al.

    The phase field d, 0 for intact material and 1 for a crack, is the
    solution of

        (1 + A) d - div((2 l)^2 grad(d)) = A,
        A = 4 l (1 - k)(H + <Wp - w0>)/Gc

    with zero flux on the boundary, where l is the characteristic length, Gc
    the critical fracture energy, H the history of the tensile elastic
    energy, Wp the plastic work and w0 its threshold. The solid model solves
    this equation once per outer iteration (it is the implicit-gradient
    field phaseField, of source A and coefficient 1 + A), and this law reads
    its solution as an input.

    The stress is degraded by g(d) = (1 - k)(1 - d)^2 + k, with k the
    residual stiffness, except for its volumetric part in compression:

        sigma = g 2 mu dev(epsilonE) + g_v K tr(epsilonE) I,
        g_v   = g if tr(epsilonE) >= 0, and 1 otherwise

    and the yield stress is degraded alike, so that the return mapping is
    that of linearElasticMisesPlastic. The tensile elastic energy, which
    drives the crack, is

        psi+ = mu dev(epsilonE):dev(epsilonE) + K <tr(epsilonE)>^2/2

    and its history is H = max(H0, psi+) over the step; the plastic work
    grows by sigmaY dEpsilonPEq, with the intact yield stress, as Whelan et
    al. write it, where Borden et al. use the degraded stress. Neither
    depends on d, so the crack driving force is a function of the
    kinematics and the old-time state alone.

    With a yield stress out of reach, the law is the brittle phase field
    model of Miehe et al., with the volumetric-deviatoric split of Amor et
    al.

    The law implements the small-strain update and the elastic predictor
    interface, so it is used at finite strain as the sub-law of logStrain,
    which is the formulation of Whelan et al.

    References:
      - M.J. Borden, T.J.R. Hughes, C.M. Landis, A. Anvari & I.J. Lee,
        Comput. Methods Appl. Mech. Eng. 312 (2016) 130-166.
      - C. Miehe, M. Hofacker & F. Welschinger, Comput. Methods Appl. Mech.
        Eng. 199 (2010) 2765-2778.
      - H. Amor, J.-J. Marigo & C. Maurini, J. Mech. Phys. Solids 57 (2009)
        1209-1229.
      - A. Whelan, T. Tang, V. Pakrashi & P. Cardiff, Int. J. Numer. Meth.
        Eng. (2025), doi:10.1002/nme.7640.

SourceFiles
    linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.
    Based on the phase field laws of Andrew Whelan, UCD.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw_H
#define linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw
                                  Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
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

        //- Critical fracture energy, per unit area
        scalar Gc_;

        //- Characteristic length
        scalar l_;

        //- Plastic work threshold
        scalar w0_;

        //- Residual stiffness of a crack
        const scalar k_;

        //- Newton loop tolerance
        const scalar loopTol_;

        //- Maximum Newton iterations
        const label maxNewtonIter_;

        //- Finite-difference step for the hardening slope
        const scalar finiteDiff_;


    // Private Member Functions

        //- Yield stress as a function of the equivalent plastic strain
        inline scalar yieldStress(const scalar epsilonPEq) const;

        //- Radial return from a trial elastic strain: the plastic multiplier
        //  increment, which the degradation does not change
        scalar plasticMultiplier
        (
            const scalar qTrial,
            const scalar epsilonPEq0,
            const scalar maxMagDEpsilon
        ) const;

        //- No copy construct
        linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw
        (
            const linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElasticMisesPlasticPhaseFieldFracture");


    // Static Data Members

        //- The name of the phase field
        static const word phaseFieldName;


    // Constructors

        //- Construct from dictionary
        explicit
        linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual
    ~linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw()
        = default;


    // Member Functions

    // Access

        //- Return the density
        virtual dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Return the bulk modulus of the intact material
        virtual dimensionedScalar kappa() const override
        {
            return kappa_;
        }

        //- The degradation function g(d)
        scalar degradation(const scalar d) const;

    // Edit

        //- Declare the plastic history and the crack driving force
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

        //- The phase field
        virtual wordList requiredScalarInputs() const override;

        //- The phase field, of source A and coefficient 1 + A
        virtual List<implicitGradientField>
        implicitGradientFields() const override;

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

        //- Count the yielding points of a step, and the largest crack
        //  driving force
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
