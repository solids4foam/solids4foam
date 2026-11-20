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
    Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw

Description
    Small-strain Gurson-Tvergaard-Needleman (GTN) plasticity: a porous
    plastic material whose void volume fraction f grows with the plastic
    dilatation and by nucleation, in local or non-local form.

    With p = tr(sigma)/3 the mean stress, q the von Mises equivalent stress
    and sigmaY(epsilonPEq) the yield stress of the matrix, the yield
    function is

        Phi = (q/sigmaY)^2 + 2 q1 f* cosh(3 q2 p/(2 sigmaY)) - 1 - q3 f*^2

    with associated flow, dEpsilonP = dLambda dPhi/dsigma. The matrix
    equivalent plastic strain follows from the equivalence of plastic work,

        (1 - f) sigmaY dEpsilonPEq = sigma : dEpsilonP

    and the void fraction grows by

        df = (1 - f) tr(dEpsilonP) + A dEpsilonPEq + kw f omega dEpsilonQ

    the terms being growth, strain-controlled nucleation (Chu and
    Needleman, with A = fN/(sN sqrt(2 pi)) exp(-((epsilonPEq - epsilonN)/sN)^2
    /2) under a tensile mean stress) and the shear term of Nahshon and
    Hutchinson, with omega = 1 - (27 J3/(2 q^3))^2 and dEpsilonQ the
    equivalent deviatoric plastic strain increment. The shear term is
    written with f, as Nahshon and Hutchinson write it, where Whelan's
    implementation used f*. With coalescence, the
    effective void fraction is

        f* = f                                for f <= fC
        f* = fC + (fu - fC)(f - fC)/(fF - fC) for f >  fC

    and f* = f without it, fu being the f* at which the yield surface
    shrinks to a point, 1/q1 for the usual q3 = q1^2. When the void
    fraction would take f* past 0.99 fu, the point fails: f is held there,
    where, for q3 = q1^2, the material keeps about 1% of its strength.

    The return mapping and the void fraction update are solved together at
    each point, implicitly in time: for a trial void fraction the plastic
    strain increments are found by Newton's method, and the void fraction
    is then the root of f - f(f0, the increments) = 0, found by regula
    falsi. This is the converged state of the scheme of Whelan et al.,
    which lags the void fraction by one outer iteration, but it is a
    function of the kinematics and the old-time state alone, as a law here
    must be.

    Non-local form. With nonLocal classic, the yield function, the work
    equivalence and the stress are written instead with the implicit
    non-local void fraction fbar, the solution of

        fbar - div(l^2 grad(fbar)) = f

    with l the characteristic length and zero flux on the boundary, which
    the solid model solves once per outer iteration (it is the
    implicit-gradient field voidFractionNonLocal) and this law reads as an
    input. The local void fraction is then updated explicitly, and is held
    at the failed void fraction once it passes it.

    The law implements the small-strain update and the elastic predictor
    interface, so it is used at finite strain as the sub-law of logStrain,
    which is the formulation of Whelan et al.

    References:
      - A.L. Gurson, J. Eng. Mater. Technol. 99 (1977) 2-15.
      - V. Tvergaard & A. Needleman, Acta Metall. 32 (1984) 157-169.
      - C.C. Chu & A. Needleman, J. Eng. Mater. Technol. 102 (1980)
        249-256.
      - K. Nahshon & J.W. Hutchinson, Eur. J. Mech. A Solids 27 (2008) 1-17.
      - N. Aravas, Int. J. Numer. Meth. Eng. 24 (1987) 1395-1416.
      - A. Whelan, T. Tang, V. Pakrashi & P. Cardiff, Int. J. Numer. Meth.
        Eng. (2025), doi:10.1002/nme.7640.

SourceFiles
    linearElasticGTNPlasticMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.
    Based on the GTN laws of Andrew Whelan, UCD.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticGTNPlasticMechanicalConstitutiveLaw_H
#define linearElasticGTNPlasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "interpolationTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class linearElasticGTNPlasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearElasticGTNPlasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
public:

    // Public Classes

        //- The result of the return mapping at one point
        struct pointResult
        {
            //- Elastic strain at the end of the step
            symmTensor epsilonE;

            //- Plastic strain increment
            symmTensor DEpsilonP;

            //- Stress
            symmTensor sigma;

            //- Volumetric plastic strain increment, tr(DEpsilonP)
            scalar DEpsilonV;

            //- Equivalent deviatoric plastic strain increment
            scalar DEpsilonQ;

            //- Matrix equivalent plastic strain increment
            scalar DEpsilonPEq;

            //- Mean stress at the end of the step
            scalar p;

            //- Equivalent stress at the end of the step
            scalar q;

            //- Void fraction at the end of the step: the local one
            scalar f;

            //- Void fraction the yield function is written with: f in the
            //  local form, the non-local void fraction in the non-local form
            scalar fYield;

            //- Yield stress of the matrix at the end of the step
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

        //- Table of matrix plastic strain versus yield stress
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Tvergaard parameters
        const scalar q1_;
        const scalar q2_;
        const scalar q3_;

        //- Initial void fraction
        const scalar f0_;

        //- Volume fraction of void-nucleating particles
        const scalar fN_;

        //- Standard deviation of the nucleation strain
        const scalar sN_;

        //- Mean nucleation strain
        const scalar epsilonN_;

        //- Shear damage parameter of Nahshon and Hutchinson
        const scalar kw_;

        //- Whether voids coalesce
        const bool coalescence_;

        //- Critical void fraction, at which coalescence starts
        scalar fC_;

        //- Void fraction at failure, at which f* reaches fu
        scalar fF_;

        //- Effective void fraction at which the yield surface vanishes
        scalar fu_;

        //- Void fraction of a failed point
        scalar fFailed_;

        //- Whether the void fraction is non-local
        bool nonLocal_;

        //- Characteristic length of the non-local form
        scalar characteristicLength_;

        //- Newton loop tolerance
        const scalar loopTol_;

        //- Maximum Newton iterations
        const label maxNewtonIter_;

        //- Finite-difference step for the hardening slope
        const scalar finiteDiff_;

        //- Absolute tolerance on the void fraction of the coupled solution
        const scalar voidFractionTol_;

        //- Maximum iterations of the coupled solution
        const label maxVoidFractionIter_;


    // Private Member Functions

        //- Yield stress of the matrix
        inline scalar yieldStress(const scalar epsilonPEq) const;

        //- The effective void fraction f* of a void fraction
        scalar effectiveVoidFraction(const scalar f) const;

        //- Return mapping from a trial elastic strain, at a given void
        //  fraction
        void returnMap
        (
            const symmTensor& epsilonETrial,
            const scalar epsilonPEq0,
            const scalar f,
            const scalar maxMagDEpsilon,
            pointResult& result
        ) const;

        //- The void fraction at the end of a returned step: the solution of
        //  the void growth equation for f, given the plastic increments
        scalar voidFractionUpdate
        (
            const pointResult& result,
            const scalar f0,
            const scalar epsilonPEq0
        ) const;

        //- The coupled return mapping and void fraction update at one point
        void localUpdate
        (
            const symmTensor& epsilonETrial,
            const scalar epsilonPEq0,
            const scalar f0,
            const scalar maxMagDEpsilon,
            pointResult& result
        ) const;

        //- The return mapping at a given non-local void fraction, and the
        //  local void fraction update
        void nonLocalUpdate
        (
            const symmTensor& epsilonETrial,
            const scalar epsilonPEq0,
            const scalar f0,
            const scalar fNonLocal,
            const scalar maxMagDEpsilon,
            pointResult& result
        ) const;

        //- No copy construct
        linearElasticGTNPlasticMechanicalConstitutiveLaw
        (
            const linearElasticGTNPlasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticGTNPlasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElasticGTNPlastic");


    // Static Data Members

        //- The name of the non-local void fraction field
        static const word nonLocalFieldName;


    // Constructors

        //- Construct from dictionary
        explicit linearElasticGTNPlasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearElasticGTNPlasticMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- Return the density
        virtual dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Return the bulk modulus of the matrix
        virtual dimensionedScalar kappa() const override
        {
            return kappa_;
        }

        //- Return the void fraction of a failed point
        scalar fFailed() const
        {
            return fFailed_;
        }

    // Edit

        //- Declare the plastic and void history
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

        //- The non-local void fraction, in the non-local form
        virtual wordList requiredScalarInputs() const override;

        //- The non-local void fraction, in the non-local form, whose source
        //  is the local void fraction
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

        //- Count the yielding and failed points of a step
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
