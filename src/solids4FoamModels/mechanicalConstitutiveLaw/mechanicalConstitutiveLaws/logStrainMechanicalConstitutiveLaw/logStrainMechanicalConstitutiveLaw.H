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
    Foam::logStrainMechanicalConstitutiveLaw

Description
    A small-strain law used at finite strain, through the logarithmic
    elastic strain.

    This is a composite: it owns a sub-law, written for small strain, and
    drives it with the elastic logarithmic (Hencky) strain, following the
    exponential-map return mapping of Eterovic and Bathe (1990) and Simo
    (1992). Over a step, with f = F & inv(F0) the relative deformation
    gradient and epsilonE0 the elastic strain at the start of the step,

        bETrial       = f & exp(2*epsilonE0) & f^T
        epsilonETrial = 0.5*log(bETrial)

    The sub-law takes epsilonETrial as its trial elastic strain, performs
    its return mapping on it exactly as it would at small strain, and gives
    back the elastic strain epsilonE the step ends at and the stress tau
    conjugate to it, which is the Kirchhoff stress. The Cauchy stress is

        sigma = tau/J

    The yield function of a plastic sub-law is therefore written on the
    Kirchhoff stress, and its hardening curve is read as Kirchhoff stress
    against equivalent logarithmic plastic strain.

    The multiplicative split F = Fe & Fp is exact here, not approximated,
    for an isotropic sub-law with scalar-valued history whose plastic flow
    keeps the elastic strain coaxial with its trial value, which is what
    providesElasticPredictor() declares. A law with kinematic hardening or
    a preferred direction does not qualify and is refused. With an elastic
    sub-law, linearElastic, this is Hencky's hyperelastic law.

    The elastic left Cauchy-Green tensor bE = exp(2*epsilonE) is this law's
    history, in its own state, alongside epsilonE itself; the sub-law's
    scalar history (equivalent plastic strain, yield stress, ...) lives in a
    child state, so the sub-law looks it up by its own unqualified names.
    Since the return mapping keeps epsilonE coaxial with epsilonETrial, bE is
    updated from the eigenvectors of bETrial, and each evaluation costs one
    eigen-decomposition per point.

    Restrictions:
      - finite strain only, so the law is used with a nonlinear geometry
        solid model;
      - no plane stress, whose out-of-plane stretch would have to be solved
        for at each point;
      - no prestress, since the elastic history starts at zero rather than
        at the prestressed state;
      - no mixed displacement-pressure split, and no analytical fourth-order
        tangent (fourthOrderFiniteDifference is available). The scalar
        tangent is the sub-law's small-strain one, which steers the
        segregated solvers but is not a consistent tangent.

    References:
      - G. Weber & L. Anand, Comput. Methods Appl. Mech. Eng. 79 (1990)
        173-202.
      - A.L. Eterovic & K.-J. Bathe, Int. J. Numer. Meth. Eng. 30 (1990)
        1099-1114.
      - J.C. Simo, Comput. Methods Appl. Mech. Eng. 99 (1992) 61-112.
      - A. Whelan, P. Cardiff et al., Int. J. Numer. Meth. Eng. (2025).

SourceFiles
    logStrainMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.
    Based on the logarithmic strain laws of Andrew Whelan, UCD.

\*---------------------------------------------------------------------------*/

#ifndef logStrainMechanicalConstitutiveLaw_H
#define logStrainMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class logStrainMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class logStrainMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- The small-strain law driven by the elastic logarithmic strain
        autoPtr<mechanicalConstitutiveLaw> subLawPtr_;

        //- Name of the child state the sub-law's history lives in
        const word subStateName_;


    // Private Member Functions

        //- The trial elastic logarithmic strain at each point, from the
        //  deformation over the step and the elastic left Cauchy-Green
        //  tensor it started at, and the eigenvectors of both, as the rows
        //  of a tensor
        void trialElasticStrain
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawState& state,
            Field<symmTensor>& epsilonETrial,
            Field<tensor>& eigenVectors
        ) const;

        //- No copy construct
        logStrainMechanicalConstitutiveLaw
        (
            const logStrainMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=(const logStrainMechanicalConstitutiveLaw&) = delete;


public:

    //- Runtime type information
    TypeName("logStrain");


    // Static Member Functions

        //- The logarithm of a symmetric positive-definite tensor
        static symmTensor logm(const symmTensor& A);

        //- The exponential of a symmetric tensor
        static symmTensor expm(const symmTensor& A);


    // Constructors

        //- Construct from dictionary
        explicit logStrainMechanicalConstitutiveLaw(const dictionary& dict);


    //- Destructor
    virtual ~logStrainMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- The density is the sub-law's
        virtual dimensionedScalar rho() const override
        {
            return subLawPtr_->rho();
        }

        //- As is the bulk modulus
        virtual dimensionedScalar kappa() const override
        {
            return subLawPtr_->kappa();
        }

        //- The coupling inputs are the sub-law's
        virtual wordList requiredScalarInputs() const override
        {
            return subLawPtr_->requiredScalarInputs();
        }

        //- As are the ones it could read but has not been asked to
        virtual wordList optionalScalarInputs() const override
        {
            return subLawPtr_->optionalScalarInputs();
        }

        //- A bound on the largest trial elastic strain over these points,
        //  which a plastic sub-law normalises its Newton residual by, as it
        //  would the largest strain at small strain
        virtual scalar localConvergenceScale
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            const mechanicalConstitutiveLawState& state
        ) const override;

        //- Whether the sub-law uses a convergence scale
        virtual bool usesConvergenceScale() const override
        {
            return subLawPtr_->usesConvergenceScale();
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

        //- Declare the elastic left Cauchy-Green tensor, this law's history,
        //  and the elastic logarithmic strain
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
