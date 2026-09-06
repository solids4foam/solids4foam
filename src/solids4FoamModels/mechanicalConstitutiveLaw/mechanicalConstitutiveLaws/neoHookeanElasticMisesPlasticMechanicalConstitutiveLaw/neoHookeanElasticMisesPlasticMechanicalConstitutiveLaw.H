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
    Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw

Description
    Finite-strain elasto-plastic constitutive law with compressible
    neo-Hookean elasticity and isotropic J2 (von Mises) plasticity.

    The formulation is the multiplicative one of Simo & Hughes (1998): the
    isochoric elastic left Cauchy-Green tensor bEbar is carried as history and
    updated by an exponential-map radial return in the intermediate
    configuration,

        relFbar    = relJ^(-1/3) relF,     relF = F F0^-1,   relJ = J/J0
        bEbarTrial = relFbar bEbar0 relFbar^T
        sTrial     = mu dev(bEbarTrial)
        fTrial     = |sTrial| - sqrt(2/3) J sigmaY

    with, on yielding, s = sTrial - 2 mu DEpsilonP and the Cauchy stress

        sigma = ( p I + s ) / J,    p = K (J^2 - 1)/2.

    Both perfect plasticity and isotropic hardening are supported through a
    yield-stress versus equivalent plastic strain curve. A two-point curve is
    treated as linear hardening with a closed-form return; three or more points
    are treated as nonlinear and solved with a local Newton iteration.

    Note on the hydrostatic stress. The legacy neoHookeanElasticMisesPlastic
    law passes its hydrostatic stress through updateSigmaHyd, which optionally
    solves a Laplacian smoothing equation over the whole field. As for
    MooneyRivlinElastic, that is deliberately not done here: a
    mechanicalConstitutiveLaw is a pure function of integration-point
    kinematics and old-time state, and pressure smoothing stabilises the
    discretisation rather than describing the material, so it belongs to the
    solid model. Cases setting smoothPressure will therefore not reproduce
    their legacy results through this law until the solid model provides the
    smoothing.

    The legacy law also under-relaxes the plastic strain increment between
    outer iterations, which is solver-loop behaviour with no equivalent here,
    and is likewise absent.

    This class implements only the constitutive response at a collection
    of integration points. Storage and management of history variables
    (bEbar, equivalent plastic strain, yield stress) are handled via
    mechanicalConstitutiveLawState objects owned by the
    mechanicalConstitutiveLawManager.

    Elastic parameters may be specified as:
      - Young's modulus (E) and Poisson's ratio (nu), or
      - Shear modulus (mu) and bulk modulus (K).

    References:
      - J.C. Simo & T.J.R. Hughes, Computational Inelasticity, Springer, 1998.
      - M.B. Rubin & A. Attia, Int. J. Numer. Meth. Eng., 39:309-320, 1996.
      - M. Hollenstein, M. Jabareen & M.B. Rubin, Comput. Mech., 52:649-667,
        2013.

SourceFiles
    neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw_H
#define neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"
#include "interpolationTable.H"
#include "Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw
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

        //- Enforce det(bEbar) == 1 when updating bEbar, by solving a cubic
        //  for its spherical part rather than taking the trial value
        const Switch updateBEbarConsistent_;

    // Private Member Functions

        //- Kirchhoff yield stress as a function of equivalent plastic strain.
        //  The curve is given as Cauchy stress against true strain, so it is
        //  scaled by J to give the Kirchhoff yield stress the return mapping
        //  is expressed in
        inline scalar yieldStress
        (
            const scalar epsilonPEq,
            const scalar J
        ) const;

        //- Yield function value
        inline scalar yieldFunction
        (
            const scalar epsilonPEqOld,
            const scalar magSTrial,
            const scalar dLambda,
            const scalar muBar,
            const scalar J
        ) const;

        //- Newton iteration for plastic multiplier
        inline void newtonLoop
        (
            scalar& dLambda,
            scalar& curSigmaY,
            const scalar epsilonPEqOld,
            const scalar magSTrial,
            const scalar muBar,
            const scalar J,
            const scalar maxMagDEpsilon
        ) const;

        //- Spherical part of bEbar that enforces det(bEbar) == 1, found by
        //  solving a cubic. Taking the trial value instead, as in Simo &
        //  Hughes, does not preserve the determinant
        inline scalar IbarConsistent(const symmTensor& devBEbar) const;

        //- No copy construct
        neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw
        (
            const neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("neoHookeanElasticMisesPlastic");

    // Constructors

        //- Construct from dictionary
        explicit neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw() = default;


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

        //- Initialise constitutive state variables
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

        //- Finite-strain constitutive evaluation.
        //  There is deliberately no small-strain overload: this law is
        //  formulated multiplicatively and has no small-strain form
        virtual void evaluate
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
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
