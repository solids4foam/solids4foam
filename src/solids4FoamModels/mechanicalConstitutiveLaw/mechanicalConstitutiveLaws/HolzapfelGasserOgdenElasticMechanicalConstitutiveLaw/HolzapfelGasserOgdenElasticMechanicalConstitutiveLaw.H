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
    Foam::HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw

Description
    Holzapfel-Gasser-Ogden anisotropic hyperelastic law: a neo-Hookean ground
    substance with two families of fibres, wound symmetrically about the
    circumferential direction.

        W = 0.5*mu*(I1bar - 3)
          + (k1/(2*k2))*sum_i (exp(k2*(I4i - 1)^2) - 1)
          + U(J)

    where the invariants are taken on the isochoric deformation
    Fbar = J^(-1/3)*F, and the two fibre families are

        N4 =  cos(fibreAngle)*Ec + sin(fibreAngle)*Ea
        N6 = -cos(fibreAngle)*Ec + sin(fibreAngle)*Ea

    The volumetric response is the same penalty the other hyperelastic laws
    here use,

        U(J) = 0.25*K*(J^2 - 1 - 2*ln(J)),   dU/dJ = 0.5*K*(J^2 - 1)/J

    which is what lets this law answer a mixed displacement-pressure
    formulation. The legacy HolzapfelGasserOgdenElastic has no volumetric term
    at all: it is written for exact incompressibility and reads the solid
    model's pressure out of the object registry to supply the whole spherical
    part of its stress. A law that reaches for a field the solid model happens
    to have registered cannot be evaluated for a tangent, restarted, or asked
    what it would do at a perturbed deformation, so that is not reproduced
    here. Instead the pressure enters the way it does for every other law in
    this framework: the solid model asks for the isochoric stress and the
    volumetric response separately and replaces the second with what it
    solved. A large bulk modulus approaches the incompressible limit the
    legacy law assumes outright.

    Nothing selects this law yet. coupledPressureDisplacementSolid is the only
    consumer of the legacy one, and it calls mechanical().correct() with no
    framework switch, so it cannot reach this. It does already carry a penalty
    term in its pressure equation - fvm::Sp(rKappa_, Dp_) with rKappa_ = 1/K -
    which is a reason to expect the shapes to match when that solid model is
    ported rather than evidence that they do. DESIGN-state-io.md section 25.2
    records the linearisation question that comes with it.

    The local basis is prescribed state. A case supplies Ec and Ea fields -
    calcLocCoordinates writes them - and a case whose fibres are the same
    everywhere may give single vectors in this dictionary instead. They are
    normalised here, so they need not be given normalised, and a zero-length
    basis vector is refused rather than divided by.

    The fibre term is not gated on I4 > 1. The legacy law is not gated either,
    and gating it is a change of material rather than of implementation, so it
    is left as it is rather than quietly improved.

    Dictionary entries:

    \verbatim
        type            HolzapfelGasserOgdenElastic;
        rho             rho [1 -3 0 0 0 0 0] <value>;
        mu              mu  [1 -1 -2 0 0 0 0] <value>;
        k1              k1  [1 -1 -2 0 0 0 0] <value>;
        k2              <value>;                    // dimensionless
        fibreAngle      <value>;                    // degrees
        bulkModulus     bulkModulus [1 -1 -2 0 0 0 0] <value>;

        // Optional: fibres the same everywhere, rather than Ec and Ea fields
        uniformLocalBasis yes;
        Ec              (1 0 0);
        Ea              (0 1 0);
    \endverbatim

SourceFiles
    HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw_H
#define HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- Density
        const dimensionedScalar rho_;

        //- Shear modulus of the ground substance
        const dimensionedScalar mu_;

        //- Fibre stiffness
        const dimensionedScalar k1_;

        //- Fibre stiffening exponent, dimensionless
        const scalar k2_;

        //- Angle of the fibre families from the circumferential direction,
        //  in radians. Read in degrees, as the legacy law reads it
        const scalar fibreAngle_;

        //- Bulk modulus, a penalty enforcing near incompressibility
        const dimensionedScalar bulkModulus_;

        //- The local basis to use where a case supplies no Ec and Ea fields.
        //  Zero when the case is expected to supply them, so that their
        //  absence is a zero-length direction and refused rather than
        //  quietly turning the material isotropic
        const vector EcDefault_;
        const vector EaDefault_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw
        (
            const HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=
        (
            const HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("HolzapfelGasserOgdenElastic");

    // Constructors

        //- Construct from dictionary
        explicit HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw() = default;


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
            return bulkModulus_;
        }

        //- This law separates its isochoric and volumetric responses.
        //  Every invariant is taken on the isochoric deformation, so the
        //  stress it produces depends on shape alone and the volumetric
        //  response is the penalty term and nothing else
        virtual bool providesVolumetricSplit() const override
        {
            return true;
        }

        //- Declare the local basis this law is given
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

#endif

// ************************************************************************* //
