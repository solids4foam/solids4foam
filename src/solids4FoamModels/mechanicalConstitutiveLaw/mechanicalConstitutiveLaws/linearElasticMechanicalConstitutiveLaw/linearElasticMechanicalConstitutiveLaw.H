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
    Foam::linearElasticMechanicalConstitutiveLaw

Description
    Small-strain linear elastic mechanical constitutive law.

    This class implements an isotropic linear elastic constitutive
    relation in the small-strain regime, computing the Cauchy stress
    as a linear function of the symmetric displacement gradient.

    The law keeps no constitutive history. It reads an optional residual
    stress, sigma0, from the dictionary or, failing that, from a field of the
    same name in the case. All material parameters (e.g. Young's modulus and
    Poisson's ratio) are read from the supplied dictionary at
    construction time.

    Stress evaluation is performed using the
    smallStrainMechanicalConstitutiveLawKinematics interface. Finite-
    strain formulations are not supported by this law.

    Usage:
    The linear elastic law is selected via the mechanicalConstitutiveLaw
    run-time selection mechanism and used by the
    mechanicalConstitutiveLawManager during stress updates.

    Notes:
    - This law assumes infinitesimal strains and rotations.
    - It provides the scalar tangent 2*mu + lambda, and (4/3)*mu for the
      deviatoric surrogate of the mixed formulation, and the fourth-order
      tangent both analytically and by finite differences.


SourceFiles
    linearElasticMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMechanicalConstitutiveLaw_H
#define linearElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class linearElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMechanicalConstitutiveLaw
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

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- First Lamé parameter
        dimensionedScalar lambda_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;

        //- Uniform initial stress, taken from the law's dictionary. A sigma0
        //  field, if the case supplies one, overrides this point by point,
        //  unless the dictionary gave a value, which is final
        symmTensor sigma0_;

        //- Whether the dictionary gave sigma0 explicitly
        bool sigma0FromDict_;


    // Private Member Functions

        //- No copy construct
        linearElasticMechanicalConstitutiveLaw
        (
            const linearElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElastic");

    // Constructors

        //- Construct from dictionary
        explicit linearElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearElasticMechanicalConstitutiveLaw() = default;


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

        //- Return Young's modulus.
        //  For a solid model whose own formulation is written in terms of E
        //  and nu, such as a Kirchhoff plate's bending stiffness. Both are
        //  set whichever pair the case gave, (E, nu) or (mu, K)
        dimensionedScalar E() const
        {
            return E_;
        }

        //- Return Poisson's ratio. See E()
        dimensionedScalar nu() const
        {
            return nu_;
        }

        //- Return the shear modulus. For a solid model that assembles the
        //  Lame constants into its own matrix rather than asking for a stress
        dimensionedScalar mu() const
        {
            return mu_;
        }

        //- Return the first Lame parameter. See mu()
        dimensionedScalar lambda() const
        {
            return lambda_;
        }

    // Edit

        //- Declare the initial stress, which is added to every evaluation and
        //  is zero unless a sigma0 field is supplied
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

        //- The stress is a deviatoric part plus K*tr(epsilon), so the two
        //  can be handed back separately
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
