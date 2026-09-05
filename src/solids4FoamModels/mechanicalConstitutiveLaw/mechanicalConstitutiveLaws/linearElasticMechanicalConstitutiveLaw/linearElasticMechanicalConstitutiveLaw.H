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

    The law is stateless and does not require any constitutive history
    variables. All material parameters (e.g. Young's modulus and
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
    - No constitutive state variables are required.
    - Tangent-related quantities are not computed at this stage, but
      may be added later if required.


SourceFiles
    linearElasticMechanicalConstitutiveLaw.C

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

        //- Poisson's ration
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

    // Edit

        //- Declare the initial stress, which is added to every evaluation and
        //  is zero unless a sigma0 field is supplied
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

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
