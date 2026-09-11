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
    Foam::OgdenElasticMechanicalConstitutiveLaw

Description
    Ogden hyperelasticity with three terms.

    The principal stresses come from the principal stretches,

        p_i = sum_k mu_k * lambda_i^alpha_k

    and are rotated back with the eigenvectors of the right Cauchy-Green
    tensor. The deviatoric part of that, less the linearised shear modulus,
    plus a volumetric term, gives the Cauchy stress.

SourceFiles
    OgdenElasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef OgdenElasticMechanicalConstitutiveLaw_H
#define OgdenElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class OgdenElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class OgdenElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Ogden shear moduli
        dimensionedScalar mu1_;
        dimensionedScalar mu2_;
        dimensionedScalar mu3_;

        //- Ogden exponents
        dimensionedScalar alpha1_;
        dimensionedScalar alpha2_;
        dimensionedScalar alpha3_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Linearised shear modulus, mu1 + mu2 + mu3
        dimensionedScalar mu_;


    // Private Member Functions

        //- No copy construct
        OgdenElasticMechanicalConstitutiveLaw
        (
            const OgdenElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=(const OgdenElasticMechanicalConstitutiveLaw&) = delete;


public:

    //- Runtime type information
    TypeName("OgdenElastic");

    // Constructors

        //- Construct from dictionary
        explicit OgdenElasticMechanicalConstitutiveLaw(const dictionary& dict);


    //- Destructor
    virtual ~OgdenElasticMechanicalConstitutiveLaw() = default;


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
            return K_;
        }

        //- The initial stress, as for linearElastic
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

    // Edit

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
