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
    Foam::linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw

Description
    Mohr-Coulomb plasticity with linear elasticity, in the return-mapping form
    of the legacy law: a trial stress from the strain increment, then a mapping
    back onto the yield surface, its apex, or one of its edges.

    The stress is carried as a variation from the initial stress,

        sigma = deltaSigma + sigma0

    with deltaSigma held as history. The trial stress is built on the previous
    step's deltaSigma rather than on the incoming stress, which is what makes
    this law work underneath poroMechanicalLaw: there the incoming stress is
    the effective stress, and the pore pressure is subtracted after this law
    has finished.

SourceFiles
    linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw_H
#define linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- First Lame parameter
        dimensionedScalar lambda_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Friction angle, in degrees
        dimensionedScalar varPhi_;

        //- Cohesion
        dimensionedScalar c_;

        //- Dilation angle, in degrees
        dimensionedScalar varPsi_;

        //- Derived plasticity parameters
        scalar k_;
        scalar m_;
        vector a_;
        vector b_;
        symmTensor C_;
        symmTensor invC_;
        vector rp_;
        vector r_lf1_;
        vector r_lf2_;
        vector r_lg1_;
        vector r_lg2_;
        vector sigma_a_;


    // Private Member Functions

        //- Map a trial stress back onto the yield surface, and record whether
        //  the point is yielding. Carried over from the legacy law unchanged:
        //  it is a function of one point's stress and the material constants,
        //  and touches no fields
        void calculateStress(symmTensor& sigma, scalar& activeYield) const;

        //- Principal stresses and their directions. Carried over from the
        //  legacy law unchanged
        void calculateEigens
        (
            vector& sigma_prin,
            tensor& ev,
            const symmTensor sigma
        ) const;

        //- No copy construct
        linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw
        (
            const linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElasticMohrCoulombPlastic");

    // Constructors

        //- Construct from dictionary
        explicit linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw()
        = default;


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

        //- The stress variation and the yield flag are history; the initial
        //  stress is supplied
        virtual void declareState
        (
            mechanicalConstitutiveLawStateSpec& spec
        ) const override;

    // Edit

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
