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
    Foam::isotropicFungElasticMechanicalConstitutiveLaw

Description
    Isotropic Fung Elastic model for biological soft tissues.

    The strain energy is exponential in a quadratic form of the Green-Lagrange
    strain:

         W = 0.5*c*(exp(Q) - 1)

         Q = a*tr(Ebar & Ebar)

    where Ebar is built from Fbar = J^(-1/3)*F. That is what lets a mixed
    displacement-pressure formulation replace the volumetric part with a solved
    pressure and still describe this material rather than a neighbouring one.

    This allows splitting of the energy into isochoric and volumetric components
        W = 0.5*c*(exp(Q) - 1) + U(J)

        U(J) = 0.25*K*(J^2 - 1 - 2*ln(J))

    SourceFiles
    isotropicFungElasticMechanicalConstitutiveLaw.C

Author
    Conor Moore, UCD, based on GuccioneElastic by Philip Cardiff

\*---------------------------------------------------------------------------*/

#ifndef isotropicFungElasticMechanicalConstitutiveLaw_H
#define isotropicFungElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class isotropicFungElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class isotropicFungElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- Density
        const dimensionedScalar rho_;

        //- Stress scale of the material multiplied by the exponential
        const dimensionedScalar c_;

        //- Stiffness coefficient of the exponential
        const scalar a_;

        //- Bulk modulus, a penalty enforcing near incompressibility
        const dimensionedScalar bulkModulus_;

        //- Small-strain shear modulus implied by c and a, c*a/2.
        //  Used for the implicit stiffness and nothing else
        const dimensionedScalar mu_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        isotropicFungElasticMechanicalConstitutiveLaw
        (
            const isotropicFungElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=
        (
            const isotropicFungElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("isotropicFungElastic");

    // Constructors

        //- Construct from dictionary
        explicit isotropicFungElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~isotropicFungElasticMechanicalConstitutiveLaw() = default;


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
        //  Q is written on the isochoric strain, so the deviatoric stress it
        //  produces depends on shape alone and the volumetric response is the
        //  penalty term and nothing else
        virtual bool providesVolumetricSplit() const override
        {
            return true;
        }

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
