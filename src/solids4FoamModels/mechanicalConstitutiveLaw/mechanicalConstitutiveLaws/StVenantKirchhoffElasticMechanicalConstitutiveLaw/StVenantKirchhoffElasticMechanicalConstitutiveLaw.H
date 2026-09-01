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
    Foam::StVenantKirchhoffElasticMechanicalConstitutiveLaw

Description
    Finite-strain isotropic St Venant-Kirchhoff elastic constitutive law.

    The second Piola-Kirchhoff stress is linear in the Green-Lagrange strain,

        E = (F^T F - I)/2,   S = 2 mu E + lambda tr(E) I,

    and the Cauchy stress follows by the push-forward sigma = F S F^T / J.

    This is the simplest finite-strain elastic law: it is exact for rigid body
    rotation, which is what distinguishes it from small-strain elasticity, but
    it is not polyconvex and should not be used at large compressive strain.

    The law is history-independent and does not store any internal state,
    but it is fully compatible with the mechanicalConstitutiveLaw framework
    and may be used alongside history-dependent laws in multi-material
    simulations.

    Material parameters are provided in terms of Young's modulus and
    Poisson's ratio, specified as dimensioned quantities and validated at
    construction time.

    The implementation is intended as a simple, robust finite-strain
    reference law and as a demonstrator for the finite-strain constitutive
    infrastructure within solids4foam.

SourceFiles
    StVenantKirchhoffElasticMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef StVenantKirchhoffElasticMechanicalConstitutiveLaw_H
#define StVenantKirchhoffElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class StVenantKirchhoffElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class StVenantKirchhoffElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private Data

        //- Density
        const dimensionedScalar rho_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ration
        dimensionedScalar nu_;

        //- First Lamé parameter
        dimensionedScalar lambda_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;


    // Private Member Functions

        //- No copy construct
        StVenantKirchhoffElasticMechanicalConstitutiveLaw
        (
            const StVenantKirchhoffElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const StVenantKirchhoffElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("StVenantKirchhoffElastic");

    // Constructors

        //- Construct from dictionary
        explicit StVenantKirchhoffElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~StVenantKirchhoffElasticMechanicalConstitutiveLaw() = default;


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

        //- Small-strain constitutive evaluation
        virtual void evaluate
        (
            const finiteStrainMechanicalConstitutiveLawKinematics& kin,
            mechanicalConstitutiveLawState& state,
            mechanicalConstitutiveLawResponse& response
        ) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
