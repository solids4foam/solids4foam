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
    Foam::neoHookeanElasticMechanicalConstitutiveLaw

Description
    Finite-strain isotropic neo-Hookean elastic constitutive law.

    This class implements a compressible neo-Hookean material model for
    finite-strain solid mechanics. The constitutive response is formulated
    in terms of the deformation gradient and associated finite-strain
    kinematic measures, and returns the true (Cauchy) stress.

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
    neoHookeanElasticMechanicalConstitutiveLaw.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElasticMechanicalConstitutiveLaw_H
#define neoHookeanElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"
#include "dimensionedScalar.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class neoHookeanElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private Data

        //- Density
        const dimensionedScalar rho_;

        //- Young's modulus
        const dimensionedScalar E_;

        //- Poisson's ration
        const dimensionedScalar nu_;

        //- First Lamé parameter
        dimensionedScalar lambda_;

        //- Second Lamé parameter (shear modulus)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar kappa_;


    // Private Member Functions

        //- No copy construct
        neoHookeanElasticMechanicalConstitutiveLaw
        (
            const neoHookeanElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const neoHookeanElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("linearElastic");

    // Constructors

        //- Construct from dictionary
        explicit neoHookeanElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~neoHookeanElasticMechanicalConstitutiveLaw() = default;


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
