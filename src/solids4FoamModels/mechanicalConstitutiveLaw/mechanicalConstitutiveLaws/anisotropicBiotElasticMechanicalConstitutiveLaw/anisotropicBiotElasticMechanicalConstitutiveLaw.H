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
    Foam::anisotropicBiotElasticMechanicalConstitutiveLaw

Description
    Orthotropic linear elasticity, sigma = A : epsilon, with the stiffness
    given as Young's moduli, Poisson's ratios and shear moduli in the
    coordinate directions.

    Two forms, selected by the mesh rather than by the dictionary. The
    directions the mesh solves in are injected by the manager as solutionD,
    because a mechanical constitutive law is constructed from a dictionary and
    cannot ask the mesh itself.

SourceFiles
    anisotropicBiotElasticMechanicalConstitutiveLaw.C

Author
    Original author: Tian Tang, DTU.
    Ported to solids4foam by: Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef anisotropicBiotElasticMechanicalConstitutiveLaw_H
#define anisotropicBiotElasticMechanicalConstitutiveLaw_H

#include "mechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class anisotropicBiotElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class anisotropicBiotElasticMechanicalConstitutiveLaw
:
    public mechanicalConstitutiveLaw
{
    // Private data

        //- Density
        dimensionedScalar rho_;

        //- True when the z direction is not solved, i.e. a 2-D case
        bool model2d_;

        //- Orthotropic stiffness coefficients
        scalar A11_;
        scalar A22_;
        scalar A33_;
        scalar A44_;
        scalar A55_;
        scalar A66_;
        scalar A12_;
        scalar A21_;
        scalar A31_;
        scalar A23_;


    // Private Member Functions

        //- No copy construct
        anisotropicBiotElasticMechanicalConstitutiveLaw
        (
            const anisotropicBiotElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const anisotropicBiotElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("anisotropicBiotElastic");

    // Constructors

        //- Construct from dictionary
        explicit anisotropicBiotElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~anisotropicBiotElasticMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- Return the density
        virtual dimensionedScalar rho() const override
        {
            return rho_;
        }

        //- Return the bulk modulus.
        //  An orthotropic material has no single bulk modulus; the largest
        //  direct stiffness is what its impK uses and is the closest thing
        //  available
        virtual dimensionedScalar kappa() const override
        {
            return dimensionedScalar
            (
                "kappa", dimPressure, max(A11_, max(A22_, A33_))
            );
        }

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
