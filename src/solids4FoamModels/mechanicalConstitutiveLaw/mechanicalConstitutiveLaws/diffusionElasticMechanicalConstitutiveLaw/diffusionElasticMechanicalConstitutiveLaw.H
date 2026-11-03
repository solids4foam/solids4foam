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
    Foam::diffusionElasticMechanicalConstitutiveLaw

Description
    Small-strain linear elastic law whose stiffness is scaled by the
    stiffnessScaleFactor field, for use with a solid model as a mesh motion
    solver (solids4foamSolidModel):

        sigma = d*(2*mu*epsilon + lambda*tr(epsilon)*I)

    where d is the stiffnessScaleFactor, a coupling input gathered from the
    field of that name, which the solids4foamSolidModel mesh motion solver
    calculates from a mesh motion diffusivity, e.g. quadratic inverseDistance,
    so that the mesh is stiffer, and deforms less, near the moving
    boundaries.

    The material parameters are those of linearElastic; kappa is accepted in
    place of K. The scalar tangent and the fourth-order tangent are those of
    linearElastic scaled by d, and the finite difference fourth-order tangent
    is that of the scaled law.

    The law does not separate its volumetric response, as kappa() returns the
    unscaled bulk modulus, so it cannot be used with a mixed
    displacement-pressure formulation, and it takes no initial stress, sigma0.

SourceFiles
    diffusionElasticMechanicalConstitutiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionElasticMechanicalConstitutiveLaw_H
#define diffusionElasticMechanicalConstitutiveLaw_H

#include "linearElasticMechanicalConstitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class diffusionElasticMechanicalConstitutiveLaw Declaration
\*---------------------------------------------------------------------------*/

class diffusionElasticMechanicalConstitutiveLaw
:
    public linearElasticMechanicalConstitutiveLaw
{
    // Private Member Functions

        //- Return the law dictionary with kappa, if given, as K, checking
        //  for entries that the law does not take
        static dictionary linearElasticDict(const dictionary& dict);

        //- No copy construct
        diffusionElasticMechanicalConstitutiveLaw
        (
            const diffusionElasticMechanicalConstitutiveLaw&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const diffusionElasticMechanicalConstitutiveLaw&
        ) = delete;


public:

    //- Runtime type information
    TypeName("diffusionElastic");


    // Constructors

        //- Construct from dictionary
        explicit diffusionElasticMechanicalConstitutiveLaw
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~diffusionElasticMechanicalConstitutiveLaw() = default;


    // Member Functions

    // Access

        //- The volumetric response is not separated
        virtual bool providesVolumetricSplit() const override
        {
            return false;
        }

    // Edit

        //- The stiffness scale factor is a coupling input
        virtual wordList requiredScalarInputs() const override;

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
