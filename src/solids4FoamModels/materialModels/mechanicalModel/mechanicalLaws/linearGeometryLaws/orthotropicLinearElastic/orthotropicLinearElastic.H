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
    orthotropicLinearElastic

Description
    Orthtropic Hookean linear elastic mechanical law, where the stress is
    related to the strain by nine independent elastic parameters: three
    stiffnesses, three shear modulii, and three Poisson's ratios.

    Properties are specified for local coordinate system and are then
    rotated to the global cartesian coordinates based on the specified
    material direction fields.

    In foam-extend, we use the symmTensor4thOrder class to store the 4th order
    stiffness tensor (9 components). In OpenFOAM, we use a tensor to store these
    9 components with the following ordering assumed:
        XXXX    ->    XX    ->    0
        XXYY    ->    XY    ->    1
        XXZZ    ->    XZ    ->    2
        YYYY    ->    YX    ->    3
        YYZZ    ->    YY    ->    4
        ZZZZ    ->    YZ    ->    5
        XYXY    ->    ZX    ->    6
        YZYZ    ->    ZY    ->    7
        ZXZX    ->    ZZ    ->    8

    For more details, see:
    P. Cardiff, A. Karać, A. Ivanković: A Large Strain Finite Volume Method
    for Orthotropic Bodies with General Material Orientations. Computer Methods
    in Applied Mechanics and Engineering, 01/2014, 268(1):318-335.
    DOI: 10.1016/j.cma.2013.09.008.

SourceFiles
    orthotropicLinearElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef orthotropicLinearElastic_H
#define orthotropicLinearElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class orthotropicLinearElastic Declaration
\*---------------------------------------------------------------------------*/

class orthotropicLinearElastic
:
    public mechanicalLaw
{
    // Private data

        //- Modulus of elasticity in the local direction 1
        dimensionedScalar E1_;

        //- Modulus of elasticity in the local direction 2
        dimensionedScalar E2_;

        //- Modulus of elasticity in the local direction 3
        dimensionedScalar E3_;

        //- Poisson's ratio
        //  Contraction in the 2 direction when extended in the 1 direction
        dimensionedScalar nu12_;

        //- Poisson's ratio
        //  Contraction in the 3 direction when extended in the 2 direction
        dimensionedScalar nu23_;

        //- Poisson's ratio
        //  Contraction in the 1 direction when extended in the 3 direction
        dimensionedScalar nu31_;

        //- Poisson's ratio
        //  Contraction in the 1 direction when extended in the 2 direction
        dimensionedScalar nu21_;

        //- Poisson's ratio
        //  Contraction in the 3 direction when extended in the 2 direction
        dimensionedScalar nu32_;

        //- Poisson's ratio
        //  Contraction in the 1 direction when extended in the 3 direction
        dimensionedScalar nu13_;

        //- Shear modulus in local 12 plane
        dimensionedScalar G12_;

        //- Shear modulus in local 23 plane
        dimensionedScalar G23_;

        //- Shear modulus in local 31 plane
        dimensionedScalar G31_;

#ifdef FOAMEXTEND
        //- Fourth order elasticity tensor volField
        //  Relative to the global Cartesian axis
        mutable autoPtr<volSymmTensor4thOrderField> elasticCPtr_;

        //- Fourth order elasticity tensor surfaceField
        //  Relative to the global Cartesian axis
        mutable autoPtr<surfaceSymmTensor4thOrderField> elasticCfPtr_;
#else
        //- Fourth order elasticity tensor volField
        //  Relative to the global Cartesian axis
        mutable autoPtr<volTensorField> elasticCPtr_;

        //- Fourth order elasticity tensor surfaceField
        //  Relative to the global Cartesian axis
        mutable autoPtr<surfaceTensorField> elasticCfPtr_;
#endif

        //- Material local x direction field
        volVectorField matDirX_;

        //- Material local y direction field
        volVectorField matDirY_;

        //- Material local z direction field
        volVectorField matDirZ_;


    // Private Member Functions

        //- Make elasticC field
        void makeElasticC() const;

        //- Make elasticCf field
        void makeElasticCf() const;

#ifdef FOAMEXTEND
        //- Return elasticC field
        const volSymmTensor4thOrderField& elasticC() const;

        //- Return elasticCf field
        const surfaceSymmTensor4thOrderField& elasticCf() const;
#else
        //- Return elasticC field
        const volTensorField& elasticC() const;

        //- Return elasticCf field
        const surfaceTensorField& elasticCf() const;
#endif

        //- Disallow default bitwise copy construct
        orthotropicLinearElastic(const orthotropicLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const orthotropicLinearElastic&);

public:

    //- Runtime type information
    TypeName("orthotropicLinearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        orthotropicLinearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const label lawI = -1,
            const solidSubMeshes* solidSubMeshes = NULL
        );


    // Destructor

        virtual ~orthotropicLinearElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
