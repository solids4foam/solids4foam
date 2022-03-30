/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    orthotropicLinearElastic

Description
    Orthtropic Hookean linear elastic mechanical law, where the stress is
    related to the strain by nine independent elastic parameters: three
    stiffnesses, three shear modulii, and three Poisson's ratios.

    Properties are specified for local coordinate system and are then
    rotated to the global cartesian coordinates based on the specified
    material direction fields.

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

        //- Fourth order elasticity tensor volField
        //  Relative to the global Cartesian axis
        mutable volSymmTensor4thOrderField* elasticCPtr_;

        //- Fourth order elasticity tensor surfaceField
        //  Relative to the global Cartesian axis
        mutable surfaceSymmTensor4thOrderField* elasticCfPtr_;

        //- Material local x direction field
        volVectorField matDirX_;

        //- Material local y direction field
        volVectorField matDirY_;

        //- Material local z direction field
        volVectorField matDirZ_;


    // Private Member Functions

        //- Make elasticC field
        void makeElasticC() const;

        //- Return elasticC field
        const volSymmTensor4thOrderField& elasticC() const;

        //- Make elasticCf field
        void makeElasticCf() const;

        //- Return elasticCf field
        const surfaceSymmTensor4thOrderField& elasticCf() const;

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
            const nonLinearGeometry::nonLinearType& nonLinGeom
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
