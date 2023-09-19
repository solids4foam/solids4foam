/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    rotationInvariantLinearElastic

Description
    Rotation-invarient linear elastic where the stress tensor is calculated as
    the so-called "Biot strain" tensor is calculated as

        biotStrain = U - I

    and the stress as

        sigma = 2*mu*biotStrain + lambda*tr(biotStrain)*I

    where
    U        right stretch tensor
    I        indentity tensor
    mu       shear modulus, aka second Lame parameter
    lambda   first Lame parameter

    Reference
    G Dhondt, H-P Hackenberg (2021) Use of a rotation-invariant linear strain
    measure for linear elastic crack propagation calculations, Engineering
    Fracture Mechanics, 10.1016/j.engfracmech.2021.107634.


SourceFiles
    rotationInvariantLinearElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef rotationInvariantLinearElastic_H
#define rotationInvariantLinearElastic_H

#include "mechanicalLaw.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "surfaceFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class rotationInvariantLinearElastic
:
    public mechanicalLaw
{
    // Private data

        //- Biot strain, defined as (U - I)
        //  This is a small strain measure strain measure, which is rotation
        //  invariant
        volSymmTensorField biotStrain_;

         //- Shear modulus (second Lame parameter)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- First Lame parameter
        dimensionedScalar lambda_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        rotationInvariantLinearElastic(const rotationInvariantLinearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const rotationInvariantLinearElastic&);

public:

    //- Runtime type information
    TypeName("rotationInvariantLinearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        rotationInvariantLinearElastic
        (
            const word& name,
            const fvMesh& mesh,
            dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~rotationInvariantLinearElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
