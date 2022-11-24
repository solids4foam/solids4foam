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
    YeohElastic

Description
    Yeoh hyperelastic model with three parameters, c1, c2 and c3, with the
    following strain energy function for the isochoric part given by:

        Psi = c1*(isoI1 - 3) + c2(isoI1 - 3)^2 + c3(isoI1 - 3)^3

    The Cauchy stress tensor for the compressible case is then given by:

        J*sigma = 0.5*K*(J^2 - 1)*I
                  + 2*(
                          c1
                        + 2*c2*(isoI1 - 3)
                        + 3*c3*(isoI1 - 3)
                     )*dev(isoB)

    where:

    J        Jacobian of the deformation gradient
    F        deformation gradient
    isoB     isochoric left Cauchy-Green deformation tensor = J^(-2/3)*(F & F.T)
    isoI1    first principal invariant of isoB

    Includes support for pre-stress (sigma0) field and heterogeneous material
    properties.

SourceFiles
    YeohElastic.C

Author
    Iago Oliveira, based on code by Philip Cardiff, UCD.
    Ported by Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef YeohElastic_H
#define YeohElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class YeohElastic Declaration
\*---------------------------------------------------------------------------*/

class YeohElastic
:
    public mechanicalLaw
{
    // Private data

        // First material parameter
        const volScalarField c1_;

        // First material parameter
        const volScalarField c2_;

        // Third material parameter
        const volScalarField c3_;

        // First material parameter for surface field
        const surfaceScalarField c1f_;

        // Second material parameter for surface field
        const surfaceScalarField c2f_;

        // Third material parameter for surface field
        const surfaceScalarField c3f_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YeohElastic(const YeohElastic&);

        //- Disallow default bitwise assignment
        void operator=(const YeohElastic&);

public:

    //- Runtime type information
    TypeName("YeohElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        YeohElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~YeohElastic()
        {}


    // Member Functions

        //- Return the implicit stiffness
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
