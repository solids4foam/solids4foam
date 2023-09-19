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
    MooneyRivlinElastic

Description
    Three-parameter Mooney-Rivlin hyperelastic solid material, model proposed in
    Rivlin and Saunders (1951) with strain-energy function for the isochoric
    part given by:

        Psi = c10*(isoI1 - 3) + c01*(ispI2 - 3) + c11*(isoI1 - 3)(isoI2 - 3)

    The Cauchy stress tensor for the compressible case is then given by:

        J*sigma = 0.5*K*(J^2 - 1)*I
                  + 2*(
                          c10
                        + c11*(isoI2 - 3)
                     )*dev(isoB)
                  -  2*(
                          c01
                        + c11*(isoI1 - 3)
                     )*dev(isoB)^-1

    where:

    J         Jacobian of the deformation gradient
    F         deformation gradient
    isoB      isochoric left Cauchy-Green deformation tensor = J^(-2/3)*(F & F.T)
    c01, c10,
    c11       Material parameters of the solid
    isoI1     First principal invariant of the left Cauchy-Green deformation
    isoI2     Second principal invariant of the left Cauchy-Green deformation

    Includes support for pre-stress (S0) field and heterogeneous material
    properties.

    Reference
    Rivlin, R. S., & Saunders, D. W. (1951).
    Large elastic deformations of isotropic materials VII. Experiments on the
        deformation of rubber.
    Philosophical Transactions of the Royal Society of London A: Mathematical,
        Physical and Engineering Sciences

SourceFiles
    MooneyRivlinElastic.C

Author
    Iago Oliveira, based on code of Philip Cardiff, UCD.
    Ported by Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MooneyRivlinElastic_H
#define MooneyRivlinElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class MooneyRivlinElastic Declaration
\*---------------------------------------------------------------------------*/

class MooneyRivlinElastic
:
    public mechanicalLaw
{
    // Private data

        // First material parameter
        volScalarField c10_;

        // Second material parameter
        volScalarField c01_;

        // Third material parameter
        volScalarField c11_;

        // First material parameter for surface field
        surfaceScalarField c10f_;

        // Second material parameter for surface field
        surfaceScalarField c01f_;

        // Third material parameter for surface field
        surfaceScalarField c11f_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MooneyRivlinElastic(const MooneyRivlinElastic&);

        //- Disallow default bitwise assignment
        void operator=(const MooneyRivlinElastic&);

public:

    //- Runtime type information
    TypeName("MooneyRivlinElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        MooneyRivlinElastic
        (
            const word& name,
            const fvMesh& mesh,
            dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~MooneyRivlinElastic()
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
