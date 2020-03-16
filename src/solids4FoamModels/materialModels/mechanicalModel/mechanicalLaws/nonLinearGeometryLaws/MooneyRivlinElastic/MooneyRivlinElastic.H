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
    MooneyRivlinElastic

Description
    Mooney Rivlin hyper-elasticity, for example, as described at: https://www.\
    brown.edu/Departments/Engineering/Courses/En221/Notes/Elasticity/\
    Elasticity.htm. The Cauchy stress is given as:

        sigma =
            K*(J - 1)*I
          + (mu1/J^(5/3))*(b - tr(b)*I/3)
          + (mu2/J^(7/3))*(tr(b)*b - (tr(b)^2)*I/3 - (b & b) + (b && b)*I/3)

    where
    K        material property
    mu1      material property
    mu2      material property
    J        Jacobian of the deformation gradient
    I        indentity tensor
    b        left Cauchy–Green deformation tensor = F & F.T()
    F        deformation gradient

    In the limit of small strains, the shear modulus is mu = mu1 + mu2 and K is
    the bulk modulus.

SourceFiles
    MooneyRivlinElastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef MooneyRivlinElastic_H
#define MooneyRivlinElastic_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class MooneyRivlinElastic
:
    public mechanicalLaw
{
    // Private data

        // Initial density
        const dimensionedScalar rho_;

        // Shear modulus parameter 1
        const dimensionedScalar mu1_;

        // Shear modulus paramter 2
        const dimensionedScalar mu2_;

        // Bulk modulus parameter
        const dimensionedScalar K_;

        // Total deformation gradient volField
        volTensorField* FPtr_;

        // Total deformation gradient surfaceField
        surfaceTensorField* FfPtr_;

    // Private Member Functions

        //- Make the F field
        void makeF();

        //- Return a reference to the F field
        volTensorField& F();

        //- Make the Ff field
        void makeFf();

        //- Return a reference to the Ff field
        surfaceTensorField& Ff();

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
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~MooneyRivlinElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

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
