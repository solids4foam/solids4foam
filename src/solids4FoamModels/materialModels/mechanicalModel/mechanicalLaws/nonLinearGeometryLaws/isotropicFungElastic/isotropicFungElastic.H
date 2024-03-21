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
    isotropicFungElastic

Description
    Isotropic Fung-like model with two parameters, c1 and c2, with
    strain energy function for the isochoric part given by:

        Psi = (c1/c2)*(exp(0.5*c2*(isoI1 - 3) - 1)

    The Cauchy stress tensor for the compressible case is then given by:

        J*sigma =   0.5*K*(J^2 - 1)*I
                  + c1*exp(c2*(isoI1 - 3))*dev(isoB)

    where:

    J        Jacobian of the deformation gradient
    F        deformation gradient
    isoB     isochoric left Cauchy-Green deformation tensor = J^(-2/3)*(F & F.T)
    isoI1    first principal invariant of isoB

    Includes support for pre-stress (S0) field and heterogeneous c1.

    Reference


SourceFiles
    isotropicFungElastic.C

Author
    Iago Oliveira, based on code by Philip Cardiff, UCD.
    Minor modifications by Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef isotropicFungElastic_H
#define isotropicFungElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class isotropicFungElastic Declaration
\*---------------------------------------------------------------------------*/

class isotropicFungElastic
:
    public mechanicalLaw
{
    // Private data

        // Shear modulus parameter
        volScalarField c1_;

        // Shear modulus parameter surface field
        surfaceScalarField c1f_;

        // First material parameter associated with I2
        const dimensionedScalar c2_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        isotropicFungElastic(const isotropicFungElastic&);

        //- Disallow default bitwise assignment
        void operator=(const isotropicFungElastic&);

public:

    //- Runtime type information
    TypeName("isotropicFungElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        isotropicFungElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const label lawI = -1
        );


    // Destructor

        virtual ~isotropicFungElastic()
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
