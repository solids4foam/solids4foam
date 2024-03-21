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
    OgdenElastic

Description
    Ogden hyper-elasticity. This is a third order Ogden model.

    The Cauchy stress is given as:

    sigma =  s - (mu1_ + mu2_ + mu3_)*I - p*I

    and

    pbi = mu1*sqrt(lambdai)^alpha1 + mu2*sqrt(lambdai)^alpha2 + mu3*sqrt(lambdai)^alpha3

    where
    p             hydrostatic pressure (sigmaHyd == -p)
    I             indentity tensor
    K             bulk modulus
    mu1           shear modulus parameter 1
    mu2           shear modulus parameter 2
    mu3           shear modulus parameter 3
    s             rotated principal stress tensor
    prinStress    principal stress tensor
    pbi           principal stresses i = x, y, z
    lambdai       principal stretches i = x, y, z
    J             Jacobian of the deformation gradient

    Includes support for pre-stress (S0) field and heterogeneous material
    properties.

SourceFiles
    OgdenElastic.C

Author
    Joszef Nagy, Eulerian Solutions.
    Philip Cardiff, UCD.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef OgdenElastic_H
#define OgdenElastic_H

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class OgdenElastic Declaration
\*---------------------------------------------------------------------------*/

class OgdenElastic
:
    public mechanicalLaw
{
    // Private data

        // Shear modulus parameter 1
        const dimensionedScalar mu1_;

        // Shear modulus paramter 2
        const dimensionedScalar mu2_;

        // Shear modulus paramter 3
        const dimensionedScalar mu3_;

        // Exponent 1
        const dimensionedScalar alpha1_;

        // Exponent 2
        const dimensionedScalar alpha2_;

        // Exponent 3
        const dimensionedScalar alpha3_;

        // Bulk modulus
        const dimensionedScalar K_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OgdenElastic(const OgdenElastic&);

        //- Disallow default bitwise assignment
        void operator=(const OgdenElastic&);

public:

    //- Runtime type information
    TypeName("OgdenElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        OgdenElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const label lawI = -1
        );


    // Destructor

        virtual ~OgdenElastic()
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
