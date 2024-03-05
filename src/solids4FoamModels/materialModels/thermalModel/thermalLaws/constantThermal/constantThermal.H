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
    constantThermal

Description
    Constant thermal properties

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    constantThermal.C

\*---------------------------------------------------------------------------*/

#ifndef constantThermal_H
#define constantThermal_H

#include "thermalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantThermal Declaration
\*---------------------------------------------------------------------------*/

class constantThermal
:
    public thermalLaw
{
    // Private data

        //- Specific heat capacity
        dimensionedScalar C_;

        //- Thermal conductivity
        dimensionedScalar k_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantThermal(const constantThermal&);

        //- Disallow default bitwise assignment
        void operator=(const constantThermal&);


public:

    //- Runtime type information
    TypeName("constant");

    // Static data members


    // Constructors

        //- Construct from dictionary
        constantThermal
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~constantThermal();


    // Member Functions

        //- Return specific heat capacity
        virtual tmp<volScalarField> C() const;

        //- Return thermal conductivity
        virtual tmp<volScalarField> k() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
