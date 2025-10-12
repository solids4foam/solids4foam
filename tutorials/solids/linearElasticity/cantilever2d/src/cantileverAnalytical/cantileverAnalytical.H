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

Function
    cantileverAnalytical

Description
    Class for calculating analytical displacements and traction for a bending
    cantilever test case, as described in:

    C.E. Augarde, A.J. Deeks, The use of Timoshenko’s exact solution for a
    cantilever beam in adaptive analysis. Finite Elements in Analysis and
    Design, 44, 2008, 595–601, 10.1016/j.finel.2008.01.010.

    It is assumed that the origin is at the centre of the fixed patch.

    The beam is L long and D deep (thick).

SourceFile
    cantileverAnalytical.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cantileverAnalytical_H
#define cantileverAnalytical_H

#include "vectorField.H"
#include "dictionary.H"
#include "symmTensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cantileverAnalytical Declaration
\*---------------------------------------------------------------------------*/

class cantileverAnalytical
{
    // Private data

        //- Load applied in the minus y direction at the other end of the beam
        const scalar P_;

        //- Young's modulus
        const scalar E_;

        //- Poisson's ratio
        const scalar nu_;

        //- Length of the beam
        const scalar L_;

        //- Depth of the beam
        const scalar D_;

        //- Second moment of area
        const scalar I_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const cantileverAnalytical&);

public:

    //- Runtime type information
    TypeName("cantileverAnalytical");


    // Constructors

        //- Empty constructor
        cantileverAnalytical();

        //- Copy constructor
        cantileverAnalytical(const cantileverAnalytical& ca);

        //- Construct from a dictionary
        cantileverAnalytical(const dictionary& dict);


    // Destructor

        virtual ~cantileverAnalytical()
        {}


    // Member Functions

        // Access

            //- Return the displacements at the given location
            vector displacement(const vector& location) const;

            //- Return the displacements at the given locations
            tmp<vectorField> displacement(const vectorField& locations) const;

            //- Return the tractions at the given location
            vector traction(const vector& location) const;

            //- Return the tractions at the given locations
            tmp<vectorField> traction(const vectorField& locations) const;

            //- Return the stress tensor at the given location
            symmTensor stress(const vector& location) const;

    // Write

        //- Write parameter dict to the given stream
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
