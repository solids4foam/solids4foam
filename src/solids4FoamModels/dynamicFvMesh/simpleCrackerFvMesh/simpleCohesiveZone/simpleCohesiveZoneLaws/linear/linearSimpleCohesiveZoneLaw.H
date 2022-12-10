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
    linearSimpleCohesiveZoneLaw

Description
    Linear simple cohesive law.

SourceFiles
    linearSimpleCohesiveZoneLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearSimpleCohesiveZoneLaw_H
#define linearSimpleCohesiveZoneLaw_H

#include "simpleCohesiveZoneLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearSimpleCohesiveZoneLaw Declaration
\*---------------------------------------------------------------------------*/

class linearSimpleCohesiveZoneLaw
:
    public simpleCohesiveZoneLaw
{

private:

    // Private data

        //- Critical separation distance
        dimensionedScalar deltaC_;

public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from components
        linearSimpleCohesiveZoneLaw
        (
            const word& cohesiveLawName,
            const dictionary& dict
        );

        //- Construct as copy
        linearSimpleCohesiveZoneLaw(const linearSimpleCohesiveZoneLaw&);

        //- Construct and return a clone
        virtual autoPtr<simpleCohesiveZoneLaw> clone() const
        {
            return autoPtr<simpleCohesiveZoneLaw>
            (
                new linearSimpleCohesiveZoneLaw(*this)
            );
        }


    // Destructor

        virtual ~linearSimpleCohesiveZoneLaw();


    // Member Functions

        //- Return reference to critical separation distance
        const dimensionedScalar& deltaC() const
        {
            return deltaC_;
        }

        //- Return current holding traction
        scalar traction(scalar delta) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
