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
    Foam::immersedBodyMotions::staticMotion

Description
    The immersed body does not move.

    Usage:
    \verbatim
    motion
    {
        type    static;
    }
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    staticMotion.C

\*---------------------------------------------------------------------------*/

#ifndef staticMotion_H
#define staticMotion_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                        Class staticMotion Declaration
\*---------------------------------------------------------------------------*/

class staticMotion
:
    public immersedBodyMotion
{
public:

    //- Runtime type information
    TypeName("static");


    // Constructors

        //- Construct from dictionary
        explicit staticMotion(const dictionary& dict);


    //- Destructor
    virtual ~staticMotion() = default;


    // Member Functions

        //- The body does not move
        virtual bool moving() const
        {
            return false;
        }

        //- Positions at time t of the surface points
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const;

        //- Velocity of the body at time t at the points x
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBodyMotions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
