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
    Foam::barycentric2D

Description
    Minimal three-component barycentric coordinates for foam-extend, which is
    the only OpenFOAM variant that does not provide this class.
    Other OpenFOAM variants use their native implementation.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef foamExtendBarycentric2D_H
#define foamExtendBarycentric2D_H

#ifdef FOAMEXTEND

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class barycentric2D Declaration
\*---------------------------------------------------------------------------*/

class barycentric2D
{
    // Private data

        //- First coordinate
        scalar a_;

        //- Second coordinate
        scalar b_;

        //- Third coordinate
        scalar c_;


public:

    // Constructors

        //- Construct with all coordinates set to zero
        barycentric2D()
        :
            a_(0.0),
            b_(0.0),
            c_(0.0)
        {}

        //- Construct from three coordinates
        barycentric2D
        (
            const scalar a,
            const scalar b,
            const scalar c
        )
        :
            a_(a),
            b_(b),
            c_(c)
        {}


    // Member Functions

        //- Return the first coordinate
        const scalar& a() const
        {
            return a_;
        }

        //- Return the second coordinate
        const scalar& b() const
        {
            return b_;
        }

        //- Return the third coordinate
        const scalar& c() const
        {
            return c_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
#endif

// ************************************************************************* //
