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
    Foam::barycentric

Description
    Minimal four-component barycentric coordinates for foam-extend, which is
    the only supported OpenFOAM variant that does not provide this class.
    Other supported OpenFOAM variants use their native implementation.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef foamExtendBarycentric_H
#define foamExtendBarycentric_H

#ifdef FOAMEXTEND

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class barycentric Declaration
\*---------------------------------------------------------------------------*/

class barycentric
{
    // Private data

        //- First coordinate
        scalar a_;

        //- Second coordinate
        scalar b_;

        //- Third coordinate
        scalar c_;

        //- Fourth coordinate
        scalar d_;


public:

    // Constructors

        //- Construct with all coordinates set to zero
        barycentric()
        :
            a_(0.0),
            b_(0.0),
            c_(0.0),
            d_(0.0)
        {}

        //- Construct from four coordinates
        barycentric
        (
            const scalar a,
            const scalar b,
            const scalar c,
            const scalar d
        )
        :
            a_(a),
            b_(b),
            c_(c),
            d_(d)
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

        //- Return the fourth coordinate
        const scalar& d() const
        {
            return d_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
#endif

// ************************************************************************* //
