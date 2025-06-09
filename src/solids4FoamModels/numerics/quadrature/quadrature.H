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
    quadrature

Description
    Base quadrature class

SourceFile
    quadrature.C

Author
    Ivan Batistic, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef quadrature_H
#define quadrature_H

#include "List.H"
#include "vector.H"
#include "scalar.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class quadrature
{

protected:

    //- Order of integration
    const label order_;

public:

    // Constructors
    quadrature(const label& n);

    // Destructor
    virtual ~quadrature();

    // Member Functions

    // Access

        //- Order of integration
        virtual label order() const;

        //- Get reference to quadrature points in global Cartesian coordinates
        virtual const List<point>& points() const = 0;

        //- Get quadrature points in global Cartesian coordinates
        virtual const List<point> points() = 0;

        //- Get reference to quadrature points weights
        virtual const List<scalar>& weights() const = 0;

        //- Get quadrature points weights
        virtual const List<scalar> weights() = 0;

        //- Number of quadrature points
        virtual label nPoints() const = 0;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
