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
    Abstract base class providing a common interface for numerical
    quadrature rules on reference geometric entities.

    The quadrature class defines storage and access for quadrature
    points and associated weights, expressed in the local coordinate
    system of the reference element. It supports integration of
    scalar- and vector-valued quantities by evaluating them at the
    quadrature points and forming the weighted sum.

    Derived classes implement specific quadrature rules for particular
    reference geometries (e.g. line, triangle, tetrahedron), including
    the placement of quadrature points and the corresponding weights.
    The quadrature order is specified at construction time and controls
    the accuracy of the numerical integration.

    This class is intended for use in finite-volume and finite-element
    style integrations within solids4foam, including high-order
    reconstruction, projection, and weak-form evaluations.

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

        //- Get reference to quadrature points weights
        virtual const List<scalar>& weights() const = 0;

        //- Number of quadrature points
        virtual label nPoints() const = 0;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
