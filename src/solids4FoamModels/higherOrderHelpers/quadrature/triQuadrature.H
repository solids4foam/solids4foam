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
    triQuadrature

Description
    Numerical quadrature rules for integration over a two-dimensional
    reference triangular element.

    This class provides Gaussian-type quadrature rules for triangles,
    defining quadrature points in barycentric (area) coordinates along
    with the associated integration weights. The selected quadrature
    order determines the accuracy of polynomial integration over the
    reference triangle.

    The triQuadrature class is intended for surface and face-based
    integrations, including high-order reconstruction, projection, and
    weak-form evaluations on triangular elements.

    Source:
    D. A. Dunavant,
    High degree efficient symmetrical Gaussian quadrature rule for the triangle.
    International Journal for Numerical Methods in Engineering,
    21(6):1129–1148, jun 1985.

SourceFile
    triQuadrature.C

Author
    Ivan Batistic, UCD. All rights reserved.
    Pablo Castrillo, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef triQuadrature_H
#define triQuadrature_H

#include "triangle.H"
#include "quadrature.H"
#include "barycentric2D.H"
#include "Map.H"
#include "Field.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triQuadrature
:
    public quadrature
{

    // Private data

        //- Store triangle using triangle template
        triangle<point, point> tri_;

        //- Quadrature points in Cartesian coordinates
        List<point> points_;

        //- Quadrature points weights
        List<scalar> weights_;

    // Private Member Functions

        //- Map integration points from barycentric to cartesian coordinates
        tmp<Field<point>> barycentricToPoint
        (
             const List<barycentric2D>& localPts
        ) const;

    // Static data

        //- Maximum order of accuracy when integrating polynomials
        static constexpr label maxSupportedOrder = 10;

        //- Quadrature rules storage
        struct quadratureRule
        {
            List<barycentric2D> points;
            List<scalar> weights;
        };

        //- Static storage for quadrature rules
        static Map<quadratureRule> rules_;

    // Static Member Functions

        //- Construct rules
        static void constructRules();

        //- Choose order (the same as requested or the higher one)
        static label chooseOrder(const label requestedOrder);

        //- Return integration rules
        static const Map<quadratureRule>& rules();

public:

    // Constructors

        //- Construct from three points
        triQuadrature
        (
            const point& a,
            const point& b,
            const point& c,
            const label& order
        );

        //- Construct from a standard triangle object
        triQuadrature
        (
            const triangle<point, const point&>& tri,
            const label& order
        );

    // Member Functions

        //- Get reference to quadrature points in global Cartesian coordinates
        virtual const List<point>& points() const override;

        //- Get reference to quadrature points weights
        virtual const List<scalar>& weights() const override;

        //- Number of quadrature points
        virtual label nPoints() const override;

        //- Return number of quadrature points without making instance
        static label nPoints(label order);
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
