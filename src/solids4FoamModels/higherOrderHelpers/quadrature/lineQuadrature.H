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
    lineQuadrature

Description
    Numerical quadrature rules for integration over a one-dimensional
    reference line element.

    This class implements Gaussian quadrature rules on the unit line,
    providing quadrature points and weights suitable for exact
    integration of polynomials up to a specified order. The quadrature
    points are defined in the local coordinate system of the reference
    line element.

    The lineQuadrature class is typically used for edge-based
    integrations, one-dimensional projections, and lower-dimensional
    integration tasks arising in high-order finite-volume and
    finite-element formulations.

    Gauss-Legendre quadrature for 1D interval [-1,1].

    Table data is taken from slide on this presentation:
    https://www.to.infn.it/~mignone/Numerical_Algorithms/
    ch03_numerical_quadrature.pdf

SourceFile
    lineQuadrature.C

Author
    Ivan Batistic, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef lineQuadrature_H
#define lineQuadrature_H

#include "quadrature.H"
#include "Map.H"
#include "line.H"
#include "Field.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lineQuadrature
:
    public quadrature
{

    // Private data

        //- Store triangle using triangle template
        line<point, point> line_;

        //- Quadrature points in Cartesian coordinates
        List<point> points_;

        //- Quadrature points weights
        List<scalar> weights_;

    // Private Member Functions

        //- Map integration points from barycentric to cartesian coordinates
        tmp<Field<point>> parametricToPoint
        (
            const List<scalar>& x
        ) const;

        //- Normalise weights to make sum of weights equal to 1
        tmp<Field<scalar>> normaliseWeights
        (
            const List<scalar>& w
        ) const;

    // Static data

        //- Maximum order of accuracy when integrating polynomials
        static constexpr label maxSupportedOrder = 19;

        //- Quadrature rules storage
        struct quadratureRule
        {
            List<scalar> points;
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

        //- Construct from two points
        lineQuadrature
        (
            const point& start,
            const point& end,
            const label& order
        );

        //- Construct from a standard line object
        lineQuadrature
        (
            const line<point,
            const point&>& l,
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
