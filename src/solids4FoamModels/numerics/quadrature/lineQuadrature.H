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
    Gauss-Legendre quadrature for 1D interval [-1,1].

    Table data is taken from slide on this presentatin:
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

#include "triPoints.H"
#include "quadrature.H"
#include "Tuple2.H"
#include "vectorList.H"
#include "scalarList.H"
#include "Map.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lineQuadrature
:
    public quadrature,
    public linePointRef
{

    // Private data

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

        //- Return integration rules
        static const Map<quadratureRule>& rules();

public:

    // Constructors
    lineQuadrature(const linePointRef& l, const label& order);

    // Member Functions

        //- Get reference to quadrature points in global Cartesian coordinates
        virtual const List<point>& points() const;

        //- Get quadrature points in global Cartesian coordinates
        virtual const List<point> points();

        //- Get reference to quadrature points weights
        virtual const List<scalar>& weights() const;

        //- Get quadrature points weights
        virtual const List<scalar> weights();

        //- Number of quadrature points
        virtual label nPoints() const;

        //- Return number of quadrature points without making instance
        static label nPoints(label order);
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
