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
    tetQuadrature

Description
    Numerical quadrature rules for integration over a three-dimensional
    reference tetrahedral element.

    This class implements Gaussian-type quadrature rules for tetrahedra,
    providing quadrature points in the local reference coordinate system
    together with the corresponding weights. The quadrature order
    controls the degree of polynomial fields that can be integrated
    accurately.

    The tetQuadrature class is primarily intended for volume
    integrations in three-dimensional finite-volume and finite-element
    formulations, including high-order reconstruction and weak-form
    assembly.

    Source:
    Lee Shunn, Frank Ham,
    Symmetric quadrature rules for tetrahedra based on a cubic close-packed
    lattice arrangement, Journal of Computational and Applied Mathematics, 2012

SourceFile
    tetQuadrature.C

Author
    Ivan Batistic, UCD. All rights reserved.
    Pablo Castrillo, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef tetQuadrature_H
#define tetQuadrature_H

#include "tetrahedron.H"
#include "barycentric.H"
#include "quadrature.H"
#include "Map.H"
#include "Field.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tetQuadrature
:
    public quadrature
{

    // Private data

        //- Store tetrahedron using tetrahedron template
        tetrahedron<point, point> tet_;

        //- Quadrature points in Cartesian coordinates
        List<point> points_;

        //- Quadrature points weights
        List<scalar> weights_;

    // Private Member Functions

        //- Map integration points from barycentric to cartesian coordinates
        tmp<Field<point>> barycentricToPoint
        (
            const List<barycentric>& localPts
        ) const;

    // Static data

        //- Maximum order of accuracy when integrating polynomials
        static constexpr label maxSupportedOrder = 6;

        //- Quadrature rules storage
        struct quadratureRule
        {
            List<barycentric> points;
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

        //- Construct from four points
        tetQuadrature
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d,
            const label& order
        );

        //- Construct from a standard tetrahedron object
        tetQuadrature
        (
            const tetrahedron<point, const point&>& tet,
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
