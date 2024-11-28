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
    triangleQuadrature

Description
    Gauss quadrature for triangle

    Source:
    D. A. Dunavant,
    High degree efficient symmetrical Gaussian quadrature rule for the triangle.
    International Journal for Numerical Methods in Engineering,
    21(6):1129–1148, jun 1985.

SourceFile
    triangleQuadrature.C

Author
    Ivan Batistic, FAMENA. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef triangleQuadrature_H
#define triangleQuadrature_H

#include "triPoints.H"
#include "Tuple2.H"
#include "vectorList.H"
#include "scalarList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triangleQuadrature
:
    public triPoints
{
    // Private data

        //- Quadrature rules storage
        struct quadratureRule
        {
            List<vector> gaussPoints;
            List<scalar> weights;
        };

        //- Static storage for quadrature rules
        static Map<quadratureRule> rules_;

        //- Construct rules
        static void constructRules();

        //- Number of Gauss points
        const label n_;

        //- Gaussian quadrature points in Cartesian coordinates
        List<point> gaussPoints_;

        //- Gaussian quadrature points weights
        List<scalar> weights_;

    // Private Member Functions

        //- Map integration points from barycentric to cartesian coordinates
        List<point> barycentricToPoint(const List<point>& localGP);

public:

    // Constructors
    triangleQuadrature(const triPoints& pts, const label& n);

    // Member Functions

    // Access

        //- Predefined integration rules
        static const Map<quadratureRule>& rules();

        //- Get reference to Gaussian quadrature points in global Cartesian
        //  coordinates
        const List<point>& gaussPoints() const;

        //- Get reference to Gaussian quadrature points weights
        const List<scalar>& weights() const;

        //- Get Gaussian quadrature points in global Cartesian coordinates
        const List<point> gaussPoints();

        //- Get Gaussian quadrature points weights
        const List<scalar> weights();
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
