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
    Gaussian quadrature integration for triangle

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triangleQuadrature
:
    public triPoints
{
    // Private data

        //- Number of Gauss points
        const label n_;

        //- Mapped Gauss points and corresponding weights
        List<Tuple2<point, scalar>> gaussPointsWeights_;

    // Private Member Functions

        //- Map integration points from local natural to cartesian coordinates
        void mapGaussPoints();

public:

    // Constructors
    triangleQuadrature(const triPoints& pts, const label& n);

    // Member Functions

    // Properties

        //- Triangle quadrature points and weights in local natural coordinates
        static const List<Tuple2
            <
                label,
                Tuple2<vectorList, scalarList>
            >> TRI_QUADRATURE;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
