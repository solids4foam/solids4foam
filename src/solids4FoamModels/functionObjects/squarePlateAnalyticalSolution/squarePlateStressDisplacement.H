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
    squarePlateAnalyticalSolution

Description
    Functions that calculate the analytical deflection for a
    thin square plate with uniform transverse pressure.

    The solution is taken from:
    Timoshenko, S., & Woinowsky-Krieger, S. (1959). Theory of plates and shells.

    Solution for clamped and simply supported page can be found on page 202
    and 197.

    The maximal deflection of a clamped plate in Timoshenko’s book is attributed
    to T. H. Evans (Journal of Applied Mechanics, vol. 6, p. A-7, 1939).
    However, because the accuracy of this coefficient is limited,
    more precise coefficient values are taken from the following sources:

    Taylor, Robert L., and Sanjay Govindjee. “Solution of clamped rectangular
    plate problems.” Communications in Numerical Methods in Engineering 20.10
    (2004): 757–765.

    Imrak, C. E., and I. Gerdemeli. “An exact solution for the deflection of
    a clamped rectangular plate under uniform load.” Applied Mathematical
    Sciences 1.43 (2007): 2129–2137.

    It is assumed that:
    - the load (pressure) is uniform
    - the plate dimensions are axb and plate is d deep/thick
    - the origin is at the centre of plate

Author
    Ivan Batistic & Philip Cardiff, UCD. All rights reserved.

SourceFiles
    squarePlateStressDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef squarePlateStressDisplacement_H
#define squarePlateStressDisplacement_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    //- Analytical displacement solution (max value)
    scalar allEdgesClampedMaxDisplacement
    (
        const vector& C, // Position vector (not used)
        const scalar p,  // Pressure
        const scalar E,  // Young's modulus
        const scalar nu, // Poisson's ratio
        const scalar a,  // Length of the plate
        const scalar b,  // Width of the plate
        const scalar h   // Thickness of the plate
    );

    //- Analytical displacemen solution
    scalar allEdgesSupportedDisplacement
    (
        const vector& C, // Position vector
        const scalar p,  // Pressure
        const scalar E,  // Young's modulus
        const scalar nu, // Poisson's ratio
        const scalar a,  // Length of the plate
        const scalar b,  // Width of the plate
        const scalar h   // Thickness of the plate
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
