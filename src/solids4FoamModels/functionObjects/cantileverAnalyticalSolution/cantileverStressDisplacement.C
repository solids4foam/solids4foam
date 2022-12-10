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

\*----------------------------------------------------------------------------*/

#include "cantileverStressDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::symmTensor Foam::cantileverStress
(
    const vector& C, // Position vector
    const scalar P, // Load applied in the minus y direction
    const scalar E, // Young's modulus
    const scalar nu, // Poisson's ratio
    const scalar L, // Length of the beam
    const scalar D, // Depth of the beam
    const scalar I // Second moment of area
)
{
    // Initialise to zero
    symmTensor sigma = symmTensor::zero;

    // Take x and y coordinates
    const scalar x = C.x();
    const scalar y = C.y();

    // Axial stress
    sigma.xx() = P*(L - x)*y/I;

    // Shear stress
    sigma.xy() = -(P/(2*I))*(D*D/4 - y*y);

    return sigma;
}


Foam::vector Foam::cantileverDisplacement
(
    const vector& C, // Position vector
    const scalar P, // Load applied in the minus y direction
    const scalar E, // Young's modulus
    const scalar nu, // Poisson's ratio
    const scalar L, // Length of the beam
    const scalar D, // Depth of the beam
    const scalar I // Second moment of area
)
{
    // Initialise to zero
    vector disp = vector::zero;

    // Augarde and Deeks

    // Take x and y coordinates
    const scalar x = C.x();
    const scalar y = C.y();

    // X displacement
    disp.x() = (P*y/(6*E*I))*((6*L - 3*x)*x + (2 + nu)*(y*y - D*D/4));

    // Y displacement
    disp.y() =
        -(P/(6*E*I))
        *(
            3*nu*y*y*(L - x) + (4 + 5*nu)*D*D*x/4 + (3*L - x)*x*x
        );

    // Timoshenko and Goodier
    // Note: I assume origin is at centre of fixed end and positive y points up
    // Take x and y coordinates
    // const scalar x = L - C.x();
    // const scalar y = -C.y();

    // // Half the depth
    // const scalar c = D/2.0;

    // // Shear modulus
    // const scalar G = E/(2*(1 + nu));

    // // X displacement
    // disp.x() = -P*
    // (
    //   - (x*x*y/2/E)
    //   - (nu*y*y*y/6/E)
    //   + (y*y*y/6/G)
    //   + (L*L/2/E - c*c/2/G)*y
    // )/I;

    // // Y displacement
    // disp.y() = -P*
    // (
    //     (nu*x*y*y/2)
    //   + (x*x*x/6)
    //   - (L*L*x/2)
    //   + (L*L*L/3)
    // )/E/I;

    return disp;
}


// ************************************************************************* //
