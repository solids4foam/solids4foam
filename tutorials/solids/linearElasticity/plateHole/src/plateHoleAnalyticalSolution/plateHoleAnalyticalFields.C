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

\*---------------------------------------------------------------------------*/

#include "plateHoleAnalyticalFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::symmTensor Foam::plateHoleAnalyticalFields::stress
(
    const vector& C,
    const scalar farFieldTractionX,
    const scalar holeRadius
)
{
    symmTensor sigma = symmTensor::zero;

    const scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    const scalar theta = atan2(C.y(), C.x());

    sigma.xx() =
        farFieldTractionX
       *(
           1.0
         - (sqr(holeRadius)/sqr(r))*(3*cos(2*theta)/2 + cos(4*theta))
         + (3*pow(holeRadius, 4)/(2*pow(r, 4)))*cos(4*theta)
        );

    sigma.yy() =
        farFieldTractionX
       *(
         - (sqr(holeRadius)/sqr(r))*(cos(2*theta)/2 - cos(4*theta))
         - (3*pow(holeRadius, 4)/(2*pow(r, 4)))*cos(4*theta)
        );

    sigma.xy() =
        farFieldTractionX
       *(
         - (sqr(holeRadius)/sqr(r))*(sin(2*theta)/2 + sin(4*theta))
         + (3*pow(holeRadius, 4)/(2*pow(r, 4)))*sin(4*theta)
        );

    return sigma;
}


Foam::vector Foam::plateHoleAnalyticalFields::displacement
(
    const vector& C,
    const scalar farFieldTractionX,
    const scalar holeRadius,
    const scalar mu,
    const scalar kappa
)
{
    const scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    const scalar theta = atan2(C.y(), C.x());

    return vector
    (
        (holeRadius*farFieldTractionX/(8*mu))
       *(
           (r/holeRadius)*(kappa + 1)*cos(theta)
         + (2*holeRadius/r)*((1 + kappa)*cos(theta) + cos(3*theta))
         - (2*pow(holeRadius, 3)/pow(r, 3))*cos(3*theta)
        ),
        (holeRadius*farFieldTractionX/(8*mu))
       *(
           (r/holeRadius)*(kappa - 3)*sin(theta)
         + (2*holeRadius/r)*((1 - kappa)*sin(theta) + sin(3*theta))
         - (2*pow(holeRadius, 3)/pow(r, 3))*sin(3*theta)
        ),
        0.0
    );
}


Foam::vector Foam::plateHoleAnalyticalFields::displacement
(
    const vector& C,
    const scalar farFieldTractionX,
    const scalar holeRadius,
    const scalar E,
    const scalar nu,
    const bool planeStress
)
{
    const scalar mu = E/(2*(1 + nu));
    scalar kappa = 3 - 4*nu;

    if (planeStress)
    {
        kappa = (3.0 - nu)/(1.0 + nu);
    }

    return displacement(C, farFieldTractionX, holeRadius, mu, kappa);
}


Foam::scalar Foam::plateHoleAnalyticalFields::hydPressure
(
    const vector& C,
    const scalar farFieldTractionX,
    const scalar holeRadius,
    const scalar nu
)
{
    const symmTensor sigma = stress(C, farFieldTractionX, holeRadius);
    const scalar sigmaZZ = nu*(sigma.xx() + sigma.yy());

    return -(sigma.xx() + sigma.yy() + sigmaZZ)/3;
}


Foam::vector Foam::plateHoleAnalyticalFields::traction
(
    const vector& C,
    const vector& n,
    const scalar farFieldTractionX,
    const scalar holeRadius
)
{
    return n & stress(C, farFieldTractionX, holeRadius);
}


// ************************************************************************* //
