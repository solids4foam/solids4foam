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

#include "squarePlateStressDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "interpolationTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


Foam::scalar Foam::allEdgesClampedMaxDisplacement
(
    const vector& C, // Position vector
    const scalar p,  // Pressure applied
    const scalar E,  // Young's modulus
    const scalar nu, // Poisson's ratio
    const scalar a,  // Length of the plate
    const scalar b,  // Width of the plate
    const scalar h   // Thickness of the plate
)
{
    // Plate bending stiffness
    const scalar D = E*Foam::pow(h, 3)/(12*(1 - Foam::pow(nu, 2)));

    const scalar beta = b / a;

    // Coefficients from the literature
    List<Tuple2<scalar, scalar>> coeffsList;

    coeffsList.append(Tuple2<scalar, scalar>(1.0, 0.001265319087));
    coeffsList.append(Tuple2<scalar, scalar>(1.2, 0.001724870503));
    coeffsList.append(Tuple2<scalar, scalar>(1.4, 0.002068143209));
    coeffsList.append(Tuple2<scalar, scalar>(1.6, 0.002299966977));
    coeffsList.append(Tuple2<scalar, scalar>(1.8, 0.002446162656));
    coeffsList.append(Tuple2<scalar, scalar>(2.0, 0.002532955769));
    coeffsList.append(Tuple2<scalar, scalar>(20.0, 0.002604166667));

    interpolationTable<scalar> Ctable
    (
        coeffsList,
        interpolationTable<scalar>::CLAMP,
        ""
    );

    const scalar coeff = Ctable(beta);

    return coeff*p*Foam::pow(a, 4) / D;
}



Foam::scalar Foam::allEdgesSupportedDisplacement
(
    const vector& C, // Position vector
    const scalar p,  // Pressure applied
    const scalar E,  // Young's modulus
    const scalar nu, // Poisson's ratio
    const scalar a,  // Length of the plate
    const scalar b,  // Width of the plate
    const scalar h   // Thickness of the plate
)
{
    // Take x and y coordinates
    const scalar x = C.x();
    const scalar y = C.y();

    // Plate bending stiffness
    const scalar D = E*Foam::pow(h, 3)/(12*(1 - Foam::pow(nu, 2)));

#ifdef OPENFOAM_NOT_EXTEND
    const scalar pi = constant::mathematical::pi;
#else
    const scalar pi = mathematicalConstant::pi;
#endif

    // Number of terms in series
    const label N = 20;

    scalar sum = 0.0;
    for (label m=1; m<= N*2; m+=2)
    {
        const scalar alphaM = (m*pi*b)/(2*a);
        const scalar t2 =
            ((alphaM * Foam::tanh(alphaM) + 2.0) / (2.0 * Foam::cosh(alphaM)))
          * Foam::cosh(m*pi*y/a);

        const scalar t3 =
            (1.0 / (2.0 * Foam::cosh(alphaM))) * (m*pi*y/a)
          * Foam::sinh(m*pi*y/a);

        const scalar t =
            (Foam::pow(-1, 0.5*(m-1)) / Foam::pow(m, 5))*Foam::cos(m*pi*x/a);

        sum += t*(1.0-t2+t3);
    }

    return ((4*p*Foam::pow(a, 4)) / (Foam::pow(pi, 5)*D))*sum;
}

// ************************************************************************* //
