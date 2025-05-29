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

#include "sphericalCavityStressDisplacement.H"
#include "coordinateSystem.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::KroneckerDelta(const int i, const int j)
{
    return (i == j) ? 1.0 : 0.0;
}

Foam::scalar Foam::sphericalCavityDisplacementComponent
(
    const int i,
    const scalar nu,
    const scalar sigma_0,
    const scalar E,
    const scalar a,
    const vector& x
)
{
    const scalar R = mag(x);
    const scalar x_3 = x[vector::Z];

    const scalar term1 =
        2
      + (5*(5 - 4*nu)*pow(a, 3))/((7 - 5*nu)*pow(R, 3))
      + (6*pow(a, 5))/((7 - 5*nu)*pow(R, 5));
    const scalar term2 =
        (-2*nu)/(1 + nu)
      + (5*nu - 6)*pow(a, 3)/((7 - 5*nu)*pow(R, 3))
      + 3*pow(a, 5)/((7 - 5*nu)*pow(R, 5));
    const scalar common_factor = (1 + nu)*sigma_0/(2*E);

    return
        common_factor
       *(
           term1*x_3*KroneckerDelta(i, 2)
         + term2*x[i]*(1 - 5*pow(x_3, 2)/pow(R, 2))
       );
}


Foam::vector Foam::sphericalCavityDisplacement
(
    const scalar nu,
    const scalar T,
    const scalar E,
    const scalar a,
    const vector& x
)
{
    // Bower 2009 approach - not correct
    // vector disp = vector::zero;
    // for (int i = 0; i < 3; ++i)
    // {
    //     disp[i] =
    //         sphericalCavityDisplacementComponent(i, nu, T, E, a, x);
    // }

    // Goodier 1933

    // Shear modulus
    const scalar mu = E/(2*(1 + nu));

    // Constants for a spherical cavity
    const scalar A = -pow3(a)*(T/(8*mu))*(13 - 10*nu)/(7 - 5*nu);
    const scalar B = pow5(a)*(T/(8*mu))*(1/(7 - 5*nu));
    const scalar C = pow3(a)*(T/(8*mu))*((5*(1 - 2*nu))/(7 - 5*nu));

    // Calculate the spherical polar coordinates
    const scalar r = mag(x);
    if (r < SMALL)
    {
        FatalError
            << "r is less than SMALL. "
            << "Are you sure the cavity is centred on the origin?"
            << abort(FatalError);
    }
    const scalar theta = Foam::acos(x.z()/r);
    const scalar phi = Foam::atan2(x.y(), x.x());

    // Radial displacement
    const scalar ur =
       - A/sqr(r) - 3*B/pow4(r)
       + (
             ((5 - 4*nu)/(1 - 2*nu))*C/sqr(r) - 9*B/pow4(r)
         )*cos(2*theta);

    // Theta displacement
    const scalar utheta = -(2*C/sqr(r) + 6*B/pow4(r))*sin(2*theta);

    // Spherical polar displacement vector (r, theta, phi)
    const vector dispPolar(ur, utheta, 0);

    // Calculate the Cartesian displacement vector
    const vector dispCart
    (
        ur*sin(theta)*cos(phi) + utheta*cos(theta)*cos(phi),
        ur*sin(theta)*sin(phi) + utheta*cos(theta)*sin(phi),
        ur*cos(theta) - utheta*sin(theta)
    );

    // Calculate the displacement field called by a uniform uniaxial T when
    // there is no cavity
    const vector dispUni(-nu*T*x.x()/E, -nu*T*x.y()/E, T*x.z()/E);

    // Return the total Cartesian displacement vector
    return dispCart + dispUni;
}


Foam::symmTensor Foam::sphericalCavityStress
(
    const scalar T,
    const scalar nu,
    const scalar a,
    const vector& x
)
{
    // Stress formulae from Southwell et al, 1926, On the concentration of
    // stress in the neighbourhood of a small spherical flaw; and on the
    // propagation of fatigue fractures in Statistically Isotropic materials

    // Polar radial coordinate
    // Take care: in general, R != r
    const scalar R = mag(x);

    // Cylindrical radial coordinate
    const scalar r = sqrt(sqr(x.x()) + sqr(x.y()));

    // Cylindrical axial coordinate
    const scalar z = x.z();

    // Initialise stress tensor
    // We will place rr in xx, thetaTheta in yy, zz in zz, and zr in xz
    // We will late rotate this stress tensor to the Cartesian coordinate system
    symmTensor sigmaCylindrical = symmTensor::zero;

    // Radial (rr) component
    sigmaCylindrical.xx() =
        (T/(14 - 10*nu))*(pow3(a)/pow3(R))
       *(
            9 - 15*nu - 12*(sqr(a)/sqr(R))
          - (sqr(r)/sqr(R))*(72 - 15*nu - 105*(sqr(a)/sqr(R)))
          + 15*(pow4(r)/pow4(R))*(5 - 7*(sqr(a)/sqr(R)))
       );

    // Hoop (thetaTheta) component
    sigmaCylindrical.yy() =
        (T/(14 - 10*nu))*(pow3(a)/pow3(R))
       *(
            9 - 15*nu - 12*(sqr(a)/sqr(R))
          - 15*(sqr(r)/sqr(R))*(1 - 2*nu - (sqr(a)/sqr(R)))
       );

    // Axial (zz) component
    sigmaCylindrical.zz() =
        T
       *(
            1
          - (1/(14 - 10*nu))*(pow3(a)/pow3(R))
           *(
                38 - 10*nu - 24*sqr(a)/sqr(R)
              - (sqr(r)/sqr(R))*(117 - 15*nu - 120*sqr(a)/sqr(R))
              + 15*(pow4(r)/pow4(R))*(5 - 7*sqr(a)/sqr(R))
            )
        );

    // Radial-axial shear (rz) component
    sigmaCylindrical.xz() =
        (T/(14 - 10*nu))*(pow3(a)*z*r/pow5(R))
       *(
          - 3*(19 - 5*nu) + 60*(sqr(a)/sqr(R))
          + 15*(sqr(r)/sqr(R))*(5 - 7*(sqr(a)/sqr(R)))
        );

    // Rotate the cylindrical stress to Cartesian coordinates
    const coordinateSystem cs("cylindrical", x, vector(0, 0, 1), x/mag(x));
#ifdef OPENFOAM_COM
    const symmTensor sigmaCartesian = transform(cs.R(), sigmaCylindrical);
#else
    const symmTensor sigmaCartesian = transform(cs.R().R(), sigmaCylindrical);
#endif

    return sigmaCartesian;
}


// ************************************************************************* //
