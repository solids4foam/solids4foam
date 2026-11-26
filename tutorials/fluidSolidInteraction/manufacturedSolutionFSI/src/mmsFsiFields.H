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

Namespace
    Foam::mmsFsiFields

Description
    Evaluate the manufactured FSI solution of the generated header
    mmsFsiExact.H (see ../derivation) at OpenFOAM points.

    The fluid pressure, pressure gradient and momentum source are kinematic,
    i.e. divided by the fluid density, as used by the solids4foam fluid
    models. The solid body force is per unit reference volume.

\*---------------------------------------------------------------------------*/

#ifndef mmsFsiFields_H
#define mmsFsiFields_H

#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "mmsFsiExact.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mmsFsiFields
{

//- Exact fluid state at a point of the current configuration
struct fluidState
{
    //- Velocity
    vector U;

    //- Kinematic pressure
    scalar p;

    //- Momentum source per unit mass
    vector f;

    //- Velocity gradient, gradU(i, j) = dU_j/dx_i
    tensor gradU;

    //- Kinematic pressure gradient
    vector gradp;
};


//- Exact solid state at a point of the reference configuration
struct solidState
{
    //- Displacement
    vector D;

    //- Velocity
    vector U;

    //- Body force per unit reference volume
    vector b0;

    //- Cauchy stress
    symmTensor sigma;
};


//- Return the exact fluid state at point c and time t
inline fluidState fluid(const vector& c, const scalar t)
{
    double U[2];
    double p[1];
    double f[2];
    double G[4];
    double gp[2];
    mmsFsi::fluid(c.x(), c.y(), t, U, p, f, G, gp);

    // G[2*i + j] = dU_i/dx_j, transposed to the OpenFOAM convention
    fluidState s;
    s.U = vector(U[0], U[1], 0);
    s.p = p[0];
    s.f = vector(f[0], f[1], 0);
    s.gradU = tensor(G[0], G[2], 0, G[1], G[3], 0, 0, 0, 0);
    s.gradp = vector(gp[0], gp[1], 0);

    return s;
}


//- Return the exact solid state at reference point C and time t
inline solidState solid(const vector& C, const scalar t)
{
    double D[2];
    double U[2];
    double b0[2];
    double sigma[6];
    mmsFsi::solid(C.x(), C.y(), t, D, U, b0, sigma);

    solidState s;
    s.D = vector(D[0], D[1], 0);
    s.U = vector(U[0], U[1], 0);
    s.b0 = vector(b0[0], b0[1], 0);
    s.sigma =
        symmTensor
        (
            sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5]
        );

    return s;
}


//- Return the exact interface data at reference abscissa X and time t: the
//  deformed position x, the fluid outward unit normal n and the Cauchy
//  traction tr = sigma & n, so the traction applied to the solid is -tr
inline void interface
(
    const scalar X,
    const scalar t,
    vector& x,
    vector& tr,
    vector& n
)
{
    double xx[2];
    double tt[2];
    double nn[2];
    mmsFsi::interface(X, t, xx, tt, nn);

    x = vector(xx[0], xx[1], 0);
    tr = vector(tt[0], tt[1], 0);
    n = vector(nn[0], nn[1], 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mmsFsiFields
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
