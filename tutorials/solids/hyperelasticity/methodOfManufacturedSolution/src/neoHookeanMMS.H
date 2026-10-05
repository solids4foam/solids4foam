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
    Foam::neoHookeanMMS

Description
    Manufactured solution for a compressible neo-Hookean solid in the total
    Lagrangian (reference configuration) form of the momentum equation:

        rho0*d2u/dt2 = div0(P) + f0

    The manufactured displacement is

        u(X, t) = a*q(X, t),
        q(X, t) = S(X)*T(t),
        S(X)    = sin(omegaX*X)*sin(omegaY*Y)*sin(omegaZ*Z),

    with the constant amplitude vector a and omega_i = n_i*pi/L, so that
    u = 0 on the boundary of the cube [0, L]^3. A zero wave number in a
    direction removes that factor from S, which gives the plane-strain
    embedding of the 2-D solution.

    The time factor T(t) is either

        ramp:   T = min(t/rampTime, 1),  d2T/dt2 = 0

    for a quasi-static load ramp, or

        cosine: T = 1 - cos(omegaT*t),   d2T/dt2 = omegaT^2*cos(omegaT*t)

    for a dynamic problem that starts from rest, i.e. u = du/dt = 0 at t = 0,
    with a nonzero initial acceleration, or

        cosineSquared: T = 0.5*(1 - cos(omegaT*t))^2,
                       d2T/dt2 = omegaT^2*(1 + c - 2*c^2),  c = cos(omegaT*t)

    for a dynamic problem whose displacement, velocity and acceleration are all
    zero at t = 0, which is consistent with the zero initial state of the
    time schemes. Both reach the peak T = 2 at omegaT*t = pi.

    The constitutive law is the neo-Hookean model with isochoric-volumetric
    split used by the solids4foam neoHookeanElastic law,

        Psi = 0.5*mu*(J^(-2/3)*I1 - 3) + 0.25*kappa*(J^2 - 1 - 2*ln(J)),

    whose stresses are, with m = mu*J^(-2/3) and
    c = 0.5*kappa*(J^2 - 1) - m*I1/3,

        P     = m*F + c*F^-T,
        sigma = (m*B + c*I)/J.

    The body force per unit reference volume that balances the manufactured
    displacement is

        f0 = rho0*a*d2q/dt2 - div0(P(u)),

    evaluated in closed form from the rank-one structure of grad0(u).
    The derivation is given in the tutorial documentation.

    All parameters are read from constant/neoHookeanMMSProperties.

Author
    Ivan Batistic, FSB Zagreb. Derivation by Pablo Castrillo.

SourceFiles
    neoHookeanMMS.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMMS_H
#define neoHookeanMMS_H

#include "fvMesh.H"
#include "volFields.H"
#include "symmTensor.H"
#include "tensor.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class neoHookeanMMS Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMMS
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Amplitude vector
        const vector a_;

        //- Spatial wave numbers omega_i = n_i*pi/L
        const vector omega_;

        //- Time function: ramp, cosine or cosineSquared
        const word timeFunction_;

        //- Ramp duration for the ramp time function
        const scalar rampTime_;

        //- Angular frequency for the cosine time function
        const scalar omegaT_;

        //- Reference density
        const scalar rho_;

        //- Shear modulus
        const scalar mu_;

        //- Bulk modulus
        const scalar kappa_;

        //- Body force field at the cell centres
        mutable autoPtr<volVectorField> bodyForcesPtr_;

        //- Time index at which the body force field was last evaluated
        mutable label bodyForcesTimeIndex_;


    // Private Member Functions

        //- Return the time factor and its first and second time derivatives
        void timeFactor
        (
            const scalar t,
            scalar& T,
            scalar& dTdt,
            scalar& d2Tdt2
        ) const;

        //- Return the spatial factor S(X)
        scalar spatialFactor(const vector& X) const;

        //- Evaluate q, g = grad0(q), H = grad0(grad0(q)), the Laplacian of q
        //  and d2q/dt2 at a point and time
        void evaluate
        (
            const vector& X,
            const scalar t,
            scalar& q,
            vector& g,
            symmTensor& H,
            scalar& lapq,
            scalar& d2qdt2
        ) const;

        //- Evaluate the body force field at the cell centres
        void calcBodyForces() const;

        //- Disallow default bitwise copy construct
        neoHookeanMMS(const neoHookeanMMS&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanMMS&);


public:

    //- Runtime type information
    TypeName("neoHookeanMMS");


    // Constructors

        //- Construct from the mesh and a dictionary
        neoHookeanMMS(const fvMesh& mesh, const dictionary& dict);

        //- Construct from the mesh, reading constant/neoHookeanMMSProperties
        explicit neoHookeanMMS(const fvMesh& mesh);


    //- Destructor
    virtual ~neoHookeanMMS()
    {}


    // Member Functions

        //- Return the manufactured displacement at a point and time
        vector displacement(const vector& X, const scalar t) const;

        //- Return the manufactured velocity at a point and time
        vector velocity(const vector& X, const scalar t) const;

        //- Return the deformation gradient at a point and time
        tensor deformationGradient(const vector& X, const scalar t) const;

        //- Return the Cauchy stress at a point and time
        symmTensor cauchyStress(const vector& X, const scalar t) const;

        //- Return the body force per unit reference volume at a point and
        //  time
        vector bodyForce(const vector& X, const scalar t) const;

        //- Return the body force per unit reference volume at the cell
        //  centres for the current time
        const volVectorField& bodyForces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
