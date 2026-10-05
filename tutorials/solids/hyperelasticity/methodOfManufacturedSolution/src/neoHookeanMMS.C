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

#include "neoHookeanMMS.H"
#include "IOdictionary.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanMMS, 0);
}


namespace
{

Foam::dictionary readPropertiesDict(const Foam::fvMesh& mesh)
{
    // Not registered, so that each user of the solution may read it
    return Foam::IOdictionary
    (
        Foam::IOobject
        (
            "neoHookeanMMSProperties",
            mesh.time().constant(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE,
            false
        )
    );
}


Foam::vector readOmega(const Foam::dictionary& dict)
{
#ifdef FOAMEXTEND
    const Foam::scalar pi = Foam::mathematicalConstant::pi;
#else
    const Foam::scalar pi = Foam::constant::mathematical::pi;
#endif

    const Foam::scalar L = Foam::readScalar(dict.lookup("L"));
    const Foam::vector nWaves(dict.lookup("nWaves"));

    if (L < Foam::SMALL)
    {
        FatalIOErrorIn("readOmega(const dictionary&)", dict)
            << "L must be positive" << Foam::exit(Foam::FatalIOError);
    }

    return nWaves*pi/L;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::neoHookeanMMS::timeFactor
(
    const scalar t,
    scalar& T,
    scalar& dTdt,
    scalar& d2Tdt2
) const
{
    if (timeFunction_ == "ramp")
    {
        T = min(t/rampTime_, 1.0);
        dTdt = (t < rampTime_) ? 1.0/rampTime_ : 0.0;
        d2Tdt2 = 0.0;
    }
    else if (timeFunction_ == "cosine")
    {
        T = 1.0 - Foam::cos(omegaT_*t);
        dTdt = omegaT_*Foam::sin(omegaT_*t);
        d2Tdt2 = sqr(omegaT_)*Foam::cos(omegaT_*t);
    }
    else
    {
        // cosineSquared: zero displacement, velocity and acceleration at t = 0
        const scalar c = Foam::cos(omegaT_*t);
        const scalar s = Foam::sin(omegaT_*t);
        T = 0.5*sqr(1.0 - c);
        dTdt = omegaT_*(1.0 - c)*s;
        d2Tdt2 = sqr(omegaT_)*(1.0 + c - 2.0*sqr(c));
    }
}


Foam::scalar Foam::neoHookeanMMS::spatialFactor(const vector& X) const
{
    scalar S = 1.0;
    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        if (mag(omega_[cmpt]) > SMALL)
        {
            S *= Foam::sin(omega_[cmpt]*X[cmpt]);
        }
    }

    return S;
}


void Foam::neoHookeanMMS::evaluate
(
    const vector& X,
    const scalar t,
    scalar& q,
    vector& g,
    symmTensor& H,
    scalar& lapq,
    scalar& d2qdt2
) const
{
    // Sines and cosines in each direction. A zero wave number removes the
    // factor from S, i.e. sin -> 1 and the derivative -> 0
    vector s = vector::one;
    vector c = vector::zero;
    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        if (mag(omega_[cmpt]) > SMALL)
        {
            s[cmpt] = Foam::sin(omega_[cmpt]*X[cmpt]);
            c[cmpt] = Foam::cos(omega_[cmpt]*X[cmpt]);
        }
    }

    scalar T = 0;
    scalar dTdt = 0;
    scalar d2Tdt2 = 0;
    timeFactor(t, T, dTdt, d2Tdt2);

    const scalar S = s.x()*s.y()*s.z();

    q = S*T;
    d2qdt2 = S*d2Tdt2;

    g = T*vector
    (
        omega_.x()*c.x()*s.y()*s.z(),
        omega_.y()*s.x()*c.y()*s.z(),
        omega_.z()*s.x()*s.y()*c.z()
    );

    H = symmTensor
    (
       -sqr(omega_.x())*q,
        T*omega_.x()*omega_.y()*c.x()*c.y()*s.z(),
        T*omega_.x()*omega_.z()*c.x()*s.y()*c.z(),
       -sqr(omega_.y())*q,
        T*omega_.y()*omega_.z()*s.x()*c.y()*c.z(),
       -sqr(omega_.z())*q
    );

    lapq = -(omega_ & omega_)*q;
}


void Foam::neoHookeanMMS::calcBodyForces() const
{
    if (!bodyForcesPtr_.valid())
    {
        bodyForcesPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "neoHookeanMMSBodyForce",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("zero", dimForce/dimVolume, vector::zero),
                "zeroGradient"
            )
        );
    }

    const scalar t = mesh_.time().value();
    const vectorField& C = mesh_.C();
    vectorField& bodyForcesI = bodyForcesPtr_();

    forAll(bodyForcesI, cellI)
    {
        bodyForcesI[cellI] = bodyForce(C[cellI], t);
    }

    bodyForcesPtr_().correctBoundaryConditions();

    bodyForcesTimeIndex_ = mesh_.time().timeIndex();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neoHookeanMMS::neoHookeanMMS
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    a_(dict.lookup("a")),
    omega_(readOmega(dict)),
    timeFunction_(dict.lookup("timeFunction")),
    rampTime_(dict.lookupOrDefault<scalar>("rampTime", 1.0)),
    omegaT_(dict.lookupOrDefault<scalar>("omegaT", 0.0)),
    rho_(readScalar(dict.lookup("rho"))),
    mu_(0.0),
    kappa_(0.0),
    bodyForcesPtr_(),
    bodyForcesTimeIndex_(-1)
{
    const scalar E = readScalar(dict.lookup("E"));
    const scalar nu = readScalar(dict.lookup("nu"));

    if (E < SMALL || nu < 0.0 || nu >= 0.5)
    {
        FatalIOErrorIn("neoHookeanMMS::neoHookeanMMS(...)", dict)
            << "E must be positive and 0 <= nu < 0.5"
            << exit(FatalIOError);
    }

    const_cast<scalar&>(mu_) = E/(2.0*(1.0 + nu));
    const_cast<scalar&>(kappa_) = E/(3.0*(1.0 - 2.0*nu));

    if (timeFunction_ == "ramp")
    {
        if (rampTime_ < SMALL)
        {
            FatalIOErrorIn("neoHookeanMMS::neoHookeanMMS(...)", dict)
                << "rampTime must be positive" << exit(FatalIOError);
        }
    }
    else if (timeFunction_ == "cosine" || timeFunction_ == "cosineSquared")
    {
        if (omegaT_ < SMALL)
        {
            FatalIOErrorIn("neoHookeanMMS::neoHookeanMMS(...)", dict)
                << "omegaT must be positive" << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorIn("neoHookeanMMS::neoHookeanMMS(...)", dict)
            << "Unknown timeFunction " << timeFunction_
            << ". Valid options are: ramp cosine cosineSquared"
            << exit(FatalIOError);
    }
}


Foam::neoHookeanMMS::neoHookeanMMS(const fvMesh& mesh)
:
    neoHookeanMMS(mesh, readPropertiesDict(mesh))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::neoHookeanMMS::displacement
(
    const vector& X,
    const scalar t
) const
{
    scalar q = 0;
    vector g = vector::zero;
    symmTensor H = symmTensor::zero;
    scalar lapq = 0;
    scalar d2qdt2 = 0;
    evaluate(X, t, q, g, H, lapq, d2qdt2);

    return a_*q;
}


Foam::vector Foam::neoHookeanMMS::velocity
(
    const vector& X,
    const scalar t
) const
{
    scalar T = 0;
    scalar dTdt = 0;
    scalar d2Tdt2 = 0;
    timeFactor(t, T, dTdt, d2Tdt2);

    return a_*spatialFactor(X)*dTdt;
}


Foam::tensor Foam::neoHookeanMMS::deformationGradient
(
    const vector& X,
    const scalar t
) const
{
    scalar q = 0;
    vector g = vector::zero;
    symmTensor H = symmTensor::zero;
    scalar lapq = 0;
    scalar d2qdt2 = 0;
    evaluate(X, t, q, g, H, lapq, d2qdt2);

    return tensor(I) + (a_*g);
}


Foam::symmTensor Foam::neoHookeanMMS::cauchyStress
(
    const vector& X,
    const scalar t
) const
{
    scalar q = 0;
    vector g = vector::zero;
    symmTensor H = symmTensor::zero;
    scalar lapq = 0;
    scalar d2qdt2 = 0;
    evaluate(X, t, q, g, H, lapq, d2qdt2);

    const scalar J = 1.0 + (a_ & g);
    const scalar A = a_ & a_;
    const scalar G = g & g;
    const scalar I1 = 3.0 + 2.0*(a_ & g) + A*G;

    const scalar m = mu_*Foam::pow(J, -2.0/3.0);
    const scalar c = 0.5*kappa_*(sqr(J) - 1.0) - m*I1/3.0;

    // Left Cauchy-Green tensor B = I + a*g + g*a + G*a*a
    const symmTensor B = symmTensor(I) + twoSymm(a_*g) + G*sqr(a_);

    return (m*B + c*symmTensor(I))/J;
}


Foam::vector Foam::neoHookeanMMS::bodyForce
(
    const vector& X,
    const scalar t
) const
{
    scalar q = 0;
    vector g = vector::zero;
    symmTensor H = symmTensor::zero;
    scalar lapq = 0;
    scalar d2qdt2 = 0;
    evaluate(X, t, q, g, H, lapq, d2qdt2);

    const scalar J = 1.0 + (a_ & g);

    if (J < SMALL)
    {
        FatalErrorIn("Foam::neoHookeanMMS::bodyForce(...)")
            << "Non-positive Jacobian J = " << J << " at X = " << X
            << ", t = " << t << ": reduce the amplitude a"
            << abort(FatalError);
    }

    const scalar A = a_ & a_;
    const scalar G = g & g;
    const scalar I1 = 3.0 + 2.0*(a_ & g) + A*G;

    const scalar m = mu_*Foam::pow(J, -2.0/3.0);
    const scalar c = 0.5*kappa_*(sqr(J) - 1.0) - m*I1/3.0;
    const scalar e = c/J;

    // h = H & a = grad0(J), r = H & g
    const vector h = H & a_;
    const vector r = H & g;

    const vector gradI1 = 2.0*h + 2.0*A*r;
    const vector gradm = -(2.0*m/(3.0*J))*h;
    const vector gradc = kappa_*J*h - (I1/3.0)*gradm - (m/3.0)*gradI1;
    const vector gradd = gradm + gradc;
    const vector grade = gradc/J - (c/sqr(J))*h;

    // div0(P) with P = d*I + m*a*g - e*g*a
    const vector divP =
        gradd
      + a_*((g & gradm) + m*lapq)
      - g*(a_ & grade)
      - e*h;

    return rho_*a_*d2qdt2 - divP;
}


const Foam::volVectorField& Foam::neoHookeanMMS::bodyForces() const
{
    if
    (
        !bodyForcesPtr_.valid()
     || bodyForcesTimeIndex_ != mesh_.time().timeIndex()
    )
    {
        calcBodyForces();
    }

    return bodyForcesPtr_();
}

// ************************************************************************* //
