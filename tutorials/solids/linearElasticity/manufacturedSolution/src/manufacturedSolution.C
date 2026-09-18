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

#include "manufacturedSolution.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manufacturedSolution, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::manufacturedSolution::calcBodyForces() const
{
    if (bodyForcesPtr_.valid())
    {
        FatalErrorIn("void Foam::manufacturedSolution::calcBodyForces() const")
            << "Pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = mesh_;

    bodyForcesPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "bodyForces",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimForce/dimVolume, vector::zero),
            "zeroGradient"
        )
    );

    // Set the body force term
    const vectorField& C = mesh.C();
    const scalar pi = constant::mathematical::pi;
    vectorField& bodyForcesI = bodyForcesPtr_();
    forAll(bodyForcesI, cellI)
    {
        const scalar x = C[cellI].x();
        const scalar y = C[cellI].y();
        const scalar z = C[cellI].z();

        bodyForcesI[cellI][vector::X] =
            lambda_
           *(
                8*ay_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z)
              + 4*az_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y)
              - 16*ax_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)
            )
          + mu_
           *(
                8*ay_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z)
              + 4*az_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y)
              - 5*ax_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)
            )
          - 32*ax_*mu_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);

        bodyForcesI[cellI][vector::Y] =
            lambda_
           *(
                8*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z)
              + 2*az_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x)
              - 4*ay_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)
            )
          + mu_
           *(
                8*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(2*pi*y)*Foam::sin(pi*z)
              + 2*az_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x)
              - 17*ay_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)
           )
          - 8*ay_*mu_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);

        bodyForcesI[cellI][vector::Z] =
            lambda_
           *(
               4*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y)
              + 2*ay_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x)
              - az_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)
            )
          + mu_
           *(
               4*ax_*pi*pi*Foam::cos(4*pi*x)*Foam::cos(pi*z)*Foam::sin(2*pi*y)
              + 2*ay_*pi*pi*Foam::cos(2*pi*y)*Foam::cos(pi*z)*Foam::sin(4*pi*x)
              - 20*az_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z)
            )
          - 2*az_*mu_*pi*pi*Foam::sin(4*pi*x)*Foam::sin(2*pi*y)*Foam::sin(pi*z);
    }

    bodyForcesPtr_().correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manufacturedSolution::manufacturedSolution
(
    const fvMesh& mesh,
    const scalar ax,
    const scalar ay,
    const scalar az,
    const scalar E,
    const scalar nu
)
:
    mesh_(mesh),
    ax_(ax),
    ay_(ay),
    az_(az),
    E_(E),
    nu_(nu),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_((E_*nu_)/((1.0 + nu_)*(1.0 - 2.0*nu_)))
{
    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn("manufacturedSolution::manufacturedSolution(...)")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}


Foam::manufacturedSolution::manufacturedSolution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    ax_(readScalar(dict.lookup("ax"))),
    ay_(readScalar(dict.lookup("ay"))),
    az_(readScalar(dict.lookup("az"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_((E_*nu_)/((1.0 + nu_)*(1.0 - 2.0*nu_)))
{
    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn("manufacturedSolution::manufacturedSolution(...)")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::symmTensor Foam::manufacturedSolution::calculateStress
(
    const vector& point
)
{
    symmTensor sigma = symmTensor::zero;

    // Shear modulus
    const scalar mu = E_/(2.0*(1.0 + nu_));

    // Lambda parameter
    const scalar lambda = (E_*nu_)/((1.0 + nu_)*(1.0 - 2.0*nu_));

    // pi
    const scalar pi = Foam::constant::mathematical::pi;

    sigma.xx() =
        lambda*(4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
      + 2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
      + az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))
      + 8*ax_*mu*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    sigma.yy() =
        lambda*(4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
      + 2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
      + az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))
      + 4*ay_*mu*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    sigma.zz() =
        lambda*(4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
      + 2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
      + az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))
      + 2*az_*mu*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y());

    sigma.xy() =
        2*ax_*mu*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z())
     +  4*ay_*mu*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    // sigma.yx() = sigma.xy();

    sigma.yz() =
        ay_*mu*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())
      + 2*az_*mu*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    // sigma.zy() = sigma.yz();

    sigma.xz() =
        ax_*mu*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())
      + 4*az_*mu*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    // sigma.zx() = sigma.xz();

    return sigma;
}


Foam::symmTensor Foam::manufacturedSolution::calculateStrain
(
    const vector& point
)
{
    symmTensor epsilon = symmTensor::zero;

    //pi
    const scalar pi = Foam::constant::mathematical::pi;

    epsilon.xx() =
        4*ax_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    epsilon.yy() =
        2*ay_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    epsilon.zz() =
        az_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y());

    epsilon.xy() =
        ax_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z()) + 2*ay_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    // epsilon.yx() = epsilon.xy();

    epsilon.yz() =
        (ay_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))/2 + az_*pi*Foam::cos(2*pi*point.y())*Foam::sin(4*pi*point.x())*Foam::sin(pi*point.z());

    // epsilon.zy() = epsilon.yz();

    epsilon.xz() =
        (ax_*pi*Foam::cos(pi*point.z())*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y()))/2 + 2*az_*pi*Foam::cos(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z());

    // epsilon.zx() = epsilon.xz();

    return epsilon;
}

Foam::vector Foam::manufacturedSolution::calculateDisplacement
(
    const vector& point
)
{
    //pi
    const scalar pi = Foam::constant::mathematical::pi;

    return vector
    (
        ax_*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z()),
        ay_*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z()),
        az_*Foam::sin(4*pi*point.x())*Foam::sin(2*pi*point.y())*Foam::sin(pi*point.z())
    );
}


const Foam::volVectorField& Foam::manufacturedSolution::bodyForces() const
{
    if (!bodyForcesPtr_.valid())
    {
        calcBodyForces();
    }

    return bodyForcesPtr_();
}

// ************************************************************************* //
