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

#include "GuccioneElastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GuccioneElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GuccioneElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::GuccioneElastic::GuccioneElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom,
    const label lawI
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom, lawI),
    bulkModulus_(dict.lookup("bulkModulus")),
    k_(dict.lookup("k")),
    cf_(readScalar(dict.lookup("cf"))),
    ct_(readScalar(dict.lookup("ct"))),
    cfs_(readScalar(dict.lookup("cfs"))),
    mu_(0.75*(cf_ - 2.0*cfs_ + 2.0*cfs_)*k_), // check: is this the equivalent of linear shear modulus?
    f0_
    (
        IOobject
        (
            "f0",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    f0f_
    (
        IOobject
        (
            "f0f",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    f0f0_("f0f0", sqr(f0_)),
    f0f0f_("f0f0f", sqr(f0f_))
{
    // Check that f0 and s0 are unit vectors

    if (min(mag(mag(f0_.primitiveField()))) < SMALL)
    {
        FatalErrorIn("GuccioneElastic::GuccioneElastic()")
            << "At least one f0 vector has a length of zero!"
            << abort(FatalError);
    }

    if (min(mag(mag(f0f_.primitiveField()))) < SMALL)
    {
        FatalErrorIn("GuccioneElastic::GuccioneElastic()")
            << "At least one f0f vector has a length of zero!"
            << abort(FatalError);
    }

    // Normalise f0
    f0_ /= mag(f0_);
    f0f_ /= mag(f0f_);

    // Re-calculate f0f0
    f0f0_ = sqr(f0_);
    f0f0f_ = sqr(f0f_);

    // Store old F
    F().storeOldTime();
    Ff().storeOldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GuccioneElastic::~GuccioneElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GuccioneElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (cf_ - 2.0*cfs_ + 2.0*cfs_)*k_ + bulkModulus_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::GuccioneElastic::bulkModulus() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            bulkModulus_
        )
    );
}


void Foam::GuccioneElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, bulkModulus_))
    {
        return;
    }

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // Calculate the right Cauchy–Green deformation tensor
    const volSymmTensorField C(symm(F.T() & F));

    // Calculate the Green-Lagrange strain
    const volSymmTensorField E(0.5*(C - I));

    // Calculate E . E
    const volSymmTensorField sqrE(symm(E & E));

    // Calculate the invariants of E
    const volScalarField I1(tr(E));
    const volScalarField I2(0.5*(sqr(tr(E)) - tr(sqrE)));
    const volScalarField I4(E && f0f0_);
    const volScalarField I5(sqrE && f0f0_);

    // Calculate Q
    const volScalarField Q
    (
        ct_*sqr(I1)
      - 2.0*ct_*I2
     + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
     + 2.0*(cfs_ - ct_)*I5
    );

    // Calculate the derivative of Q wrt to E
    const volSymmTensorField dQdE
    (
        2.0*ct_*E
      + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0_
      + 2.0*(cfs_ - ct_)*symm((E & f0f0_) + (f0f0_ & E))
    );

    // Calculate the 2nd Piola-Kirchhoff stress
    const volSymmTensorField S(dQdE*0.5*k_*exp(Q));

    // Update the hydrostatic stress
    updateSigmaHyd
    (
        0.5*bulkModulus_*(pow(J, 2.0) - 1.0)/J,
        (4.0/3.0)*k_ + bulkModulus_
    );

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and add
    // a hydrostatic pressure term
    sigma = J*symm(F & S & F.T()) + sigmaHyd()*I;
}


void Foam::GuccioneElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, bulkModulus_))
    {
        return;
    }

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const surfaceTensorField& F = this->Ff();

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(F));

    // Calculate the right Cauchy–Green deformation tensor
    const surfaceSymmTensorField C(symm(F.T() & F));

    // Calculate the Green-Lagrange strain
    const surfaceSymmTensorField E(0.5*(C - I));

    // Calculate E . E
    const surfaceSymmTensorField sqrE(symm(E & E));

    // Calculate the invariants of E
    const surfaceScalarField I1(tr(E));
    const surfaceScalarField I2(0.5*(sqr(tr(E)) - tr(sqrE)));
    const surfaceScalarField I4(E && f0f0f_);
    const surfaceScalarField I5(sqrE && f0f0f_);

    // Calculate Q
    const surfaceScalarField Q
    (
        ct_*sqr(I1)
      - 2.0*ct_*I2
     + (cf_ - 2.0*cfs_ + ct_)*sqr(I4)
     + 2.0*(cfs_ - ct_)*I5
    );

    // Calculate the derivative of Q wrt to E
    const surfaceSymmTensorField dQdE
    (
        2.0*ct_*E
      + 2.0*(cf_ - 2.0*cfs_ + ct_)*I4*f0f0f_
      + 2.0*(cfs_ - ct_)*symm((E & f0f0f_) + (f0f0f_ & E))
    );

    // Calculate the 2nd Piola-Kirchhoff stress with penalty pressure term
    const surfaceSymmTensorField S(dQdE*0.5*k_*exp(Q));

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress and add
    // a hydrostatic pressure term
    sigma = J*symm(F & S & F.T()) + 0.5*bulkModulus_*((pow(J, 2) - 1)/J)*I;
}


void Foam::GuccioneElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
