/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    K_(dict.lookup("K")),
    cf_(dict.lookup("cf")),
    ct_(dict.lookup("ct")),
    cfs_(dict.lookup("cfs")),
    mu_(max(max(cf_, ct_), cfs_)),
    A_((cf_ - 2.0*cfs_ + ct_)/2.0),
    B_((cfs_ - ct_)/4.0),
    C_(cf_/4.0),
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
    s0_
    (
        IOobject
        (
            "s0",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    )
{
    // Check that f0 and so are unit vectors

    if (max(mag(mag(f0_) - 1.0)).value() > SMALL)
    {
        FatalErrorIn("GuccioneElastic::GuccioneElastic()")
            << "f0 should be unit vectors"
            << abort(FatalError);
    }

    if (max(mag(mag(s0_) - 1.0)).value() > SMALL)
    {
        FatalErrorIn("GuccioneElastic::GuccioneElastic()")
            << "s0 should be unit vectors"
            << abort(FatalError);
    }

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
    WarningIn("uccioneElastic::impK()")
        << "Check if impK is set appropriately using an effective mu"
        << endl;

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
            (4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::GuccioneElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Take a reference to the deformation gradient to make the code easier to
    // read
    const volTensorField& F = this->F();

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F));

    // Inverse deformation gradient
    const volTensorField Finv(inv(F));

    // Calculate the co-factor tensor
    const volTensorField H(J*Finv.T());

    // Calculate invariants
    const volScalarField I4f0((F & f0_) & (F & f0_));
    const volScalarField I4s0((F & s0_) & (F & s0_));
    const volScalarField I8f0s0((F & f0_) & (F & s0_));
    const volScalarField I4Hf0((H & f0_) & (H & f0_));
    const volScalarField IIF(magSqr(F)/2.0);
    const volScalarField IIH(magSqr(H)/2.0);

    // Calculate six terms required for exponent Q
    const volScalarField QAniso1(sqr(I4f0 - 1.0));
    const volScalarField QAniso2(IIF*I4f0 - IIH + I4Hf0 - 2.0*I4f0 + 1.0);
    const volScalarField QIso(sqr(IIF) - 2.0*IIH - 2.0*IIF + 3.0);

    // Calculate exponent Q based on F, H and J
    const volScalarField Q(A_*QAniso1 + B_*QAniso2 + C_*QIso);

    // Calculate strain energy density
    const volScalarField U((K_/2.0)*(exp(Q) - 1.0));

    // Second Garcia-Blanco et al paper calculate the first Piola-Kirchoff
    // stress as
    // P = dU/dF + (d(U)/dH x F) + d(U*H)/dJ
    // I should check the two earlier papers by Gil and Ortigosa for details, or
    // The thesis of Garcia-Blanco
    // Also, the cross-product "x" for tensors is not yet defined in OpenFOAM

    // Convert P to S

    // Convert the second Piola-Kirchhoff stress to the Cauchy stress
    // sigma = J*symm(F & S & F.T());;
}


void Foam::GuccioneElastic::correct(surfaceSymmTensorField& sigma)
{
    notImplemented("GuccioneElastic::correct(surfaceSymmTensorField& sigma)");

    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }
}


void Foam::GuccioneElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
