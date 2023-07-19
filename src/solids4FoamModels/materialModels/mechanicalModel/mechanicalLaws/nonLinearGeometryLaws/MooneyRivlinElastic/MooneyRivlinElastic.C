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

#include "MooneyRivlinElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MooneyRivlinElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, MooneyRivlinElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::MooneyRivlinElastic::MooneyRivlinElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    c10_
    (
        IOobject
        (
            "c10",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("c10"))
    ),
    c01_
    (
        IOobject
        (
            "c01",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("c01"))
    ),
    c11_
    (
        IOobject
        (
            "c11",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("c11"))
    ),
    c10f_(fvc::interpolate(c10_)),
    c01f_(fvc::interpolate(c01_)),
    c11f_(fvc::interpolate(c11_))
{
    // Compute material properties
    // Shear modulus: based on pure shear stress state for hyperelasticity
    mu() = 2.0*(c10_ + c01_);

    // Bulk modulus
    // The user can specify K directly or the Poisson's ratio, nu
    if (dict.found("K"))
    {
        K() = dimensionedScalar(dict.lookup("K"));
    }
    else if (dict.found("nu") && !dict.found("K"))
    {
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));
        
        // Young's modulus
        const volScalarField E(6.0*(c10_ + c01_));

        // Compute K based on linear elasticity
        K() = E/(3.0*(1.0 - 2.0*nu));
    }
    else
    {
        FatalErrorIn
        (
            "MooneyRivlinElastic::MooneyRivlinElastic::()"
        )   << "Either K or nu elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    Info<< "Material properties "                 << endl
        << "    max(c10) = " << gMax(mag(c10_)()) << endl
        << "    max(c01) = " << gMax(mag(c01_)()) << endl
        << "    max(c11) = " << gMax(mag(c11_)()) << endl
        << "    max(mu) = "  << gMax(mag(mu())())  << endl
        << "    max(K) = "   << gMax(mag(K())())   << endl;

    // Update surface fields
    muf() = fvc::interpolate(mu());
    Kf()  = fvc::interpolate(K());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::MooneyRivlinElastic::impK() const
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
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            ((4.0/3.0)*mu() + K())
        )
    );
}


void Foam::MooneyRivlinElastic::correct
(
    volSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    if (updateF(sigma, mu(), K()))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Calculate the left Cauchy Green deformation tensor
    const volSymmTensorField isoB(pow(J, -2.0/3.0)*symm(F() & F().T()));
    const volSymmTensorField sqrB(symm(isoB & isoB));

    // Compute invariants fields
    const volScalarField I1(tr(isoB));
    const volScalarField I2(0.5*(pow(I1, 2.0) - tr(sqrB)));

    // Compute deviatoric stress
    const volSymmTensorField s
    (
        2.0*
        (
            c10_
          + c11_*(I2 - 3.0)
        )*isoB
      - 2.0*
        (
            c01_
          + c11_*(I1 - 3.0)
        )*inv(isoB)
    );

    // Update the hydrostatic stress
    updateSigmaHyd
    (
        0.5*K()*(pow(J, 2.0) - 1.0),
        (4.0/3.0)*mu() + K()
    );

    // Calculate the Cauchy stress
    sigma =
        (1.0/J)
       *(
           dev(s) + sigmaHyd()*I + symm(F() & sigma0() & F().T())
        );
}


void Foam::MooneyRivlinElastic::correct
(
    surfaceSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    if (updateF(sigma, muf(), Kf()))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Calculate the left Cauchy Green deformation tensor
    const surfaceSymmTensorField isoB(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));
    const surfaceSymmTensorField sqrB(symm(isoB & isoB));

    // Compute invariants fields
    const surfaceScalarField I1(tr(isoB));
    const surfaceScalarField I2(0.5*(pow(I1, 2.0) - tr(sqrB)));

    // Compute deviatoric stress
    const surfaceSymmTensorField s
    (
        2.0*
        (
            c10f_
          + c11f_*(I2 - 3.0)
        )*isoB
      - 2.0*
        (
            c01f_
          + c11f_*(I1 - 3.0)
        )*inv(isoB)
    );

    // Note: updateSigmaHyd is not used
    const surfaceScalarField sigmaHydf(0.5*Kf()*(pow(J, 2.0) - 1.0));

    // Calculate the Cauchy stress
    sigma =
        (1.0/J)
       *(
           dev(s) + sigmaHydf*I
         + symm(Ff() & linearInterpolate(sigma0()) & Ff().T())
        );
}


// ************************************************************************* //
