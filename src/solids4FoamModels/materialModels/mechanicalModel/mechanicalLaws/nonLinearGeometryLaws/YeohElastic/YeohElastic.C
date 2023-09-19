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

#include "YeohElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YeohElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, YeohElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::YeohElastic::YeohElastic
(
    const word& name,
    const fvMesh& mesh,
    dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    c1_
    (
        IOobject
        (
            "c1",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("c1"))
    ),
    c2_
    (
        IOobject
        (
            "c2",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("c2"))
    ),
    c3_
    (
        IOobject
        (
            "c3",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("c3"))
    ),
    c1f_(fvc::interpolate(c1_)),
    c2f_(fvc::interpolate(c2_)),
    c3f_(fvc::interpolate(c3_))
{
    // Compute material properties
    // Shear modulus: based on pure shear stress state for hyperelasticity
    mu() = 2.0*c1_;

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
        const volScalarField E(6.0*c1_);

        // Compute K based on linear elasticity
        K() = E/(3.0*(1.0 - 2.0*nu));
    }
    else
    {
        FatalErrorIn
        (
            "YeohElastic::YeohElastic::()"
        )   << "Either K or nu elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    Info<< "Material properties "                << endl
        << "    max(c1) = " << gMax(mag(c1_)())  << endl
        << "    max(c2) = " << gMax(mag(c2_)())  << endl
        << "    max(c3) = " << gMax(mag(c3_)())  << endl
        << "    max(mu) = "  << gMax(mag(mu())()) << endl
        << "    max(K) = "   << gMax(mag(K())())  << endl;

    // Update surface fields
    muf() = fvc::interpolate(mu());
    Kf()  = fvc::interpolate(K());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::YeohElastic::impK() const
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


void Foam::YeohElastic::correct(volSymmTensorField& sigma)
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

    // Compute invariants fields
    const volScalarField I1(tr(isoB));

    // Compute deviatoric stress
    const volSymmTensorField s
    (
        2.0*
        (
            c1_
          + 2.0*c2_*(I1 - 3.0)
          + 3.0*c3_*pow((I1 - 3.0), 2)
        )*isoB
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


void Foam::YeohElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    if (updateF(sigma, muf(), Kf()))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    const surfaceSymmTensorField isoB(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));

    // Compute invariants fields
    const surfaceScalarField I1(tr(isoB));

    // Compute deviatoric stress
    const surfaceSymmTensorField s
    (
        2.0*
        (
            c1f_
          + 2.0*c2f_*(I1 - 3.0)
          + 3.0*c3f_*pow((I1 - 3.0), 2)
        )*isoB
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
