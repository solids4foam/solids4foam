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

#include "isotropicFungElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isotropicFungElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, isotropicFungElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::isotropicFungElastic::isotropicFungElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom,
    const label lawI
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom, lawI),
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
    c1f_(fvc::interpolate(c1_)),
    c2_(dict.lookup("c2"))
{
    // Compute material properties
    // Shear modulus: based on pure shear stress state for hyperelasticity
    mu() = c1_;

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
        const volScalarField E(3.0*c1_);

        // Compute K based on linear elasticity
        K() = E/(3.0*(1.0 - 2.0*nu));
    }
    else
    {
        FatalErrorIn
        (
            "isotropicFungElastic::isotropicFungElastic::()"
        )   << "Either K or nu elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    Info<< "Material properties "               << nl
        << "    max(c1) = " << gMax(mag(c1_)()) << nl
        << "    c2 = "      << c2_.value()      << nl
        << "    max(mu) = " << gMax(mag(mu())()) << nl
        << "    max(K) = "  << gMax(mag(K())())  << endl;

    // Set surface fields
    muf() = fvc::interpolate(mu());
    Kf()  = fvc::interpolate(K());

    // Force the initial stress field to be created/read
    sigma0();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::isotropicFungElastic::impK() const
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


void Foam::isotropicFungElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    if (updateF(sigma, mu(), K()))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Calculate the left Cauchy-Green deformation tensor
    const volSymmTensorField isoB(pow(J, -2.0/3.0)*symm(F() & F().T()));

    // Compute invariants field
    const volScalarField I1(tr(isoB));

    // Compute first derivative of strain-energy function
    const volScalarField psi1(c1_*exp(0.5*c2_*(I1 - 3)));

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
           psi1*dev(isoB) + sigmaHyd()*I + symm(F() & sigma0() & F().T())
        );
}


void Foam::isotropicFungElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    if (updateF(sigma, muf(), Kf()))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const surfaceScalarField J(det(Ff()));

    // Calculate the left Cauchy-Green deformation tensor
    const surfaceSymmTensorField isoB(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));

    // Compute invariants fields
    const surfaceScalarField I1(tr(isoB));

    // Compute first derivative of strain-energy function
    const surfaceScalarField psi1(c1f_*exp(0.5*c2_*(I1 - 3.0)));

    // Note: updateSigmaHyd is not used
    const surfaceScalarField sigmaHydf(0.5*Kf()*(pow(J, 2.0) - 1.0));

    // Calculate the Cauchy stress
    sigma =
        (1.0/J)
       *(
           psi1*dev(isoB) + sigmaHydf*I
         + symm(Ff() & linearInterpolate(sigma0()) & Ff().T())
       );
}

// ************************************************************************* //
