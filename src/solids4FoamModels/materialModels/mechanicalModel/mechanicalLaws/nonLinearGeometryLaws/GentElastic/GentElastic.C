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

#include "GentElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GentElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, GentElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::GentElastic::GentElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom,
    const label lawI,
    const solidSubMeshes* solidSubMeshes
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom, lawI, solidSubMeshes),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    temperature_(dict.lookup("temperature")),
    ilim_(dict.lookup("ilim")),
    K_
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    Na_
    (
        "Na", dimensionSet(0, 0, 0, 0, 0, 0, 0), 6.022141e23
    ),
    kb_
    (
        "kb", dimensionSet(1, 2, -2, -1, 0, 0, 0), 1.38064852e-23
    ),
    omega_(dimensionedScalar(dict.lookup("Vs"))/Na_),
    N_(dimensionedScalar(dict.lookup("N"))/omega_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GentElastic::~GentElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GentElastic::impK() const
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
            mesh(),
            ((4.0/3.0)*mu_ + K_)
        )
    );
}


void Foam::GentElastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Right tensor product
    const volSymmTensorField B(symm(F() & F().T()));

    // Left tensor product
    const volSymmTensorField C(symm(F().T() & F()));

    // Calculate the sigma Terzaghi (Gent)
    sigma = kb_*temperature_*(N_/J)*((ilim_/(ilim_ - tr(B) + 3))*C - I);
}


void Foam::GentElastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Jacobian of the deformation gradient
    const surfaceScalarField Jf(det(Ff()));

    // Right tensor product
    const surfaceSymmTensorField Bf(symm(Ff() & Ff().T()));

    // Left tensor product
    const surfaceSymmTensorField Cf(symm(Ff().T() & Ff()));

    // Calculate the sigma Terzaghi (Gent)
    sigma = kb_*temperature_*(N_/Jf)*((ilim_/(ilim_ - tr(Bf) + 3))*Cf - I);
}


void Foam::GentElastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    Ff().writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
