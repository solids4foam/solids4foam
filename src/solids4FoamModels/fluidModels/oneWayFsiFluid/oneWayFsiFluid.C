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

#include "oneWayFsiFluid.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneWayFsiFluid, 0);
addToRunTimeSelectionTable(fluidModel, oneWayFsiFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oneWayFsiFluid::oneWayFsiFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nu_(transportProperties_.lookup("nu")),
    rho_(transportProperties_.lookup("rho"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> oneWayFsiFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tvF.ref() =
#else
    tvF() =
#endif
        rho_.value()*nu_.value()*U().boundaryField()[patchID].snGrad();

    return tvF;
}


tmp<scalarField> oneWayFsiFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tpF.ref() =
#else
    tpF() =
#endif
        rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


bool oneWayFsiFluid::evolve()
{
    Info << "Evolving fluid model" << endl;

    // Read the latest fluid mesh
    const_cast<dynamicFvMesh&>(mesh()).readUpdate();

    // Read the latest velocity and pressure fields

    IOobject Uheader
    (
        "U",
        runTime().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

    IOobject pheader
    (
        "p",
        runTime().timeName(),
        mesh(),
        IOobject::MUST_READ
    );

#ifdef OPENFOAM_NOT_EXTEND
    if
    (
        Uheader.typeHeaderOk<volVectorField>(true)
     && pheader.typeHeaderOk<volScalarField>(true)
    )
#else
    if (Uheader.headerOk() && pheader.headerOk())
#endif
    {
        U() = volVectorField(Uheader, mesh());

        p() = volScalarField(pheader, mesh());
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
