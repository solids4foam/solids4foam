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

#include "pimpleFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "fvc.H"
#include "fvm.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pimpleFluid, 0);
addToRunTimeSelectionTable(fluidModel, pimpleFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pimpleFluid::pimpleFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    laminarTransport_(U(), phi()),
    turbulence_
    (
        incompressible::turbulenceModel::New
        (
            U(), phi(), laminarTransport_
        )
    ),
    rho_(laminarTransport_.lookup("rho")),
    ddtU_(fvc::ddt(U())),
    Uf_
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U())
    ),
    aU_
    (
        IOobject
        (
            "aU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        1.0/runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    ),
    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    ),
    rAUf_
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        runTime.deltaT(),
        calculatedFvPatchScalarField::typeName
    ),
    pRefCell_(0),
    pRefValue_(0),
    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    useFoamExtend40formulation_
    (
        fluidProperties().lookupOrDefault("useFoamExtend40formulation", true)
    ),
    robin_(U(), p()),
    sumLocalContErr_(0),
    globalContErr_(0),
    cumulativeContErr_(0)
{
    UisRequired();
    pisRequired();

    phi().oldTime();
    Uf_.oldTime().oldTime();

    ddtU_.checkIn();

    setRefCell(p(), pimple().dict(), pRefCell_, pRefValue_);
    mesh().schemesDict().setFluxRequired(p().name());

    if (useFoamExtend40formulation_ && robin_.patchIDs().size() > 0)
    {
        FatalErrorIn(type())
            << "Currently, the Robin coupling approach only works with "
            << "useFoamExtend40formulation turned off" << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> pimpleFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() = rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> pimpleFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


bool pimpleFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = this->mesh();
    pimpleControl& pimple = this->pimple();
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    surfaceScalarField& phi = this->phi();
    volScalarField& aU = aU_;
    volScalarField& rAU = rAU_;
    surfaceScalarField& rAUf = rAUf_;
    const label pRefCell = pRefCell_;
    const scalar& pRefValue = pRefValue_;
    scalar& sumLocalContErr = sumLocalContErr_;
    scalar& globalContErr = globalContErr_;
    scalar& cumulativeContErr = cumulativeContErr_;
    autoPtr<incompressible::turbulenceModel>& turbulence = turbulence_;
    const Switch correctPhi = correctPhi_;
    const Switch checkMeshCourantNo = checkMeshCourantNo_;
    const Switch useFoamExtend40formulation = useFoamExtend40formulation_;

    // Fluxes should be absolute at the end of the time-step
    // Make the fluxes absolute
    fvc::makeAbsolute(phi, U);

    bool meshChanged = false;
    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = mesh.update();
    }
    reduce(meshChanged, orOp<bool>());

#   include "volContinuity.H"

    if (correctPhi && meshChanged)
    {
        // Fluxes will be corrected to absolute velocity
        if (useFoamExtend40formulation)
        {
#           include "correctPhi.pimpleFoam.foamextend40.H"
        }
        else
        {
#           include "correctPhi.pimpleFoam.foamextend41.H"
        }
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    if (meshChanged)
    {
        // Update flux in FSI interface with Robin BC
        robin_.setInterfaceFluxToZero(phi);
    }

    // CourantNo
    if (mesh.moving() && checkMeshCourantNo)
    {
#       include "meshCourantNo.H"
    }

    if (meshChanged)
    {
#       include "CourantNo.H"
    }

    // --- PIMPLE loop
    while (pimple.loop())
    {
        if (useFoamExtend40formulation)
        {
#           include "UEqn.pimpleFoam.foamextend40.H"

            // --- PISO loop
            while (pimple.correct())
            {
#               include "pEqn.pimpleFoam.foamextend40.H"
            }
        }
        else
        {
#           include "UEqn.pimpleFoam.foamextend41.H"

            // --- PISO loop
            while (pimple.correct())
            {
#               include "pEqn.pimpleFoam.foamextend41.H"
            }
        }

        turbulence_->correct();
    }

    // Make the fluxes absolute to the mesh motion
    fvc::makeAbsolute(phi, U);

    // Update velocity on faces to account for modifications in the flux
    // when using the Robin BCs
    Uf_ = fvc::interpolate(U, "interpolate(U)");

    // Get tangential velocity
    Uf_ -= (mesh.Sf() & Uf_)*mesh.Sf()/magSqr(mesh.Sf());

    // Update with normal from flux
    Uf_ += phi*mesh.Sf()/magSqr(mesh.Sf());

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
