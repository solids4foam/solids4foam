/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "pUCoupledFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pUCoupledFluid, 0);
addToRunTimeSelectionTable(physicsModel, pUCoupledFluid, fluid);
addToRunTimeSelectionTable(fluidModel, pUCoupledFluid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pUCoupledFluid::pUCoupledFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    laminarTransport_(U(), phi()),
    turbulence_
    (
        incompressible::RASModel::New(U(), phi(), laminarTransport_)
    ),
    Up_
    (
        IOobject
        (
           "Up",
           runTime.timeName(),
           mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector4("zero", dimless, vector4::zero)
    ),
    rAU_
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        runTime.deltaT()
    ),
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
    rho_(transportProperties_.lookup("rho")),
    maxResidual_(0),
    convergenceCriterion_(0)
{
    mesh().schemesDict().setFluxRequired(p().name());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<vectorField> pUCoupledFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        rho_.value()
       *(
            mesh().boundary()[patchID].nf()
          & turbulence_->devReff()().boundaryField()[patchID]
        );

    const vectorField n = mesh().boundary()[patchID].nf();
    tvF() -= (sqr(n) & tvF());

    return tvF;
}


tmp<scalarField> pUCoupledFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> pUCoupledFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
        rho_.value()*turbulence_->nuEff()().boundaryField()[patchID];

    tmp<scalarField> tMuEff
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& muEff = tMuEff();

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(pMuEff, I)
    {
        muEff[mesh().faceZones()[zoneID].whichFace(patchStart + I)] =
            pMuEff[I];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(muEff, sumOp<scalarField>());

    return tMuEff;
}


bool pUCoupledFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    const fvMesh& mesh = this->mesh();

    label pRefCell = 0;
    scalar pRefValue = 0;
    setRefCell
    (
        p(),
        mesh.solutionDict().subDict("blockSolver"),
        pRefCell,
        pRefValue
    );

    mesh.solutionDict().subDict("blockSolver").readIfPresent
    (
        "convergence",
        convergenceCriterion_
    );

    p().storePrevIter();

    // Initialize the Up block system (matrix, source and reference to Up)
    fvBlockMatrix<vector4> UpEqn(Up_);

    // Assemble and insert momentum equation
    {
        volScalarField divPhi
        (
            "divPhi",
            fvc::div(phi())
        );

        // Momentum equation
        fvVectorMatrix UEqn
        (
            fvm::div(phi(), U())
          + turbulence_->divDevReff()
        );

        rAU_ = 1.0/UEqn.A();

        // Insert the additional components. Note this will destroy the H and A

        UEqn += fvm::SuSp(-divPhi, U()) + divPhi*U();
        UEqn.relax();

        UpEqn.insertEquation(0, UEqn);
    }

    // Assemble and insert pressure equation

    // Pressure parts of the continuity equation
    surfaceScalarField presSource
    (
        "presSource",
        fvc::interpolate(rAU_)*
        (fvc::interpolate(fvc::grad(p())) & mesh.Sf())
    );

    fvScalarMatrix pEqn
    (
      - fvm::laplacian(rAU_, p())
     ==
      - fvc::div(presSource)
    );

    pEqn.setReference(pRefCell, pRefValue);

    UpEqn.insertEquation(3, pEqn);

    // Assemble and insert coupling terms
    {
        // Calculate grad p coupling matrix. Needs to be here if one uses
        // gradient schemes with limiters.  VV, 9/June/2014
        BlockLduSystem<vector, vector> pInU(fvm::grad(p()));

        // Calculate div U coupling.  Could be calculated only once since
        // it is only geometry dependent.  VV, 9/June/2014
        BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U()));

        // Last argument in insertBlockCoupling says if the column direction
        // should be incremented. This is needed for arbitrary positioning
        // of U and p in the system. This could be better. VV, 30/April/2014
        UpEqn.insertBlockCoupling(0, 3, pInU, true);
        UpEqn.insertBlockCoupling(3, 0, UInp, false);
    }

    // Solve the block matrix
    maxResidual_ = cmptMax(UpEqn.solve().initialResidual());

    // Retrieve solution
    UpEqn.retrieveSolution(0, U().internalField());
    UpEqn.retrieveSolution(3, p().internalField());

    U().correctBoundaryConditions();
    p().correctBoundaryConditions();

    phi() = (fvc::interpolate(U()) & mesh.Sf()) + pEqn.flux() + presSource;

    fluidModel::continuityErrs();

    fluidModel::boundPU(p(), U());

    p().relax();

    turbulence_->correct();

    gradU() = fvc::grad(U());

    // Check convergence
    if (maxResidual_ < convergenceCriterion_)
    {
        Info<< "reached convergence criterion: " << convergenceCriterion_ << endl;
        // For now, we use const_cast, but we can do this better
        //const_cast<Time&>(runTime()).writeAndEnd();
        Info<< "latestTime = " << runTime().timeName() << endl;
    }

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels

} // End namespace Foam

// ************************************************************************* //
