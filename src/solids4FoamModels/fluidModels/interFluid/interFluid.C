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

#include "interFluid.H"
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

defineTypeNameAndDebug(interFluid, 0);
addToRunTimeSelectionTable(physicsModel, interFluid, fluid);
addToRunTimeSelectionTable(fluidModel, interFluid, dictionary);


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //


void interFluid::correctPhi
(
    pimpleControl& pimple,
    const label pdRefCell,
    const scalar pdRefValue
)
{
    fluidModel::continuityErrs();

    // Note: we should store pcorrTypes
    wordList pcorrTypes
    (
        pd_.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i = 0; i < pd_.boundaryField().size(); i++)
    {
        if (pd_.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("pcorr", pd_.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf
    (
        "(1|A(U))",
        dimTime/rho_.dimensions(),
        runTime().deltaT().value()
    );

    phi() = (fvc::interpolate(U()) & mesh().Sf());

    adjustPhi(phi(), U(), pcorr);

    mesh().schemesDict().setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi())
        );

        pcorrEqn.setReference(pdRefCell, pdRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi() -= pcorrEqn.flux();
        }
    }

    fluidModel::continuityErrs();

    // Courant number
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;
        fluidModel::CourantNo(CoNum, meanCoNum, velMag);
    }

    // Recalculate rhoPhi from rho
    rhoPhi_ = fvc::interpolate(rho_)*phi();
}


void interFluid::solveAlphaEqnSubCycle(const pimpleControl& pimple)
{
    const label nAlphaCorr
    (
        readLabel(pimple.dict().lookup("nAlphaCorr"))
    );

    const label nAlphaSubCycles
    (
        readLabel(pimple.dict().lookup("nAlphaSubCycles"))
    );

    if (nAlphaSubCycles > 1)
    {
        const dimensionedScalar totalDeltaT = runTime().deltaT();
        surfaceScalarField rhoPhiSum = 0.0*rhoPhi_;

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1_, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphaEqn(nAlphaCorr);
            rhoPhiSum += (runTime().deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphaEqn(nAlphaCorr);
    }

    interface_.correct();

    rho_ == alpha1_*rho1_ + (scalar(1) - alpha1_)*rho2_;
}


void interFluid::solveAlphaEqn(const label nAlphaCorr)
{
    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi()/mesh().magSf());
    phic = min(interface_.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface_.nHatf();

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi(),
                alpha1_,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1_, alpharScheme),
                alpha1_,
                alpharScheme
            );

        MULES::explicitSolve(alpha1_, phi(), phiAlpha, 1, 0);

        rhoPhi_ = phiAlpha*(rho1_ - rho2_) + phi()*rho2_;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1_.weightedAverage(mesh().V()).value()
        << "  Min(alpha1) = " << min(alpha1_).value()
        << "  Max(alpha1) = " << max(alpha1_).value()
        << endl;
}


tmp<fvVectorMatrix> interFluid::solveUEqn(pimpleControl& pimple)
{
    const surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties_.muf()
      + fvc::interpolate(rho_*turbulence_->nut())
    );

    tmp<fvVectorMatrix> tUEqn
    (
        new fvVectorMatrix
        (
            fvm::ddt(rho_, U())
          + fvm::div(rhoPhi_, U())
          - fvm::laplacian(muEff, U())
          - (fvc::grad(U()) & fvc::grad(muEff))
        )
    );
    fvVectorMatrix& UEqn = tUEqn();

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_)
                  - ghf_*fvc::snGrad(rho_)
                  - fvc::snGrad(pd_)
                )*mesh().magSf()
            )
        );
    }

    return tUEqn;
}


void interFluid::solvePEqn
(
    pimpleControl& pimple,
    fvVectorMatrix& UEqn,
    const label pdRefCell,
    const scalar pdRefValue
)
{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U() = rUA*UEqn.H();

    surfaceScalarField phiU("phiU", (fvc::interpolate(U()) & mesh().Sf()));

    if (pd_.needReference())
    {
        adjustPhi(phi(), U(), pd_);
    }

    phi() = phiU +
    (
        fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_)
      - ghf_*fvc::snGrad(rho_)
    )*rUAf*mesh().magSf();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd_) == fvc::div(phi())
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.solve
        (
            mesh().solutionDict().solver(pd_.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi() -= pdEqn.flux();
        }
    }

    U() += rUA*fvc::reconstruct((phi() - phiU)/rUAf);
    U().correctBoundaryConditions();

    fluidModel::continuityErrs();

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interFluid::interFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    pd_
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    alpha1_
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    twoPhaseProperties_(U(), phi(), "alpha1"),
    rho1_(twoPhaseProperties_.rho1()),
    rho2_(twoPhaseProperties_.rho2()),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT
        ),
        alpha1_*rho1_ + (scalar(1) - alpha1_)*rho2_,
        alpha1_.boundaryField().types()
    ),
    rhoPhi_
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1_*phi()
    ),
    gh_("gh", g() & mesh().C()),
    ghf_("gh", g() & mesh().Cf()),
    interface_(alpha1_, U(), twoPhaseProperties_),
    pdRefCell_(0),
    pdRefValue_(0.0),
    pRefValue_(0.0),
    turbulence_
    (
        incompressible::turbulenceModel::New(U(), phi(), twoPhaseProperties_)
    )
{
    UisRequired();

    // Reset p dimensions
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure);
    p() = pd_ + rho_*gh_;

    rho_.oldTime();

    setRefCell(p(), fluidProperties(), pdRefCell_, pdRefValue_);
    mesh().schemesDict().setFluxRequired(pd_.name());

    if (pd_.needReference())
    {
        pRefValue_ = readScalar(pimple().dict().lookup("pRefValue"));

        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pdRefCell_)
        );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> interFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() =
        (
            mesh().boundary()[patchID].nf()
          & turbulence_->devReff()().boundaryField()[patchID]
        );

    vectorField n = mesh().boundary()[patchID].nf();
    tvF() -= (sqr(n) & tvF());

    return tvF;
}


tmp<scalarField> interFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = p().boundaryField()[patchID];

    return tpF;
}


tmp<scalarField> interFluid::faceZoneMuEff
(
    const label zoneID,
    const label patchID
) const
{
    scalarField pMuEff =
       turbulence_->nuEff()().boundaryField()[patchID];

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


bool interFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    fvMesh& mesh = fluidModel::mesh();

    // Take a reference to the pimple control
    pimpleControl& pimple = fluidModel::pimple();

    bool meshChanged = false;
    if (fluidModel::fsiMeshUpdate())
    {
        // The FSI interface is in charge of calling mesh.update()
        meshChanged = fluidModel::fsiMeshUpdateChanged();
    }
    else
    {
        meshChanged = refCast<dynamicFvMesh>(mesh).update();
        reduce(meshChanged, orOp<bool>());
    }

    {
        const Time& runTime = fluidModel::runTime();
#       include "volContinuity.H"
    }

    // Update gh fields as the mesh may have moved
    gh_ = g() & mesh.C();
    ghf_ = g() & mesh.Cf();

    //if (correctPhi && meshChanged)
    if (meshChanged)
    {
        correctPhi(pimple, pdRefCell_, pdRefValue_);
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi(), U());

    // if (checkMeshCourantNo)
    // {
    // #       include "meshCourantNo.H"
    // }

    // Pressure-velocity corrector
    while (pimple.loop())
    {
        twoPhaseProperties_.correct();

        solveAlphaEqnSubCycle(pimple);

        fvVectorMatrix UEqn = solveUEqn(pimple);

        // --- PISO loop
        while (pimple.correct())
        {
            solvePEqn(pimple, UEqn, pdRefCell_, pdRefValue_);
        }

        p() = pd_ + rho_*gh_;

        if (pd_.needReference())
        {
            p() +=
                dimensionedScalar
                (
                    "p",
                    p().dimensions(),
                    pRefValue_ - getRefCellValue(p(), pdRefCell_)
                );
        }

        gradp() = fvc::grad(p());

        gradU() = fvc::grad(U());

        turbulence_->correct();
    }

    // Make the fluxes absolute for when runTime++ is called
    fvc::makeAbsolute(phi(), U());

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
