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
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interFluid, 0);
addToRunTimeSelectionTable(fluidModel, interFluid, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interFluid::interFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    LTS_(fv::localEulerDdt::enabled(mesh())),
    trDeltaT_(),
    Uf_(),
    mixture_(U(), phi()),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT
        ),
        mixture_.alpha1()*mixture_.rho1()
      + (scalar(1) - mixture_.alpha1())*mixture_.rho2(),
        mixture_.alpha1().boundaryField().types()
    ),
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho_)*phi()
    ),
    turbulence_
    (
        new transportModelType
        (
            rho_, U(), phi(), rhoPhi_, mixture_
        )
    ),
    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    hRef_
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, 0)
    ),
    ghRef_(- mag(g())*hRef_),
    gh_("gh", (g() & mesh().C()) - ghRef_),
    ghf_("ghf", (g() & mesh().Cf()) - ghRef_),
    pRefCell_(0.0),
    pRefValue_(0.0),
    rAU_(),
    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    moveMeshOuterCorrectors_
    (
        pimple().dict().lookupOrDefault("moveMeshOuterCorrectors", false)
    ),
    cumulativeContErr_(0),
    alphaPhiUn_
    (
        IOobject
        (
            "alphaPhiUn",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(phi().dimensions(), Zero)
    ),
    alphaPhi10Header_
    (
        IOobject::groupName("alphaPhi0", mixture_.alpha1().group()),
        runTime.timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaRestart_(alphaPhi10Header_.typeHeaderOk<surfaceScalarField>()),
    alphaPhi10_
    (
        alphaPhi10Header_,
        phi()*fvc::interpolate(mixture_.alpha1())
    ),
    talphaPhi1Corr0_()
{
    UisRequired();

    // Reset p dimensions: we should allow p not to be read!
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure);
    p() = p_rgh_ + rho_*gh_;

    // Set pressure reference
    if (p_rgh_.needReference())
    {
        p() += dimensionedScalar
        (
            "p",
            p().dimensions(),
            pRefValue_ - getRefCellValue(p(), pRefCell_)
        );
        p_rgh_ = p() - rho_*gh_;
    }

    mesh().setFluxRequired(p_rgh_.name());
    mesh().setFluxRequired(mixture_.alpha1().name());

    // Store old-time rho
    rho_.oldTime();

    {
        volVectorField& U = this->U();
        volScalarField& p_rgh = p_rgh_;
        const dynamicFvMesh& mesh = this->mesh();
        pimpleControl& pimple = this->pimple();
        surfaceScalarField& phi = this->phi();    
        const bool correctPhi = correctPhi_;
        const volScalarField& rho = rho_;
        scalar& cumulativeContErr = cumulativeContErr_;

        #include "initCorrectPhi.esi.H"
    }

    if (alphaRestart_)
    {
        Info<< "Restarting alpha" << endl;
    }

    if (LTS_)
    {
        Info<< "Using LTS" << endl;

        trDeltaT_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
    else
    {
        const fvMesh& mesh = this->mesh();
        const surfaceScalarField& phi = this->phi();
        #include "CourantNo.H"
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> interFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    // Is this the easiest way to do this? Surely it should be possible from
    // the turbulence_ object...
    typedef incompressible::turbulenceModel icoTurbModel;
    const incompressible::turbulenceModel& turb =
        mesh().lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

    tvF.ref() =
        (
            mesh().boundary()[patchID].nf()
          & (
               - rho_.boundaryField()[patchID]
                *turb.devReff()().boundaryField()[patchID]
            )
          // & (-turbulence_->devReff()().boundaryField()[patchID])
        );

    return tvF;
}


tmp<scalarField> interFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(p().boundaryField()[patchID])
    );

    return tpF;
}


bool interFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = fluidModel::mesh();
    pimpleControl& pimple = fluidModel::pimple();
    volVectorField& U = this->U();
    volScalarField& p = this->p();
    volScalarField& p_rgh = p_rgh_;
    surfaceScalarField& phi = this->phi();    
    autoPtr<transportModelType>& turbulence = turbulence_;
    immiscibleIncompressibleTwoPhaseMixture& mixture = mixture_;
    volScalarField& alpha1 = mixture.alpha1();
    volScalarField& alpha2 = mixture.alpha2();
    volScalarField& rho = rho_;
    surfaceScalarField& rhoPhi = rhoPhi_;
    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();
    autoPtr<surfaceVectorField>& Uf = Uf_;
    scalar& cumulativeContErr = cumulativeContErr_;
    const label pRefCell = pRefCell_;
    const scalar pRefValue = pRefValue_;
    const bool correctPhi = correctPhi_;
    const bool checkMeshCourantNo = checkMeshCourantNo_;
    const bool moveMeshOuterCorrectors = moveMeshOuterCorrectors_;
    const bool LTS = LTS_;
    const bool alphaRestart = alphaRestart_;
    tmp<surfaceScalarField>& talphaPhi1Corr0 = talphaPhi1Corr0_;
    const uniformDimensionedVectorField g = this->g();
    dimensionedScalar& ghRef = ghRef_;
    volScalarField& gh = gh_;
    surfaceScalarField& ghf = ghf_;
    surfaceScalarField& alphaPhi10 = alphaPhi10_;
    tmp<volScalarField>& rAU = rAU_;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            if (fluidModel::fsiMeshUpdate())
            {
                // The FSI interface is in charge of calling mesh.update()
                fluidModel::fsiMeshUpdateChanged();
            }
            else
            {
                mesh.update();
            }

            if (mesh.changing())
            {
                // Do not apply previous time-step mesh compression flux
                // if the mesh topology changed
                if (mesh.topoChanging())
                {
                    talphaPhi1Corr0.clear();
                }

                gh = (g & mesh.C()) - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;

                // MRF not implemented
                // MRF.update();

                if (correctPhi)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh.Sf() & Uf();

                    #include "correctPhi.esi.H"

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }

        #include "alphaControls.esi.H"
        #include "alphaEqnSubCycle.esi.H"

        mixture.correct();

        if (pimple.frozenFlow())
        {
            continue;
        }

        #include "UEqn.esi.H"

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.esi.H"
        }

        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

// ************************************************************************* //
