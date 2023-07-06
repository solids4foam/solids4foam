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

#ifdef OPENFOAMESI

#include "compressibleInterFluid.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(compressibleInterFluid, 0);
addToRunTimeSelectionTable(fluidModel, compressibleInterFluid, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

compressibleInterFluid::compressibleInterFluid
(
    Time& runTime,
    const word& region
)
:
    fluidModel(typeName, runTime, region),
    cumulativeContErr_(0),
    adjustTimeStep_
    (
        runTime.controlDict().getOrDefault("adjustTimeStep", false)
    ),
    maxCo_(runTime.controlDict().getOrDefault<scalar>("maxCo", 1)),
    maxDeltaT_(runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT)),
    correctPhi_(pimple().dict().lookupOrDefault("correctPhi", false)),
    checkMeshCourantNo_
    (
        pimple().dict().lookupOrDefault("checkMeshCourantNo", false)
    ),
    moveMeshOuterCorrectors_
    (
        pimple().dict().lookupOrDefault("moveMeshOuterCorrectors", false)
    ),
    LTS_(fv::localEulerDdt::enabled(mesh())),
    trDeltaT_(),
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
    mixture_(U(), phi()),
    alpha1_(mixture_.alpha1()),
    alpha2_(mixture_.alpha2()),
    rho1_(mixture_.thermo1().rho()),
    rho2_(mixture_.thermo1().rho()),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT
        ),
        alpha1_*rho1_ + alpha2_*rho2_
    ),
    pMin_("pMin", dimPressure, mixture_),
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
    dgdt_(alpha1_*fvc::div(phi())),
    alphaPhi10Header_
    (
        IOobject::groupName("alphaPhi0", alpha1_.group()),
        runTime.timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaRestart_(alphaPhi10Header_.typeHeaderOk<surfaceScalarField>(true)),
    alphaPhi10_
    (
        alphaPhi10Header_,
        phi()*fvc::interpolate(alpha1_)
    ),
    talphaPhi1Corr0_(),
    turbulence_
    (
        rho_, U(), phi(), rhoPhi_, alphaPhi10_, mixture_
    ),
    K_("K", 0.5*magSqr(U())),
    MRF_(mesh()),
    Uf_
    (
        new surfaceVectorField
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
        )
    )
{
    UisRequired();

    // initContinuityErrs.H
    {
        uniformDimensionedScalarField cumulativeContErrIO
        (
            IOobject
            (
                "cumulativeContErr",
                runTime.timeName(),
                "uniform",
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(dimless, Zero)
        );
        cumulativeContErr_ = cumulativeContErrIO.value();
    }

    // Reset p dimensions: we should allow p not to be read!
    Info<< "Resetting the dimensions of p" << endl;
    p().dimensions().reset(dimPressure);
    p() = p_rgh_ + rho_*gh_;

    mesh().setFluxRequired(p_rgh_.name());
    mesh().setFluxRequired(alpha1_.name());

    // Store old-time rho
    rho_.oldTime();

    if (alphaRestart_)
    {
        Info<< "Restarting alpha" << endl;
    }

    // CourantNo.H and setInitialDeltaT.H
    {
        const fvMesh& mesh = this->mesh();
        const surfaceScalarField& phi = this->phi();
        #include "CourantNo.H"

        const bool adjustTimeStep = adjustTimeStep_;
        const scalar maxCo = maxCo_;
        const scalar maxDeltaT = maxDeltaT_;
        #include "setInitialDeltaT.H"
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> compressibleInterFluid::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    // Is this the easiest way to do this? Surely it should be possible from
    // the turbulence_ object...
    typedef compressible::turbulenceModel turbModel;
    const compressible::turbulenceModel& turb =
        mesh().lookupObject<turbModel>(turbModel::propertiesName);

    tvF.ref() =
        (
            mesh().boundary()[patchID].nf()
          & (
               - turb.devRhoReff()().boundaryField()[patchID]
            )
        );

    return tvF;
}


tmp<scalarField> compressibleInterFluid::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mixture_.p().boundaryField()[patchID])
    );

    return tpF;
}


void compressibleInterFluid::setDeltaT(Time& runTime)
{
    if (!LTS_)
    {
        const fvMesh& mesh = this->mesh();
        const surfaceScalarField& phi = this->phi();
        const bool adjustTimeStep = adjustTimeStep_;
        const scalar maxCo = maxCo_;
        const scalar maxDeltaT = maxDeltaT_;
        twoPhaseMixtureThermo& mixture = mixture_;

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
    }
}


bool compressibleInterFluid::evolve()
{
    Info<< "Evolving fluid model: " << this->type() << endl;

    // Take references
    const Time& runTime = fluidModel::runTime();
    dynamicFvMesh& mesh = fluidModel::mesh();
    pimpleControl& pimple = fluidModel::pimple();
    volVectorField& U = this->U();
    //volScalarField& p = this->p(); // see below
    volScalarField& p_rgh = p_rgh_;
    surfaceScalarField& phi = this->phi();
    compressibleInterPhaseTransportModel& turbulence = turbulence_;
    twoPhaseMixtureThermo& mixture = mixture_;
    volScalarField& alpha1 = alpha1_;
    volScalarField& alpha2 = alpha2_;
    volScalarField& rho = rho_;
    volScalarField& rho1 = rho1_;
    volScalarField& rho2 = rho2_;
    volScalarField& K = K_;
    surfaceScalarField& rhoPhi = rhoPhi_;
    surfaceVectorField& Uf = Uf_.ref();
    scalar& cumulativeContErr = cumulativeContErr_;
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
    dimensionedScalar& pMin = pMin_;
    volScalarField& dgdt = dgdt_;
    tmp<volScalarField>& trDeltaT = trDeltaT_;
    fv::options& fvOptions = options();
    IOMRFZoneList& MRF = MRF_;

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    // Store divU and divUp from the previous mesh so that it can be mapped
    // and used in correctPhi to ensure the corrected phi has the
    // same divergence
    volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));

    if (LTS)
    {
        #include "setRDeltaT.H"
    }
    // else // this happens in the setDeltaT function
    // {
    //     #include "CourantNo.H"
    //     #include "alphaCourantNo.H"
    //     #include "setDeltaT.H"
    // }

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        if (pimple.firstIter() || moveMeshOuterCorrectors)
        {
            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            mesh.update();

            if (mesh.changing())
            {
                MRF.update();

                Info<< "Execution time for mesh.update() = "
                    << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                    << " s" << endl;

                gh = (g & mesh.C()) - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;
            }

            if ((mesh.changing() && correctPhi))
            {
                // Calculate absolute flux from the mapped surface velocity
                phi = mesh.Sf() & Uf;

                #include "correctPhi.H"

                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U);

                mixture.correct();
            }

            if (mesh.changing() && checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        #include "alphaControls.H"
        #include "compressibleAlphaEqnSubCycle.H"

        turbulence.correctPhasePhi();

        #include "UEqn.H"
        #include "TEqn.H"

        // --- Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }

        if (pimple.turbCorr())
        {
            turbulence.correct();
        }
    }

    rho = alpha1*rho1 + alpha2*rho2;

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    // Set fluidModel::p() == p;
    fluidModel::p() = 1.0*p;

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidModels
} // End namespace Foam

#endif

// ************************************************************************* //
