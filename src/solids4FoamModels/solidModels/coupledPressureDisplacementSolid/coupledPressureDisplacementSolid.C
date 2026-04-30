/*---------------------------------------------------------------------------* \
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

#include "coupledPressureDisplacementSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "wedgeFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "blockSymmPlaneFvPatchVectorField.H"
// #include "blockRadialSlipFvPatchVectorField.H"
#include "tractionPressureDisplacementFvPatchVectorField.H"
// #include "tractionPressureFvPatchScalarField.H"

#include "skewCorrectionVectors.H"
#include "EulerD2dt2Scheme.H"
#include "steadyStateD2dt2Scheme.H"

// #include "trapezoidalDdtScheme.H"
// #include "trapezoidalD2dt2Scheme.H"

#include "diagVectorField.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "steadyStateDdtScheme.H"
#include "extrapolatedFvPatchFields.H"
#include "fvcGradf.H"
#include "fvcInterpolate.H"

#include "stabLeastSquaresGrad.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"

#include "correctedSnGrad.H"
#include "simpleUncorrectedSnGrad.H"
#include "gaussLaplacianScheme.H"
#include "harmonic.H"

#include "ggiFvPatchVectorNFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledPressureDisplacementSolid, 0);
// addToRunTimeSelectionTable
// (
//     physicsModel, coupledPressureDisplacementSolid, solid
// );
addToRunTimeSelectionTable
(
    solidModel, coupledPressureDisplacementSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


scalar coupledPressureDisplacementSolid::residual(const volVectorField& DD) const
{
    scalar refRes =
        gMax
        (
            mag
            (
                DD.internalField()
                // D.internalField() - D.oldTime().internalField()
            )
        );

    if (refRes < SMALL)
    {
        refRes = 10*refRes/solutionTol();
    }

    return
        gMax
        (
            mag(DD.internalField() - DD.prevIter().internalField())/refRes
            // /max
            //  (
            //      gMax
            //      (
            //          (
            //              mag
            //              (
            //                  D.internalField() - D.oldTime().internalField()
            //              )
            //          )
            //      ),
            //      SMALL
            //  )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledPressureDisplacementSolid::coupledPressureDisplacementSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    nonLinear_(true),
    DDDp_
    (
        IOobject
        (
            "DDDp",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector4("zero", dimless, vector4::zero)
    ),
    DU_
    (
        IOobject
        (
            "DU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    Dp_
    (
        IOobject
        (
            "Dp",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", Dp_.dimensions(), 0)
    ),
    Dpf_
    (
        IOobject
        (
            "Dpf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(Dp_)
    ),
    pf_
    (
        IOobject
        (
            "pf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", p_.dimensions(), 0)
    ),
    gradDp_
    (
        IOobject
        (
            "grad(Dp)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            dimPressure/dimLength,
            vector::zero
        ),
        zeroGradientFvPatchVectorField::typeName
    ),
    Df_
    (
        IOobject
        (
            "Df",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(D())
    ),
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
        mesh(),
        dimensionedVector
        (
            "zero",
            dimVelocity,
            vector::zero
        )
    ),
    DUf_
    (
        IOobject
        (
            "DUf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            dimVelocity,
            vector::zero
        )
    ),
    DDf_
    (
        IOobject
        (
            "DDf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            dimLength,
            vector::zero
        )
    ),
    curSf_
    (
        IOobject
        (
            "curSf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("curSf", dimArea, vector::zero)
    ),
    AD_
    (
        IOobject
        (
            "AD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimDensity/(dimTime*dimTime), 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    ADt_
    (
        IOobject
        (
            "ADt",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimDensity/(dimTime*dimTime), 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    phi_
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(DD()) & mesh().Sf()
    ),
    nonLinForceCorrection_
    (
        IOobject
        (
            "nonLinForceCorrection",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce, vector::zero)
    ),
    force_
    (
        IOobject
        (
            "force",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce, vector::zero)
    ),
    solidInterfaceTraction_
    (
        IOobject
        (
            "solidIntefaceTraction",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimPressure, vector::zero)
    ),
    A_
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimAcceleration, vector::zero)
    ),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    gradDDf_
    (
        IOobject
        (
            "grad(" + DD().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    tGradDDn_
    (
        IOobject
        (
            "tGradDDn",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    nGradDDn_
    (
        IOobject
        (
            "nGradDDn",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    impKf_
    (
        IOobject
        (
            "impKf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mechanical().impKf()
    ),
    rKappa_(1.0/mechanical().bulkModulus()()),
    nuCoeff_
    (
        IOobject
        (
            "nuCoeff",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1)
    ),
    debug_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "debug",
            false
        )
    ),
    K_
    (
        solidModelDict().lookupOrDefault<dimensionedScalar>
        (
            "K",
            dimensionedScalar("K", dimless/dimTime, 0)
        )
    ),
    relativeTol_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "solutionTolerance",
            solutionTol()
        )
    ),
    refPoints_
    (
        IOobject
        (
            "points",
            mesh().time().constant(),
            polyMesh::meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    consistentRhieChow_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "consistentRhieChow",
            false
        )
    ),
    refConfig_(true),
    composite_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "composite",
            false
        )
    ),
    stdDispGrad_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "stdDispGrad",
            true
        )
    ),
    coeff_("coeff", dimless, 1.0),
    delta_
    (
        IOobject
        (
            "delta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    weights_
    (
        IOobject
        (
            "weights",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimless, 0)
    ),
    convergDataFilePtr_(),
    writeConvergData_
    (
        solidModelDict().lookupOrDefault<Switch>
        (
            "writeConvergenceData",
            false
        )
    ),
    // Requred only if consisten Rhie-Chow correction is used
    interpDDfPtr_(NULL),
    interpDUfPtr_(NULL),
    interpDfPtr_(NULL),
    interpUfPtr_(NULL)
{
    DDisRequired();
    // DisRequired();

    if (runTime.timeIndex() == 0)
    {
        Info << "Initialize curSf, deltas and weights" << endl;

        curSf_ = mesh().Sf();
        weights_ = mesh().weights();
        #include "updateDelta.H"

        // AD_ should be calculated using initial impKf_
        #include "updateAD.H"
    }

    curSf_.oldTime();
    weights_.oldTime();
    delta_.oldTime();

    mesh().schemesDict().setFluxRequired(Dp_.name());

    p_.oldTime(); // Needed for incremental pressure eq.
    pf_.oldTime(); // Needed for incremental pressure eq.
    force_.oldTime();
    gradD().oldTime();
    gradDf_.oldTime();
    U().oldTime();

    // Calc. nuCoeff
    {
        const volScalarField G = mechanical().impK();

        if (gMax(rKappa_.internalField())/10 < SMALL)
        {
            rKappa_ = dimensionedScalar("0", dimless/dimPressure, 0);
        }

        const volScalarField GbyKappa = G*rKappa_;

        volScalarField nu =
            (1.0-2*GbyKappa/3)/(2*(1.0+GbyKappa/3));

        nuCoeff_ = 3*nu/(1+nu);

        Info << "nuCoeff, max: " << gMax(nuCoeff_.internalField())
             << ", avg: " << gAverage(nuCoeff_.internalField())
             << ", min: " << gMin(nuCoeff_.internalField()) << endl;
    }

    if
    (
        word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() + ')'))
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        U().oldTime().oldTime().oldTime();
        DD().oldTime().oldTime().oldTime();
        DU_.oldTime().oldTime().oldTime();
        DDf_.oldTime().oldTime().oldTime();
    }

    // If requested, create the convergence file
    if (writeConvergData_)
    {
        if (Pstream::master())
        {
            Info<< "Creating convergence data file" << endl;
            convergDataFilePtr_.set
            (
                new OFstream(runTime.path()/"convergData.dat")
            );
        }
    }

    // Check if linear model is specified
    {
        nonLinear_ =
            solidModelDict().lookupOrDefault<Switch>
            (
                "nonLinear",
                true
            );
    }

    // Prepare for consisten Rhie-Chow correction
    if(consistentRhieChow_)
    {
        interpDDfPtr_ =
            new surfaceVectorField
            (
                IOobject
                (
                    "interpDDf",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "zero",
                    dimLength,
                    vector::zero
                )
            );

        interpDfPtr_ =
            new surfaceVectorField
            (
                IOobject
                (
                    "interpDf",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "zero",
                    dimLength,
                    vector::zero
                )
            );

        interpDUfPtr_ =
            new surfaceVectorField
            (
                IOobject
                (
                    "interpDUf",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "zero",
                    dimVelocity,
                    vector::zero
                )
            );

        interpUfPtr_ =
            new surfaceVectorField
            (
                IOobject
                (
                    "interpUf",
                    runTime.timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector
                (
                    "zero",
                    dimVelocity,
                    vector::zero
                )
            );

        if
        (
            word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() + ')'))
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            interpDDf().oldTime().oldTime().oldTime();
            interpDf().oldTime().oldTime().oldTime();
            interpUf().oldTime().oldTime().oldTime();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool coupledPressureDisplacementSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    scalar initialResidual = 0;
    BlockSolverPerformance<vector4> solverPerfDp;
    blockLduMatrix::debug = 1;
    scalar res = 1.0;
    scalar maxRes = 0;
    scalar curConvergenceTolerance = solutionTol();

    scalar sumLocalContErr = 0;
    scalar globalContErr = 0;
    scalar cumulativeContErr = 0;

    enforceLinear() = false;

    dimensionedScalar backward("backward", dimless, 0);

    dimensionedScalar Cn("Cn", dimless, 3./2.);
    dimensionedScalar Co("Co", dimless, -2);
    dimensionedScalar Coo("Coo", dimless, 1./2.);

    if (runTime().timeIndex() == 1)
    {
        Cn.value()  =  1;
        Co.value()  = -1;
        Coo.value() =  0;
    }

    bool steadyState = false;

    if
    (
        word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() + ')'))
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        backward.value() = 1;
        coeff_.value() = 1;

        if (runTime().timeIndex() == 1)
        {
            // Euler
            backward.value() = 0;
            coeff_.value() = 1;
        }
    }
    else if
    (
        word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() +')'))
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        backward.value() = 0;
        coeff_.value() = 1;
    }
    else if
    (
        word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() +')'))
     == fv::CrankNicolsonDdtScheme<vector>::typeName
    )
    {
        backward.value() = 0;
        coeff_.value() = 2;

        if (runTime().timeIndex() == 1)
        {
            // Euler
            coeff_.value() = 1;
        }
    }
    else if
    (
        word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() +')'))
     == fv::steadyStateDdtScheme<vector>::typeName
    )
    {
        steadyState = true;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::evolve()\n"
        )   << "Unsupported temporal scheme: "
            << word(mesh().schemesDict().ddtScheme("ddt(" + DD().name() +')'))
            << abort(FatalError);
    }

    if (composite_)
    {
        if (runTime().timeIndex()%2)
        {
            backward.value() = 0;
            coeff_.value() = 2;

            if (runTime().timeIndex() == 1)
            {
                // Euler
                coeff_.value() = 1;
            }
        }
        else
        {
            backward.value() = 1;
            coeff_.value() = 1;
        }
    }

    dimensionedScalar rDeltaT = 1.0/runTime().deltaT();

    #include "updateADt.H"

    volScalarField rAD("rAD", nuCoeff_/(AD_ + ADt_));
    rAD.boundaryField().evaluateCoupled();

    if (consistentRhieChow_)
    {
        interpDf().oldTime().oldTime();
        interpUf().oldTime().oldTime();
        interpDUf().oldTime();
        interpDDf().oldTime();
    }

    curSf_.oldTime();
    delta_.oldTime();
    weights_.oldTime();

    do
    {
        if (blockLduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        // Store previous iteration to allow
        // under-relaxation and residual calculation
        Dp_.storePrevIter();
        DD().storePrevIter();
        nonLinForceCorrection_.storePrevIter();

        // Initialize the Dp block system (matrix,
        // source and reference to Dp)
        fvBlockMatrix<vector4> DDDpEqn(DDDp_);

        // Momentum equation
        {
            // Update gradient of displacement
            if (stdDispGrad_)
            {
                #include "stdDispGrad.H"
            }
            else
            {
                mechanical().grad(DD(), pointDD(), gradDD());
                mechanical().grad(DD(), pointDD(), gradDDf_);
            }

            gradD() = gradD().oldTime() + gradDD();
            gradDf_ = gradDf_.oldTime() + gradDDf_;

            // Calculate the stress using run-time
            // selectable mechanical law
            if (nonLinear_)
            {
                mechanical().correct(sigma());
            }

            Dp_.boundaryField().updateCoeffs();
            #include "calcFacePressureIncrement.H"

            // Calculate the stress using run-time
            // selectable mechanical law
            if (nonLinear_)
            {
                mechanical().correct(sigmaf_);
            }

            surfaceVectorField nf = mesh().Sf()/mesh().magSf();
            tGradDDn_ = ( (I-sqr(nf)) & (gradDDf_ & nf) );
            nGradDDn_ = ((nf*nf) & fvc::snGrad(DD()));

            #include "calcTraction.H"

            Info << "Assembling momentum equaton" << endl;
            momentumStabilisation().updateVector(DD(), &gradDD());

            // Construct momentum equation in total Lagrangian
            // form where gradients are calculated directly at the faces
            fvVectorMatrix DDEqn
            (
              - fvm::laplacian(mechanical().impK(), DD(), "laplacian(DDD,DD)")
              - fvc::div(impKf_*tGradDDn_*mesh().magSf())
              - fvc::div(nonLinForceCorrection_)
              - fvc::div(force_.oldTime())
              - momentumStabilisation().cellVector(&impKf_, true)
             == rho()*g()
            );

            // Add temporal terms
            if (!steadyState)
            {
                DDEqn +=
                    (1.0-backward)*sqr(coeff_)*rho()*rDeltaT*
                    (
                        fvm::Sp(rDeltaT, DD())
                      - U().oldTime()
                    )
                  + backward*rho()*rDeltaT*
                    (
                        Cn*rDeltaT*
                        (
                            fvm::Sp(Cn, DD())
                          + (Cn+Co)*D().oldTime()
                          + Coo*D().oldTime().oldTime()
                        )
                      + Co*U().oldTime()
                      + Coo*U().oldTime().oldTime()
                    )
                  - (coeff_-1.0)*fvc::div(force_.oldTime());
            }

            // Add damping
            if (K_.value() > SMALL)
            {
                DDEqn += K_*rho()*fvm::ddt(DD());
                DDEqn += K_*rho()*U().oldTime();
            }

            // Under-relax the linear system
            DDEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DDEqn);

            // Insert momentum equation
            DDDpEqn.insertEquation(0, DDEqn);

            // Add normal derivative of normal disp component term
            #include "completeMomentumEqAssembly.H"

            // Correction for blockCoupled boundary conditions
            #include "addBlockCoupledBC.H"

            // Add pressure gradient term
            BlockLduSystem<vector, vector> DpInDD(fvm::grad(Dp_));
            DpInDD *= nuCoeff_;
            DDDpEqn.insertBlockCoupling(0, 3, DpInDD, true);
        }

        // Move mesh to current configuration
        if (nonLinear_)
        {
            moveToCurConfig(false); // true = UL
        }

        surfaceScalarField rADf
        (
            "rADf",
            linear<scalar>(mesh()).interpolate(rAD)
        );

        // compute pressure gradient
        gradDp_ = fv::stabLeastSquaresGrad<scalar>(mesh()).grad(Dp_);
        gradDp_.correctBoundaryConditions();

        #include "calcPhi.H"

        Info << "Assembling presure equation" << endl;
        fvScalarMatrix DpEqn
        (
            fvm::Sp(rKappa_, Dp_)
          - fv::gaussLaplacianScheme<scalar, scalar>
            (
                mesh(),
                linear<scalar>(mesh()),
                fv::simpleUncorrectedSnGrad<scalar>(mesh())
            )
           .fvmLaplacian(rADf, Dp_)
          + fvc::div(phi_)
        );

        DDDpEqn.insertEquation(3, DpEqn);

        // Assemble and insert coupling terms
        BlockLduSystem<vector, scalar> DDInDp(fvm::UDiv(DD()));
        DDDpEqn.insertBlockCoupling(3, 0, DDInDp, false);

        // Solve the block matrix
        solverPerfDp = DDDpEqn.solve();

        // Retrieve solution
        DDDpEqn.retrieveSolution(0, DD().internalField());
        DDDpEqn.retrieveSolution(3, Dp_.internalField());

        DD().boundaryField().evaluateCoupled();
        DD().correctBoundaryConditions();
        Dp_.correctBoundaryConditions();

        // This is flux increment
        // Skewness correction for DD is already in phi
        phi_ += (linear<vector>(mesh()).interpolate(DD()) & mesh().Sf()) -
            rADf*fv::simpleUncorrectedSnGrad<scalar>(mesh()).snGrad(Dp_)*
            mesh().magSf();

        #include "solidContinuityErrs.H"

        // Total displacement and pressure
        D() = D().oldTime() + DD();
        p_ = p_.oldTime() + Dp_;

        // Calculate divergence-free face displacement
        {
            // Interpolated displacement icrement with skewness correction
            // (kewness correction is already added (see calcPhi.H))
            DDf_ += linear<vector>(mesh()).interpolate(DD());

            // Store interpolated displacement with skewness correction
            // for consisten Rhie-Chow correctin
            if (consistentRhieChow_)
            {
                interpDDf() = DDf_;
                // Calcualte total interpolated displacement
                interpDf() = interpDf().oldTime() + interpDDf();
            }

            // Conservative correction
            surfaceVectorField n = mesh().Sf()/mesh().magSf();
            surfaceScalarField DDfn = phi_/mesh().magSf();
            DDf_ += n*DDfn - n*(n & DDf_);

            // Calculate total conservative displacement
            Df_ = Df_.oldTime() + DDf_;
        }

        // Move mesh to reference configuration
        if (nonLinear_)
        {
            moveToRefConfig();
        }

        DD().relax();
        Dp_.relax();

        if (iCorr == 0)
        {
            initialResidual = mag(solverPerfDp.initialResidual());
        }

        // Calculate relative momentum residual
        res = residual(DD());

        if (res > maxRes)
        {
            maxRes = res;
        }

        curConvergenceTolerance = maxRes*relativeTol_;

        if (curConvergenceTolerance < solutionTol())
        {
            curConvergenceTolerance = solutionTol();
        }

        if
        (
            blockLduMatrix::debug
         || (iCorr % infoFrequency()) == 0
         || res < curConvergenceTolerance
         || maxIterReached() == nCorr()
        )
        {
            Info << "Corr " << iCorr << ", relative residual = "
                 << res << endl;
        }

        if (maxIterReached() == nCorr())
        {
            maxIterReached()++;
        }

        // Interpolate D to pointD
        mechanical().interpolate(DD(), pointDD(), false);

        // Total point displacement
        pointD() = pointD().oldTime() + pointDD();
    }
    while (res > curConvergenceTolerance && ++iCorr < nCorr());

    if (steadyState)
    {
        U() = dimensionedVector("zero", dimVelocity, vector::zero);
        DU_ = dimensionedVector("zero", dimVelocity, vector::zero);
        Uf_ = dimensionedVector("zero", dimVelocity, vector::zero);
        A_ = dimensionedVector("zero", dimAcceleration, vector::zero);
    }
    else
    {
        // Velocity increment
        DU_ = (1.0-backward)*
            (
                coeff_*fv::EulerDdtScheme<vector>(mesh()).fvcDdt(DD())
              - (coeff_-1.0)*DU_.oldTime()
            )
          + backward*
            (Cn*DD() + Co*DD().oldTime() + Coo*DD().oldTime().oldTime())/
            runTime().deltaT();

        // Velocity
        // U() = U().oldTime() + DU_;
        U() = (1.0-backward)*
            (
                coeff_*fv::EulerDdtScheme<vector>(mesh()).fvcDdt(D())
              - (coeff_-1.0)*U().oldTime()
            )
          + backward*
            (Cn*D() + Co*D().oldTime() + Coo*D().oldTime().oldTime())/
            runTime().deltaT();
        // U().correctBoundaryConditions();


        // Surface velocity increment
        DUf_ = (1.0-backward)*
            (
                coeff_*(DDf_ - DDf_.oldTime())/runTime().deltaT()
              - (coeff_-1.0)*DUf_.oldTime()
            )
          + backward*
            (Cn*DDf_ + Co*DDf_.oldTime() + Coo*DDf_.oldTime().oldTime())/
            runTime().deltaT();

        // Surface velocity (divergence free)
        Uf_ = (1.0-backward)*
            (
                coeff_*(Df_ - Df_.oldTime())/runTime().deltaT()
              - (coeff_-1.0)*Uf_.oldTime()
            )
          + backward*
            (Cn*Df_ + Co*Df_.oldTime() + Coo*Df_.oldTime().oldTime())/
            runTime().deltaT();

        // Interpolated velocity
        if (consistentRhieChow_)
        {
            interpDUf() =
                (1.0-backward)*
                (
                    coeff_*(interpDDf() - interpDDf().oldTime())/
                    runTime().deltaT()
                  - (coeff_-1.0)*interpDUf().oldTime()
                )
              + backward*
                (
                    Cn*interpDDf() + Co*interpDDf().oldTime() +
                    Coo*interpDDf().oldTime().oldTime()
                )
               /runTime().deltaT();

            interpUf() = //interpUf_.oldTime() + interpDUf_;
                (1.0-backward)*
                (
                    coeff_*(interpDf() - interpDf().oldTime())/
                    runTime().deltaT()
                  - (coeff_-1.0)*interpUf().oldTime()
                )
              + backward*
                (
                    Cn*interpDf() + Co*interpDf().oldTime() +
                    Coo*interpDf().oldTime().oldTime()
                )/
                runTime().deltaT();
        }

        // Acceleration
        A_ = (1.0-backward)*
            (
                coeff_*fv::EulerDdtScheme<vector>(mesh()).fvcDdt(U())
              - (coeff_-1.0)*A_.oldTime()
            )
          + backward*
            (Cn*U() + Co*U().oldTime() + Coo*U().oldTime().oldTime())/
            runTime().deltaT();
          // + backward*fv::backwardDdtScheme<vector>(mesh()).fvcDdt(U());
    }

    // Calculate the stress using run-time selectable mechanical law
    if (nonLinear_)
    {
        mechanical().correct(sigma());
    }
    else
    {
        volSymmTensorField epsilon = symm(gradD());
        volScalarField mu = mechanical().impK();
        sigma() = 2*mu*epsilon - nuCoeff_*p_*I;
    }

    // Calculate final cell-face forces
    #include "calcForce.H"

    // Print summary of residuals
    Info<< solverPerfDp.solverName() << ": Solving for " << DDDp_.name()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerfDp.initialResidual()
        << ", No outer iterations = " << iCorr << nl
        << " Max relative residual = " << maxRes
        << ", Relative residual = " << res << endl;

    blockLduMatrix::debug = 1;

    if (writeConvergData_)
    {
        if (!convergDataFilePtr_.empty())
        {
            convergDataFilePtr_() << runTime().value() << " "
                << res << " " << iCorr << endl;
        }
    }

    return true;
}


tmp<vectorField> coupledPressureDisplacementSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& mu = impKf_.boundaryField()[patchID];

    // Patch mechanical property
    const scalarField& nuCoeff = nuCoeff_.boundaryField()[patchID];

    // Patch gradient
    // const tensorField& pGradDD = gradDD().boundaryField()[patchID];
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch area vector (initial configuration)
    const vectorField& pSf = mesh().Sf().boundaryField()[patchID];

    // Patch tangential gradient of normal displacement component
    const vectorField& tGradDDn = tGradDDn_.boundaryField()[patchID];
    const vectorField& nGradDDn = nGradDDn_.boundaryField()[patchID];
    // vectorField tGradDDn = ((I - sqr(n)) & (pGradDD & n));
    // scalarField nGradDDn = ((n & pGradDD) & n);

    // Patch pressure increment
    const scalarField& Dp = Dp_.boundaryField()[patchID];

    // Old patch face forces
    const vectorField& oldPatchForce =
        force_.oldTime().boundaryField()[patchID];

    if (nonLinear_)
    {
        tensorField Ff = I + pGradD.T();
        tensorField invFf = inv(Ff);
        scalarField Jf = det(Ff);

        if (gMax(nuCoeff_.internalField()) > 0.5-SMALL)
        {
            // Incompressible case
            Jf = 1.0;
        }

        const scalarField magSf = patch.magSf();

        vectorField curSf = Jf*(invFf.T() & pSf);
        vectorField curnf = curSf/mag(curSf);

        vectorField oldSf =
            curSf_.oldTime().boundaryField()[patchID];

        // Patch nonlinear force correction
        const vectorField& nonLinForceCorrection =
            nonLinForceCorrection_.boundaryField()[patchID];

        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (traction - pressure*curnf)*mag(curSf)/(magSf*mu)
              - nGradDDn
              - tGradDDn
              + nuCoeff*Dp*oldSf/(magSf*mu)
              - nonLinForceCorrection/(magSf*mu)
              - oldPatchForce/(magSf*mu)
            )
        );
    }
    else
    {
        const scalarField magSf = patch.magSf();

        // Return patch snGrad
        return tmp<vectorField>
        (
            new vectorField
            (
                (traction - pressure*n)/mu
              - nGradDDn
              - tGradDDn
              + nuCoeff*Dp*n/mu
              - oldPatchForce/(mu*magSf)
            )
        );
    }
}


tmp<scalarField> coupledPressureDisplacementSolid::tractionBoundaryHydPressure
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& mu = impKf_.boundaryField()[patchID];

    // Patch mechanical property
    const scalarField& nuCoeff = nuCoeff_.boundaryField()[patchID];

    // Patch gradient
    // const tensorField& pGradDD = gradDD().boundaryField()[patchID];
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch area vector (initial configuration)
    const vectorField& pSf =
        mesh().Sf().boundaryField()[patchID];

    // Old patch face forces
    const vectorField& oldPatchForce =
        force_.oldTime().boundaryField()[patchID];

    if (nonLinear_)
    {
        tensorField Ff = I + pGradD.T();
        tensorField invFf = inv(Ff);
        scalarField Jf = det(Ff);

        if (gMax(nuCoeff_.internalField()) > 0.5-SMALL)
        {
            // Incompressible case
            Jf = 1.0;
        }

        vectorField curSf = Jf*(invFf.T() & pSf);
        vectorField curnf = curSf/mag(curSf);

        const scalarField magSf = patch.magSf();

        // Calc deviatoric part of the stress field
        // (subtract pressure)
        scalarField pf(n.size(), 0);
        {
            symmTensorField s =
                sigma().boundaryField()[patchID]
                // sigma() must be updated earlier
              + p_.boundaryField()[patchID]*I;

            pf = ((curnf & s) & curnf)
              - ((traction - pressure*curnf) & curnf);
        }

        return tmp<scalarField>
        (
            pf - p_.oldTime().boundaryField()[patchID]
        );
    }
    else
    {
        const scalarField magSf = patch.magSf();

        // Calculate normal derivative of normal comp of displacement
        // using conservative displacement increment at the boundary
        scalarField nGradDDn(magSf.size(), 0);
        {
            vectorField pDD =
                DDf_.boundaryField()[patchID];

            vectorField DDP =
                DD().boundaryField()[patchID].patchInternalField();

            tensorField gradDDP =
                gradDD().boundaryField()[patchID].patchInternalField();

            vectorField delta = mesh().boundary()[patchID].delta();
            vectorField k = ((I - sqr(n)) & delta);

            DDP += (k & gradDDP);

            nGradDDn = ((pDD - DDP) & n)*patch.deltaCoeffs();
        }

        // const scalarField nGradDDn =
        //     (nGradDDn_.boundaryField()[patchID] & n);

        // Return patch pressure
        return tmp<scalarField>
        (
            new scalarField
            (
                (
                    2*mu*nGradDDn
                  + ((oldPatchForce/magSf) & n)
                  - ((traction - pressure*n) & n)
                )/nuCoeff
            )
        );
    }
}


void coupledPressureDisplacementSolid::moveToCurConfig(const bool prevConfig)
{
    if (refConfig_)
    {
        Info << "Move mesh to current confguration" << endl;

        pointVectorField newPointD = pointD().oldTime();

        if (!prevConfig)
        {
            newPointD += 0.5*pointDD();
        }

        // vectorField newPoints = refPoints_;
        vectorField newPoints =
            refPoints_ + newPointD.internalField();

        mesh().movePoints(newPoints);
        mesh().V00();
        mesh().moving(false);
        mesh().changing(false);
        mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

        refConfig_ = false;

        curSf_ = mesh().Sf();
        weights_ = mesh().weights();
        #include "updateDelta.H"
    }
}

void coupledPressureDisplacementSolid::moveToRefConfig()
{
    if (!refConfig_)
    {
        Info << "Move mesh back to the initial confguration" << endl;
        mesh().movePoints(refPoints_);
        mesh().V00();
        mesh().moving(false);
        mesh().changing(false);
        mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

        refConfig_ = true;
    }
}


void coupledPressureDisplacementSolid::setTraction
(
    fvPatchVectorField& tractionPatch,
    const vectorField& traction
)
{
#ifndef OPENFOAMESIORFOUNDATION
    if
    (
        tractionPatch.type()
     == tractionPressureDisplacementFvPatchVectorField::typeName
    )
    {
        tractionPressureDisplacementFvPatchVectorField& patchD =
            refCast<tractionPressureDisplacementFvPatchVectorField>
            (
                tractionPatch
            );

        if (!patchD.frozenTraction())
        {
            patchD.traction() = traction;
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << tractionPatch.type()
            << " for patch " << tractionPatch.patch().name()
            << " should instead be type "
            << tractionPressureDisplacementFvPatchVectorField::typeName
            << abort(FatalError);
    }
#endif
}


void coupledPressureDisplacementSolid::setPressure
(
    fvPatchVectorField& pressurePatch,
    const scalarField& pressure
)
{
#ifndef OPENFOAMESIORFOUNDATION
    if
    (
        pressurePatch.type()
     == tractionPressureDisplacementFvPatchVectorField::typeName
    )
    {
        tractionPressureDisplacementFvPatchVectorField& patchD =
            refCast<tractionPressureDisplacementFvPatchVectorField>
            (
                pressurePatch
            );

        if (!patchD.frozenTraction())
        {
            patchD.pressure() = pressure;
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::setPressure\n"
            "(\n"
            "    fvPatchVectorField& pressurePatch,\n"
            "    const vectorField& pressure\n"
            ")"
        )   << "Boundary condition "
            << pressurePatch.type()
            << "for patch" << pressurePatch.patch().name()
            << " should instead be type "
            << tractionPressureDisplacementFvPatchVectorField::typeName
            << abort(FatalError);
    }
#endif
}


void coupledPressureDisplacementSolid::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    const vectorField patchTraction
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTraction)
    );

#ifdef OPENFOAMESIORFOUNDATION
    setTraction(solutionD().boundaryFieldRef()[patchID], patchTraction);
#else
    setTraction(solutionD().boundaryField()[patchID], patchTraction);
#endif

    solidInterfaceTraction_.boundaryField()[patchID] = patchTraction;
}


Foam::tmp<Foam::vectorField> coupledPressureDisplacementSolid::
faceZoneAcceleration
(
    const label interfaceI
) const
{
    // const volVectorField a(fvc::ddt(U()));
    // const volVectorField a(fvc::d2dt2(D()));

    return globalPatches()[interfaceI].patchFaceToGlobal
    (
        A_.boundaryField()[globalPatches()[interfaceI].patch().index()]
    );
}


void coupledPressureDisplacementSolid::updateTotalFields()
{
    solidModel::updateTotalFields();
    impKf_ = mechanical().impKf();

    // Update settings
    {
        K_ =
            solidModelDict().lookupOrDefault<dimensionedScalar>
            (
                "K",
                dimensionedScalar("K", dimless/dimTime, 0)
            );

        Info << "K: " << K_ << endl;
    }

    if (nonLinear_)
    {
        if (true)
        {
            vectorField newPoints =
                refPoints_ + pointD().internalField();
            mesh().movePoints(newPoints);

            if (debug_)
            {
                mesh().checkMesh(true);
            }

            curSf_ = mesh().Sf();
            weights_ = mesh().weights();
            #include "updateDelta.H"
            // #include "updateAD.H"


            // Report volume conservation error
            scalar curVolume = gSum(mesh().V().field());

            // Move mesh back to reference configuration
            mesh().movePoints(refPoints_);
            mesh().V00();
            mesh().moving(false);
            mesh().changing(false);
            mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

            scalar initVolume = gSum(mesh().V().field());
            scalar volumeRelError = (curVolume - initVolume)/initVolume;

            Info << "Initial total volume: " << initVolume << endl;
            Info << "Current total volume: " << curVolume << endl;
            Info << "Rel. volume error: " << volumeRelError << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
