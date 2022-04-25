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

#include "coupledNonLinGeomPressureDisplacementSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "fvBlockMatrix.H"
#include "skewCorrectionVectors.H"
#include "extrapolatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledNonLinGeomPressureDisplacementSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, coupledNonLinGeomPressureDisplacementSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void coupledNonLinGeomPressureDisplacementSolid::predict()
{
    Info<< "Sigma predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledNonLinGeomPressureDisplacementSolid::
coupledNonLinGeomPressureDisplacementSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimForce/dimArea,
            symmTensor::zero
        )
    ),
    gradDf_
    (
        IOobject
        (
            "grad(" + D().name() + ")f",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Ff_
    (
        IOobject
        (
            "Ff",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    Finvf_
    (
        IOobject
        (
            "Finvf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(Ff_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    Jf_
    (
        IOobject
        (
            "Jf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(Ff_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    pf_
    (
        IOobject
        (
            "pf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(p_)
    ),
    gradp_
    (
        IOobject
        (
            "grad(p)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            dimPressure/dimLength,
            vector::zero
        ),
        extrapolatedFvPatchVectorField::typeName
    ),
    Dp_
    (
        IOobject
        (
            "Dp",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector4("zero", dimless, vector4::zero)
    ),
    pressureRhieChowScaleFac_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "pressureRhieChowScaleFactor", 0.5
        )
    ),
    correction_
    (
        IOobject
        (
            "correction",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            impK_.dimensions()*dimArea,
            vector::zero
        )
    ),
    tGradDnf_
    (
        IOobject
        (
            "tGradDn_",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    nGradDnf_
    (
        IOobject
        (
            "nGradDn_",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    nonLinearForcef_
    (
        IOobject
        (
            "nonLinearForce",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce, vector::zero)
    ),
    linearForcef_
    (
        IOobject
        (
            "linearForce",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimForce, vector::zero)
    ),
    AU_
    (
        IOobject
        (
            "AU",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar
        (
            "zero",
            dimDensity/(dimTime*dimTime),
            0.0
        ),
        extrapolatedFvPatchScalarField::typeName
    ),
    HbyA_
    (
        IOobject
        (
            "HbyA",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            dimLength,
            vector::zero
        ),
        calculatedFvPatchVectorField::typeName
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
        fvc::interpolate(D()) & mesh().Sf()
    ),
    predictor_
    (
        solidModelDict().lookupOrDefault<Switch>("predictor", false)
    )
{
#   include "calcGradp.H"

    // Force p oldTime to be stored
    p_.oldTime();

    DisRequired();

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, solidProperties(), pRefCell, pRefValue);
    mesh().schemesDict().setFluxRequired(p_.name());

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // Calculate AU
    {
        const surfaceScalarField& deltaCoeffs = mesh().deltaCoeffs();
        const scalarField& deltaCoeffsI = deltaCoeffs.internalField();

        const surfaceScalarField& magS = mesh().magSf();
        const scalarField& magSI = magS.internalField();

        const unallocLabelList& own = mesh().owner();
        const unallocLabelList& nei = mesh().neighbour();

        const scalarField& impKfI = impKf_.internalField();

        scalarField& AUI = AU_.internalField();

        forAll(deltaCoeffsI, faceI)
        {
            AUI[own[faceI]] += impKfI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
            AUI[nei[faceI]] += impKfI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
        }

        forAll(deltaCoeffs.boundaryField(), patchI)
        {
            const scalarField& deltaCoeffsI =
                deltaCoeffs.boundaryField()[patchI];

            const scalarField& magSI = magS.boundaryField()[patchI];

            const labelList& faceCells =
                mesh().boundary()[patchI].faceCells();

            const scalarField& impKfI = impKf_.boundaryField()[patchI];

            forAll(deltaCoeffsI, faceI)
            {
                AUI[faceCells[faceI]] +=
                    0.5*impKfI[faceI]*deltaCoeffsI[faceI]*magSI[faceI];
            }
        }

        AUI /= mesh().V().field();

        AU_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool coupledNonLinGeomPressureDisplacementSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    if (predictor_)
    {
        predict();
    }

    // Disable default writing of linear solver residuals
    blockLduMatrix::debug = 0;

    int iCorr = 0;
    BlockSolverPerformance< VectorN<double, 4> > solverPerfDp;

    Info<< "Solving the momentum equation for D and p" << endl;

    // Loop around displacement and pressure equations
    do
    {
        p_.boundaryField().updateCoeffs();

        // Store fields for under-relaxation and residual calculation
        p_.storePrevIter();
        D().storePrevIter();
#       include "calcGradp.H"

        // Initialize the Dp block system (matrix, source and reference to
        // Dp)
        fvBlockMatrix<vector4> DpEqn(Dp_);

#       include "calcGradU.H"

        surfaceVectorField nf = mesh().Sf()/mesh().magSf();
        surfaceVectorField snGradD = fvc::snGrad(D());

        // Update components of surface normal gradient
        tGradDnf_ = (I - nf*nf) & (gradDf_ & nf);
        nGradDnf_ = snGradD & nf;

        // Compute surface force due to non-linear law
        nonLinearForcef_ = (Jf_*Finvf_.T() & mesh().Sf()) & sigmaf_;

        // Compute surface force due to linear law
        linearForcef_ =
              impKf_*snGradD*mesh().magSf()
            + impKf_*nGradDnf_*mesh().Sf()
            + impKf_*tGradDnf_*mesh().magSf()
            - pf_*mesh().Sf();

        // Momentum equation in terms of displacement
        fvVectorMatrix DEqn
        (
            fvm::d2dt2(rho(), D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(impKf_*tGradDnf_*mesh().magSf())
          + fvc::div(impKf_*nGradDnf_*mesh().Sf())
          + fvc::div(nonLinearForcef_ - linearForcef_)
          + rho()*g()
          - (gradp_ - fvc::grad(p_))
          + stabilisation().stabilisation(DD(), gradDD(), impK_)
        );

        // Under-relaxation the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Insert displacement equation into block system
        DpEqn.insertEquation(0, DEqn);

#if FOAMEXTEND
#       include "addBlockCoupledBC.H"
#endif

        // Update "flux" phi
        HbyA_ = D() + gradp_/AU_;
        HbyA_.correctBoundaryConditions();

#       include "calcPhi.H"

        // Store reciprocal of diagonal
        const volScalarField rAU = 1.0/AU_;
        const surfaceScalarField rAUf = fvc::interpolate(rAU);

        // Pressure equation
        fvScalarMatrix pEqn
        (
            fvm::Sp(rImpK_, p_)
          - fvm::laplacian(rAU, p_, "laplacian(Dp,p)")
         ==
          - fvc::div(phi_)
          + pressureRhieChowScaleFac_
           *(
                fvc::laplacian(rAUf, p_, "laplacian(Dp,p)")
              - fvc::div(rAUf*mesh().Sf() & fvc::interpolate(gradp_))
            )
        );

        label pRefCell = 0;
        scalar pRefValue = 0.0;
        pEqn.setReference(pRefCell, pRefValue);

        // Under-relaxation the linear system
        pEqn.relax();

        // Insert pressure equation into block system
        DpEqn.insertEquation(3, pEqn);

        // Insert coupling terms into block system
        {
            // Calculate grad p coupling matrix
            BlockLduSystem<vector, vector> pInD(fvm::grad(p_));

            // Calculate div D coupling
            BlockLduSystem<vector, scalar> DInp(fvm::UDiv(D()));

            DpEqn.insertBlockCoupling(0, 3, pInD, true);
            DpEqn.insertBlockCoupling(3, 0, DInp, false);
        }

        // Add laplacian of normal displacement component
#       include "finalizeMomentumEqn.H"

        // Solve the block matrix
        solverPerfDp = DpEqn.solve();

        // Retrieve solution
        DpEqn.retrieveSolution(0, D().internalField());
        DpEqn.retrieveSolution(3, p_.internalField());

        D().correctBoundaryConditions();
        p_.correctBoundaryConditions();

        // Update phi
        phi_ += (fvc::interpolate(D()) & mesh().Sf()) + pEqn.flux();

        // Fixed or adaptive field under-relaxation
        p_.relax();

        // Update increment of displacement
        DD() = D() - D().oldTime();

        mechanical().interpolate(D(), pointD());

        // Update gradient of displacement
        mechanical().grad(D(), pointD(), gradD());
        mechanical().grad(D(), pointD(), gradDf_);

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        Ff_ = I + gradDf_.T();

        // Inverse of the deformation gradient
        Finvf_ = inv(Ff_);

        // Jacobian of the deformation gradient
        Jf_ = det(Ff_);

        // Update the gradient of pressure
        gradp_ = fvc::grad(p_);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
        mechanical().correct(sigmaf_);
    }
    while
    (
        !converged
        (
            iCorr,
            cmptMax(solverPerfDp.initialResidual()),
            solverPerfDp.nIterations(),
            D()
        )
     && ++iCorr < nCorr()
    );

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    mechanical().correct(sigma());

    blockLduMatrix::debug = 1;

    return true;
}


tmp<vectorField> coupledNonLinGeomPressureDisplacementSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impKf_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch hydrostatic pressure
    const scalarField& pP = pf_.boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finvf_.boundaryField()[patchID];

    // Patch total Jacobian
    const scalarField& J = Jf_.boundaryField()[patchID];

    // Patch unit normals (deformed configuration)
    const vectorField nCurrent = J*Finv.T() & n;

    // Patch current normal
    const vectorField curSf = J*(Finv.T() & patch.Sf());
    const scalarField magCurSf = mag(curSf);

    // Patch tang. and normal comp. of surface normal gradient
    const vectorField& tGradDn = tGradDnf_.boundaryField()[patchID];
    //const scalarField& nGradDn = nGradDnf_.boundaryField()[patchID];

    const vectorField& pNonLinForce = nonLinearForcef_.boundaryField()[patchID];
    const vectorField& pLinearForce = linearForcef_.boundaryField()[patchID];

    // Difference between nonlinear and linear traction
    // Note: we are using traction rather than force
    const vectorField deltaNonlinTrac =
        pNonLinForce/magCurSf - pLinearForce/mag(patch.Sf());

    // Return patch snGrad
    return tmp<vectorField>
        (
            new vectorField
            (
                (
                    (I - sqr(n))
                  & (
                        traction - pressure*nCurrent // T we want to apply
                      - deltaNonlinTrac
                    )
                )*rImpK
              - tGradDn
              // Normal gradient component
              + n
               *(
                   pP
                 + (
                       n
                     & (
                           traction - pressure*nCurrent
                         - deltaNonlinTrac
                       )
                   )
                )/(2.0*impK)
            )
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
