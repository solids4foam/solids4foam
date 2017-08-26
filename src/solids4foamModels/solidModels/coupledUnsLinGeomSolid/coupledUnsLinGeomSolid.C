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

#include "coupledUnsLinGeomSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"
#include "BlockFvmDivSigma.H"

#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledUnsLinGeomSolid, 0);
addToRunTimeSelectionTable(physicsModel, coupledUnsLinGeomSolid, solid);
addToRunTimeSelectionTable(solidModel, coupledUnsLinGeomSolid, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledUnsLinGeomSolid::coupledUnsLinGeomSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    extendedMesh_(mesh()),
    solutionVec_
    (
        IOobject
        (
            "solutionVec",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        vectorField(extendedMesh_.nVariables(), vector::zero)
    ),
    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
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
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    muf_("muf", impKf_/3.5), // assuming a Poisson's ratio of 0.3
    lambdaf_("lambdaf", 1.5*muf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool coupledUnsLinGeomSolid::evolve()
{
    Info << "Evolving solid solver" << endl;

    int iCorr = 0;
    BlockSolverPerformance<vector> solverPerfD("undefined", "blockD");
    BlockLduMatrix<vector>::debug = 5;

    Info<< "Solving the momentum equation for D" << endl;

    // Momentum equation loop
    // In this case, a Hookean elastic constutive equation is assumed, and outer
    // corrections are performed is a different (possibily nonlinear )definition
    // of stress is used
    do
    {
        // Global coefficients are currently stored in the extendedMesh so we
        // must clear them out each time before constructing a new equation
        // Using a global patch may be a nicer solution
        //extendedMesh_.clearOut();
        extendedMesh_.clearOutGlobalCoeffs();

        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Create source vector for block matrix
        vectorField blockB(solutionVec_.size(), vector::zero);

        // Create block system
        BlockLduMatrix<vector> blockM(extendedMesh_);

        // Grab block diagonal and set it to zero
        Field<tensor>& d = blockM.diag().asSquare();
        d = tensor::zero;

        // Grab linear off-diagonal and set it to zero
        Field<tensor>& l = blockM.lower().asSquare();
        Field<tensor>& u = blockM.upper().asSquare();
        u = tensor::zero;
        l = tensor::zero;

        // Insert coefficients

        // For now we create separate matrices for each term of the three
        // diffusion terms in the momentum equation and add the contributions
        // together manually; a BlockFvMatrix class would help make this look
        // nicer

        // Laplacian
        // non-orthogonal correction is treated implicitly
        BlockLduMatrix<vector> blockMatLap =
            BlockFvm::laplacian(extendedMesh_, muf_, D(), blockB);

        // Laplacian transpose == div(mu*gradU.T())
        BlockLduMatrix<vector> blockMatLapTran =
            BlockFvm::laplacianTranspose(extendedMesh_, muf_, D(), blockB);

        // Laplacian trace == div(lambda*I*tr(gradU))
        BlockLduMatrix<vector> blockMatLapTrac =
            BlockFvm::laplacianTrace(extendedMesh_, lambdaf_, D(), blockB);

        // Add diagonal contributions
        d += blockMatLap.diag().asSquare();
        d += blockMatLapTran.diag().asSquare();
        d += blockMatLapTrac.diag().asSquare();

        // Add off-diagonal contributions
        u += blockMatLap.upper().asSquare();
        u += blockMatLapTran.upper().asSquare();
        u += blockMatLapTrac.upper().asSquare();
        l += blockMatLap.lower().asSquare();
        l += blockMatLapTran.lower().asSquare();
        l += blockMatLapTrac.lower().asSquare();

        // Add contribution for processor boundaries
        blockM.interfaces() = blockMatLap.interfaces();
        forAll(mesh().boundaryMesh(), patchI)
        {
            const word& patchType = mesh().boundaryMesh()[patchI].type();
            if (patchType == processorPolyPatch::typeName)
            {
                Field<tensor>& coupleUpper =
                    blockM.coupleUpper()[patchI].asSquare();

                coupleUpper = blockMatLap.coupleUpper()[patchI].asSquare();
                coupleUpper +=
                    blockMatLapTran.coupleUpper()[patchI].asSquare();
                coupleUpper +=
                    blockMatLapTrac.coupleUpper()[patchI].asSquare();
            }
        }

        // We manually add the boundary conditions equations
        // More thinking is required to get it to fit cleanly with the rest of
        // the block coupled machinery
        extendedMesh_.insertBoundaryConditions
        (
            blockM, blockB, muf_, lambdaf_, D()
        );


        // TESTING start
        // The implicit component assumes small strain Hookean elastic, so we
        // will explicitly remove this and then explicitly add on the divergence
        // of the actual stress
        /*
        // Store the div(sigma) coeffs, cells and boundary conditions
        vectorField AxDivSigma = vectorField(solutionVec_.size(), vector::zero);
        blockM.Amul(AxDivSigma, solutionVec_);

        // Store the div(sigma) coeffs, just cells
        vectorField AxDivSigmaCells =
            SubField<vector>(AxDivSigma, mesh().nCells(), 0);
            */
        // We will leave the boundary conditions untouched to be will set the
        // contribution to zero
        // for (int varI = mesh().nCells(); varI < Ax.size(); varI++)
        // {
        //     Ax[varI] = vector::zero;
        // }

        // Remove contribution to internal cells
        //blockB -= Ax;
        // TESTING end



        // Add terms temporal and gravity terms to the block matrix and source
        // Alos, we explicitly remove the implicit Hooke's law definition of
        // stress and explicitly add the definition of stress from the
        // mechanicalLaw
        extendedMesh_.addFvMatrix
        (
            blockM,
            blockB,
            rho()*fvm::d2dt2(D()) - rho()*g(),
            true
        );

        /*
        // Store the div(sigma) coeffs, cells and boundary conditions
        vectorField AxD2dt2RhoG =
            vectorField(solutionVec_.size(), vector::zero);
        blockM.Amul(AxD2dt2RhoG, solutionVec_);
        AxD2dt2RhoG -= AxDivSigma;

        // Store the just cell equations
        vectorField AxD2dt2RhoGCells =
            SubField<vector>(AxD2dt2RhoG, mesh().nCells(), 0);


        // TESTING start
        // Explicitly add the real definition of stress
        volVectorField newDivSigma =
            volVectorField("newDivSigma", fvc::div(mesh().Sf() & sigmaf_));
        newDivSigma.internalField() *= mesh().V();
        newDivSigma.write();
        const surfaceVectorField force = mesh().Sf() & sigmaf_;
        Info<< "Force[0] is: " << max(mag(force.boundaryField()[0])) << nl
            << "Force[1] is: " << max(mag(force.boundaryField()[1])) << endl;

        vectorField Ax = vectorField(solutionVec_.size(), vector::zero);
        blockM.Amul(Ax, solutionVec_);
        volVectorField AxField("Ax", 0.0*newDivSigma);
        AxField.internalField() = SubField<vector>(Ax, mesh().nCells(), 0);
        AxField.write();

        volVectorField blockBField("blockB", 0.0*newDivSigma);
        blockBField.internalField() =
            SubField<vector>(blockB, mesh().nCells(), 0);
        blockBField.write();

        // volVectorField AxD2dt2RhoGCellsField("AxD2dt2RhoG", 0.0*newDivSigma);
        // AxD2dt2RhoGCellsField.internalField() = AxD2dt2RhoGCells;
        // AxD2dt2RhoGCellsField.write();

        volVectorField AxDivSigmaCellsField("AxDivSigmaCells", 0.0*newDivSigma);
        AxDivSigmaCellsField.internalField() = AxDivSigmaCells;
        AxDivSigmaCellsField.write();

        volVectorField RhoTerm("RhoTerm", 0.0*newDivSigma);
        RhoTerm.internalField() = rho()*fvc::d2dt2(D()) - rho()*g();
        RhoTerm.internalField() *= mesh().V();
        RhoTerm.write();

        // Minus means we add to the right-hand side
        volVectorField contrib =
            volVectorField("contrib", AxDivSigmaCellsField - newDivSigma);
        const unallocLabelList& faceCells0 =
            mesh().boundaryMesh()[0].faceCells();
        forAll(contrib.boundaryField()[0], faceI)
        {
            contrib.internalField()[faceCells0[faceI]] = vector::zero;
        }
        const unallocLabelList& faceCells1 =
            mesh().boundaryMesh()[1].faceCells();
        forAll(contrib.boundaryField()[1], faceI)
        {
            contrib.internalField()[faceCells1[faceI]] = vector::zero;
        }
        contrib.write();
        // extendedMesh_.addFvSource
        // (
        //     blockB,
        //     contrib
        // );
*/
        // Block coupled solver call
        // solverPerfD =
        //     BlockLduSolver<vector>::New
        //     (
        //         D().name(),
        //         blockM,
        //         mesh().solutionDict().solver("blockD")
        //     )->solve(solutionVec_, blockB);

        // Under-relax the linear system
        if (mesh().solutionDict().relaxEquation("DEqn"))
        {
            FatalError
                << "Equation under-relaxation disabled"
                << abort(FatalError);

            // Info<< "Under-relaxing the equation" << endl;
            // blockM.relax
            // (
            //     solutionVec_,
            //     blockB,
            //     mesh().solutionDict().relaxationFactor("DEqn")
            // );
        }

        // Create the linear solver
        autoPtr<BlockLduSolver<vector> > solver =
            BlockLduSolver<vector>::New
            (
                D().name(),
                blockM,
                mesh().solutionDict().solver("blockD")
            );

        // Solve the linear system
        solver->solve(solutionVec_, blockB);

        // Transfer solution vector to D field
        extendedMesh_.copySolutionVector(solutionVec_, D());

        // Under-relax the field
        D().relax();

        // Update gradient of displacement
        mechanical().interpolate(D(), pointD());
        gradD() = fvc::grad(D(), pointD());
        gradDf_ = fvc::fGrad(D(), pointD());

        // We will call fvc::grad a second time as fixed displacement boundaries
        // need the patch internal field values to calculate snGrad, which
        // is then used to set snGrad on the gradD boundary
        D().correctBoundaryConditions();
        gradD() = fvc::grad(D(), pointD());
        gradDf_ = fvc::fGrad(D(), pointD());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigmaf_);

        // TESTING
        mechanical().correct(sigma());
    }
    while
    (
       !converged
        (
            iCorr,
            mag(solverPerfD.initialResidual()),
            solverPerfD.nIterations(),
            D()
        ) && ++iCorr < nCorr()
    );

    // Calculate cell stress
    mechanical().correct(sigma());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    return true;
}


tmp<vectorField> coupledUnsLinGeomSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& gradD = gradDf_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& sigma = sigmaf_.boundaryField()[patchID];

    // Patch unit normals
    const vectorField n = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (sigma - impK*gradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
