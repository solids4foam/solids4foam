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

#include "freeVibrationSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGradf.H"
#include "BlockFvmDivSigma.H"
#include "linearElastic.H"
#include "blockFixedDisplacementZeroShearFvPatchVectorField.H"
#include "BlockEigenProblemSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(freeVibrationSolid, 0);
addToRunTimeSelectionTable
(
    physicsModel, freeVibrationSolid, solid
);
addToRunTimeSelectionTable
(
    solidModel, freeVibrationSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeVibrationSolid::freeVibrationSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    nModes_(readInt(solidProperties().lookup("nModes"))),
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
    muf_
    (
        IOobject
        (
            "interpolate(mu)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    lambdaf_
    (
        IOobject
        (
            "interpolate(lambda)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    )
{
    // We will directly read the linearElastic mechanicalLaw
    const PtrList<mechanicalLaw>& mechLaws = mechanical();
    if (mechLaws.size() != 1)
    {
        FatalErrorIn
        (
            "freeVibrationSolid::"
            "freeVibrationSolid"
        )   << type() << " can currently only be used with a single material"
            << "\nConsider using one of the other solidModels."
            << abort(FatalError);
    }
    else if (!isA<linearElastic>(mechLaws[0]))
    {
        FatalErrorIn
        (
            "freeVibrationSolid::"
            "freeVibrationSolid"
        )   << type() << " can only be used with the linearElastic "
            << "mechanicalLaw" << nl
            << "Consider using one of the other linearGeometry solidModels."
            << abort(FatalError);
    }

    // Cast the mechanical law to a linearElastic mechanicalLaw
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set mu and lambda fields
    muf_ = mech.mu();
    lambdaf_ = mech.lambda();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool freeVibrationSolid::evolve()
{
    Info << "Evolving solid solver" << endl;

    BlockSolverPerformance<vector> solverPerfD("undefined", "blockD");
    BlockLduMatrix<vector>::debug = 5;

    Info<< "Solving the momentum equation for D" << endl;

    // Note: there is no momentum loop as everything is linear and we have
    // hard-coded the linear definition of stress

    // Global coefficients are currently stored in the extendedMesh so we
    // must clear them out each time before constructing a new equation
    // Using a global patch may be a nicer solution
    //extendedMesh_.clearOut();
    extendedMesh_.clearOutGlobalCoeffs();

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

    // Add terms temporal and gravity terms to the block matrix and source
    extendedMesh_.addFvMatrix
    (
        blockM,
        blockB,
        rho()*fvm::d2dt2(D()) - rho()*g(),
        true
    );

    // Create the EigenProblem linear solver
    BlockEigenProblemSolver solver
    (
        D().name(),
        blockM,
        mesh().solutionDict().solver("blockD")
    );

    // Solve the eigenproblem to calculate the eigenvalues and eigenvectors
    solver.solve(solutionVec_, blockB);

    // Store all requested modes
    for (int modeI = 0; modeI < nModes_; modeI++)
    {
        Info<< nl;

        // Retrieve the mode from the solver into the solutionVec
        solver.mode(modeI, solutionVec_);

        // Transfer solution vector to D field
        extendedMesh_.copySolutionVector(solutionVec_, D());

        // Update gradient of displacement
        mechanical().interpolate(D(), pointD());
        gradD() = fvc::grad(D(), pointD());

        // We will call fvc::grad a second time as fixed displacement boundaries
        // need the patch internal field values to calculate snGrad, which
        // is then used to set snGrad on the gradD boundary
        D().correctBoundaryConditions();
        gradD() = fvc::grad(D(), pointD());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Write the results to the current time-step
        runTime().write();

        // Move to the next time-step
        // Using const_cast is not the prettiest but vibration cases are
        // relatively exception, for now
        Info<< "    Writing mode " << (modeI + 1) << " to time "
            << runTime().value() << endl;

        if (modeI < (nModes_ - 1))
        {
            const_cast<Time&>(runTime())++;
        }
        else
        {
            FatalError
                << "done" << abort(FatalError);
        }
    }

    return true;
}


tmp<vectorField> freeVibrationSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& mu = muf_.boundaryField()[patchID];
    const scalarField& lambda = lambdaf_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField pN = patch.nf();

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - pN*pressure)
              - (pN & (pSigma - (2.0*mu + lambda)*pGradD))
            )/(2.0*mu + lambda)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
