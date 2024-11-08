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

#include "linGeomTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "symmetryFvPatchFields.H"

#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linGeomTotalDispSolid, 0);
addToRunTimeSelectionTable(solidModel, linGeomTotalDispSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void linGeomTotalDispSolid::predict()
{
    Info<< "Applying linear predictor to D" << endl;

    // Predict D using previous time steps
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


void linGeomTotalDispSolid::enforceTractionBoundaries
(
    surfaceVectorField& traction,
    const volVectorField& D,
    const surfaceVectorField& n
) const
{
    // Enforce traction conditions
    forAll(D.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            const solidTractionFvPatchVectorField& tracPatch =
                refCast<const solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            const vectorField& nPatch = n.boundaryField()[patchI];

            traction.boundaryFieldRef()[patchI] =
                tracPatch.traction() - nPatch*tracPatch.pressure();
        }
        else if
        (
            isA<fixedDisplacementZeroShearFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
         || isA<symmetryFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            // Unit normals
            const vectorField& nPatch = n.boundaryField()[patchI];

            // Set shear traction to zero
            traction.boundaryFieldRef()[patchI] =
                sqr(nPatch) & traction.boundaryField()[patchI];
        }
    }
}


bool linGeomTotalDispSolid::evolveImplicitSegregated()
{
    Info<< "Evolving solid solver using an implicit segregated approach"
        << endl;

    if (predictor_ && newTimeStep())
    {
        predict();
    }

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        SolverPerformance<vector> solverPerfD;
        SolverPerformance<vector>::debug = 0;
#else
        lduSolverPerformance solverPerfD;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the momentum equation for D" << endl;

        // Momentum equation loop
        do
        {
            // Store fields for under-relaxation and residual calculation
            D().storePrevIter();

            // Linear momentum equation total displacement form
            fvVectorMatrix DEqn
            (
                rho()*fvm::d2dt2(D())
             == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
              - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
              + fvc::div(sigma(), "div(sigma)")
              + rho()*g()
              + stabilisation().stabilisation(D(), gradD(), impK_)
              + fvOptions()(ds_, D())
            );

            // Add damping
            if (dampingCoeff().value() > SMALL)
            {
                DEqn += dampingCoeff()*rho()*fvm::ddt(D());
            }

            // Under-relaxation the linear system
            DEqn.relax();

            // Enforce any cell displacements
            solidModel::setCellDisps(DEqn);

            // Solve the linear system
            solverPerfD = DEqn.solve();

            // Fixed or adaptive field under-relaxation
            relaxField(D(), iCorr);

            // Update increment of displacement
            DD() = D() - D().oldTime();

            // Update gradient of displacement
            mechanical().grad(D(), gradD());

            // Testing new grad scheme start

            // Here we will calculate cell-centred gradient using
            // high-order disretisation approach
            // STEPS:
            // 1. Loop over all cells
            //   1.1 For each cell construct stencil
            //
            // 2. Loop over all cells
            //   2.1 For each cell make interpolation coefficients
            //
            // 3. Explicitly calculate gradient at cell centres unsing cell
            //    stencil and corresponting interpolation coeffs for each cell
            //    in stencil

            //
            // 1. Step - make stencils
            //
            //          I'm using layer approach as it has the best condition
            //          number for LS matrix. Some authors use sphere and put
            //          in stencil all cells inside the sphere but that produces
            //          higher values of condition number at structured mesh

            const label maxStencilSize = 70;
            const label nLayers = 2;
            const labelListList& cellCells = mesh().cellCells();

            List<DynamicList<label>> lsStencil(mesh().nCells());
            forAll(lsStencil, cellI)
            {
                lsStencil[cellI].setCapacity(maxStencilSize);

                DynamicList<label>& curStencil = lsStencil[cellI];
                const labelList& curCellCells = cellCells[cellI];

                labelHashSet stencilCells;
                labelHashSet prevLayer;

                // Add first layer of cells
                forAll(curCellCells, cI)
                {
                    stencilCells.insert(curCellCells[cI]);
                    prevLayer.insert(curCellCells[cI]);
                }

                // Remaining layers of cells
                for (int layerI = 1; layerI < nLayers; layerI++)
                {
                    labelList prevLayerCells(prevLayer.toc());
                    labelHashSet curLayer;

                    // Loop over previous layer and add one level of
                    // layer neighbours
                    for (const label cI : prevLayerCells)
                    {
                        const labelList& cellINei= mesh().cellCells()[cI];
                        forAll (cellINei, nei)
                        {
                            if (!stencilCells.found(cellINei[nei]))
                            {
                                curLayer.insert(cellINei[nei]);
                            }
                        }
                    }

                    // Now we have curent layer which we need to add to stencil
                    // and current layer will now be previous layer for next
                    // loop
                    prevLayer.clear();
                    prevLayer = curLayer;
                    stencilCells.merge(curLayer);
                }

                if (Pstream::parRun())
                {
                    notImplemented("not implemented for parallel run");
                }

                // Neighbours of first layer will store cellI in stencil so we
                // need to remove it
                stencilCells.erase(cellI);

                curStencil.append(stencilCells.toc());
            }

            forAll(lsStencil, cellI)
            {
                // Note: lsStencil is sorted but I do not see problem in that
                lsStencil[cellI].shrink();
            }

            //
            // 2. Step - calculate interpolation coefficients
            //

            // Order of interpolation
            const label N = 1;

            // Number of terms in Taylor expression
            // 1 for zero order, 4 for 1 order, 10 for second order
            label Np = factorial(N+3)/(factorial(N)*factorial(3));

            // Kernel shape parameter
            const label k = 6;

            // Cells condition numbers
            volScalarField condNumber
            (
                IOobject
                (
                    "condNumber",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("0", dimless, Zero)
            );

            List<DynamicList<scalar>> c(mesh().nCells());
            List<DynamicList<scalar>> cx(mesh().nCells());
            List<DynamicList<scalar>> cy(mesh().nCells());
            List<DynamicList<scalar>> cz(mesh().nCells());

            // Cell boundary faces global index
            List<DynamicList<label>> cellBoundaryFaces(mesh().nCells());

            const polyBoundaryMesh& boundaryMesh = mesh().boundaryMesh();

            forAll(boundaryMesh, patchI)
            {
                if (boundaryMesh[patchI].type() != emptyPolyPatch::typeName)
                {
                    const labelUList& faceCells = boundaryMesh[patchI].faceCells();

                    forAll(faceCells, faceI)
                    {
                        const label cellID = faceCells[faceI];
                        const label gI = faceI+boundaryMesh[patchI].start();
                        cellBoundaryFaces[cellID].append(gI);
                    }
                }
            }

            forAll(cellBoundaryFaces, cellI)
            {
                cellBoundaryFaces[cellI].shrink();
            }


            forAll(lsStencil, cellI)
            {
                DynamicList<label>& curStencil = lsStencil[cellI];

                // Find max distance in this stencil
                scalar maxDist = 0.0;
                forAll(curStencil, cI)
                {
                    const label neiCellID = curStencil[cI];
                    const scalar d =mag(mesh().C()[neiCellID]-mesh().C()[cellI]);
                    if ( d > maxDist)
                    {
                        maxDist = d;
                    }
                }

                // Loop over neighbours and construct matrix Q
                const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

                // For now I will use matrix format from Eigen/Dense library
                Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(Np, Nn);

                // Check to avoid Eigen error
                if (Nn < Np)
                {
                    FatalErrorInFunction
                        << "Interpolation stencil needs to be bigger than the "
                       "number of elements in Taylor order!"
                        << exit(FatalError);
                }

                // Loop over cells in stencil, each cell have its corresponding
                // row in Q matrix
                for(label cI=0; cI<Nn; cI++)
                //forAll(curStencil, cI)
                {
                    vector neiC = vector::zero;

                    if (cI < curStencil.size())
                    {
                        const label neiCellID = curStencil[cI];
                        neiC = mesh().C()[neiCellID];
                    }
                    else
                    {
                        // For boundary cells we need to add boundary face as
                        // neigbour
                        const label i = cI-curStencil.size();
                        const label globalFaceID = cellBoundaryFaces[cellI][i];
                        neiC = mesh().Cf()[globalFaceID];
                    }

                    const vector C = mesh().C()[cellI];

                    // Add linear interpolation part N=1
                    if ( N > 0 )
                    {
                        Q(0, cI) = 1;
                        Q(1, cI) = neiC.x() - C.x();
                        Q(2, cI) = neiC.y() - C.y();
                        Q(3, cI) = neiC.z() - C.z();
                    }
                    // Add quadratic interpolation part N=2
                    if ( N > 1 )
                    {
                        Q(4, cI) = (1.0/2.0)*pow(neiC.x() - C.x(), 2);
                        Q(5, cI) = (1.0/2.0)*pow(neiC.y() - C.y(), 2);
                        Q(6, cI) = (1.0/2.0)*pow(neiC.z() - C.z(), 2);
                        Q(7, cI) = (neiC.x() - C.x())*(neiC.y() - C.y());
                        Q(8, cI) = (neiC.x() - C.x())*(neiC.z() - C.z());
                        Q(9, cI) = (neiC.y() - C.y())*(neiC.z() - C.z());
                    }
                    // Add cubic interpolation part N=3
                    if ( N > 2)
                    {
                        // This have 20 terms
                        notImplemented("Orders higher that quadratic not implemented");
                    }
                }

                Eigen::MatrixXd W = Eigen::MatrixXd::Zero(Nn, Nn);

                //forAll(curStencil, cI)
                for(label cI=0; cI<Nn; cI++)
                {
                    vector neiC = vector::zero;

                    if (cI < curStencil.size())
                    {
                        const label neiCellID = curStencil[cI];
                        neiC = mesh().C()[neiCellID];
                    }
                    else
                    {
                        // For boundary cells we need to add boundary face as
                        // neigbour
                        const label i = cI-curStencil.size();
                        const label globalFaceID = cellBoundaryFaces[cellI][i];
                        neiC = mesh().Cf()[globalFaceID];
                    }

                    //const label neiCellID = curStencil[cI];
                    //const vector neiC = mesh().C()[neiCellID];

                    const vector C = mesh().C()[cellI];

                    const scalar d = mag(neiC-C);

                    // Smoothing length
                    const scalar dm = 2 * maxDist;

                    // Weight using radially symmetric exponential function
                    const scalar sqrK = -pow(k,2);
                    const scalar w =  (Foam::exp(pow(d/dm, 2) * sqrK) - Foam::exp(sqrK)) / (1 - exp(sqrK));

                    W(cI, cI) = w;
                }

                // Now when we have W and Q, next step is QR decomposition
                Eigen::MatrixXd sqrtW = W.cwiseSqrt();
                Eigen::MatrixXd Qhat = Q * sqrtW;
                Eigen::HouseholderQR<Eigen::MatrixXd> qr(Qhat.transpose());

                // Q and R matrices
                Eigen::MatrixXd O = qr.householderQ();
                Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();

                // B hat
                Eigen::MatrixXd Bhat = sqrtW * Eigen::MatrixXd::Identity(W.rows(), W.cols());

                // Slice Rbar and Qbar, as we do not need full matrix
                Eigen::MatrixXd Rbar = R.topLeftCorner(Np, Np);
                Eigen::MatrixXd Qbar = O.leftCols(Np);

                // Solve to get A
                Eigen::MatrixXd A = Rbar.colPivHouseholderQr().solve(Qbar.transpose() * Bhat);

                // To be aware of interpolation accuracy we need to control the
                // condition number
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(Rbar, Eigen::ComputeFullU | Eigen::ComputeFullV);
                Eigen::VectorXd singularValues = svd.singularValues();
                condNumber[cellI] = singularValues(0) / (singularValues(singularValues.size() - 1) + VSMALL);

                c[cellI].setCapacity(A.cols());
                cx[cellI].setCapacity(A.cols());
                cy[cellI].setCapacity(A.cols());
                cz[cellI].setCapacity(A.cols());

                Eigen::RowVectorXd cRow = A.row(0);
                Eigen::RowVectorXd cxRow = A.row(1);
                Eigen::RowVectorXd cyRow = A.row(2);
                Eigen::RowVectorXd czRow = A.row(3);

                for (label i=0; i < A.cols(); ++i)
                {
                    c[cellI].append(cRow(i));
                    cx[cellI].append(cxRow(i));
                    cy[cellI].append(cyRow(i));
                    cz[cellI].append(czRow(i));
                }

                c[cellI].shrink();
                cx[cellI].shrink();
                cy[cellI].shrink();
                cz[cellI].shrink();
            }

            List<DynamicList<scalar>> interpCoeffs(mesh().nCells());
            List<DynamicList<vector>> interpGradCoeffs(mesh().nCells());

            forAll(interpCoeffs, cellI)
            {
                DynamicList<label>& curStencil = lsStencil[cellI];
                const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

                for(label I=0; I<Nn; I++)
                //forAll(curStencil, nei)
                {
                    interpCoeffs[cellI].append(c[cellI][I]);
                    interpGradCoeffs[cellI].append(vector(cx[cellI][I], cy[cellI][I], cz[cellI][I]));
                }

                interpCoeffs[cellI].shrink();
                interpGradCoeffs[cellI].shrink();
            }

            // //
            // // Step 3: Interpolate D() to get gradient of D
            // //

            volTensorField hoGradD
            (
                IOobject
                (
                    "hoGrad(D)",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedTensor("0", dimless, Zero)
            );

            forAll(lsStencil, cellI)
            {
                DynamicList<label>& curStencil = lsStencil[cellI];
                const label Nn = curStencil.size() + cellBoundaryFaces[cellI].size();

                // Loop over stencil and multiply stencil cell values with
                // corresponding interpolation coefficient
                for(label cI=0; cI<Nn; cI++)
                //forAll(curStencil, cI)
                {
                    if (cI < curStencil.size())
                    {
                        hoGradD[cellI] += D()[curStencil[cI]] * interpGradCoeffs[cellI][cI];
                    }
                    else
                    {
                        const label i = cI-curStencil.size();
                        const label globalFaceID = cellBoundaryFaces[cellI][i];


                        vector boundaryD = vector::zero;

                        forAll(boundaryMesh, patchI)
                        {
                            const label start = boundaryMesh[patchI].start();
                            const label nFaces = boundaryMesh[patchI].nFaces();


                            if (globalFaceID >= start && globalFaceID < start+nFaces)
                            {
                                const label k = globalFaceID-start;
                                boundaryD = D().boundaryField()[patchI][k];
                            }
                        }

                        hoGradD[cellI] += boundaryD * interpGradCoeffs[cellI][cI];
                    }
                }
            }

            hoGradD.write();
            condNumber.write();

            scalar cellRelErrorXX = 0.0;
            scalar cellRelErrorXY = 0.0;
            scalar cellRelErrorYY = 0.0;
            scalar cellRelErrorXZ = 0.0;
            scalar cellRelErrorYZ = 0.0;
            scalar cellRelErrorZZ = 0.0;

            forAll(hoGradD, cellI)
            {
                cellRelErrorXX = abs((hoGradD[cellI].xx()-gradD()[cellI].xx())/(gradD()[cellI].xx()+SMALL))*100;
                cellRelErrorYY = abs((hoGradD[cellI].yy()-gradD()[cellI].yy())/(gradD()[cellI].yy()+SMALL))*100;
                cellRelErrorXY = abs((hoGradD[cellI].xy()-gradD()[cellI].xy())/(gradD()[cellI].xy()+SMALL))*100;
                cellRelErrorXZ = abs((hoGradD[cellI].xz()-gradD()[cellI].xz())/(gradD()[cellI].xz()+SMALL))*100;
                cellRelErrorYZ = abs((hoGradD[cellI].yz()-gradD()[cellI].yz())/(gradD()[cellI].yz()+SMALL))*100;
                cellRelErrorZZ = abs((hoGradD[cellI].zz()-gradD()[cellI].zz())/(gradD()[cellI].zz()+SMALL))*100;
            }
            cellRelErrorXX /= mesh().nCells();
            cellRelErrorYY /= mesh().nCells();
            cellRelErrorXZ /= mesh().nCells();
            cellRelErrorXY /= mesh().nCells();
            cellRelErrorYZ /= mesh().nCells();
            cellRelErrorZZ /= mesh().nCells();

            Info << "Average error for XX component" << cellRelErrorXX << endl;
            Info << "Average error for YY component" << cellRelErrorYY << endl;
            Info << "Average error for ZZ component" << cellRelErrorZZ << endl;
            Info << "Average error for XY component" << cellRelErrorXY << endl;
            Info << "Average error for XZ component" << cellRelErrorXZ << endl;
            Info << "Average error for YZ component" << cellRelErrorYZ << endl;
            Info<<endl;
            // Testing new grad schene end

            // Update gradient of displacement increment
            gradDD() = gradD() - gradD().oldTime();

            // Update the momentum equation inverse diagonal field
            // This may be used by the mechanical law when calculating the
            // hydrostatic pressure
            const volScalarField DEqnA("DEqnA", DEqn.A());

            // Calculate the stress using run-time selectable mechanical law
            mechanical().correct(sigma());
        }
        while
        (
            !converged
            (
                iCorr,
#ifdef OPENFOAM_NOT_EXTEND
                mag(solverPerfD.initialResidual()),
                cmptMax(solverPerfD.nIterations()),
#else
                solverPerfD.initialResidual(),
                solverPerfD.nIterations(),
#endif
                D()
            )
         && ++iCorr < nCorr()
        );

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), gradD(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (solidModel::mesh().update());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


bool linGeomTotalDispSolid::evolveSnes()
{
    Info<< "Solving the momentum equation for D using PETSc SNES" << endl;

    // Update D boundary conditions
    D().correctBoundaryConditions();

    // Solution predictor
    if (predictor_ && newTimeStep())
    {
        predict();

        // Map the D field to the SNES solution vector
        foamPetscSnesHelper::mapSolutionFoamToPetsc();
    }

    // Solve the nonlinear system and check the convergence
    foamPetscSnesHelper::solve();

    // Retrieve the solution
    // Map the PETSc solution to the D field
    foamPetscSnesHelper::mapSolutionPetscToFoam();

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());

    // Increment of displacement
    DD() = D() - D().oldTime();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    return true;
}


bool linGeomTotalDispSolid::evolveExplicit()
{
    if (time().timeIndex() == 1)
    {
        Info<< "Solving the solid momentum equation for D using an explicit "
            << "approach" << nl
            << "Simulation Time, Clock Time, Max Stress" << endl;
    }

    physicsModel::printInfo() = bool
    (
        time().timeIndex() % infoFrequency() == 0
     || mag(time().value() - time().endTime().value()) < SMALL
    );

    if (physicsModel::printInfo())
    {
        Info<< time().value() << " " << time().elapsedClockTime()
            << " " << max(mag(sigma())).value() << endl;

        physicsModel::printInfo() = false;
    }

    // Take references for brevity and efficiency
    const fvMesh& mesh = solidModel::mesh();
    volVectorField& D = solidModel::D();
    volTensorField& gradD = solidModel::gradD();
    volVectorField& U = solidModel::U();
    volSymmTensorField& sigma = solidModel::sigma();
    const volScalarField& rho = solidModel::rho();

    // Central difference scheme

    // Take a reference to the current and previous time-step
    const dimensionedScalar& deltaT = time().deltaT();
    //const dimensionedScalar& deltaT0 = time().deltaT0();

    // Compute the velocity
    // Note: this is the velocity at the middle of the time-step
    //pointU_ = pointU_.oldTime() + 0.5*(deltaT + deltaT0)*pointA_.oldTime();
    U = U.oldTime() + deltaT*A_.oldTime();

    // Compute displacement
    D = D.oldTime() + deltaT*U;

    // Enforce boundary conditions on the displacement field
    D.correctBoundaryConditions();

    if (solidModel::twoD())
    {
        // Remove displacement in the empty directions
        forAll(mesh.geometricD(), dirI)
        {
            if (mesh.geometricD()[dirI] < 0)
            {
                D.primitiveFieldRef().replace(dirI, 0.0);
            }
        }
    }

    // Update gradient of displacement
    mechanical().grad(D, gradD);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma);

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Calculate the traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // Solve the momentum equation for acceleration
    A_ = fvc::div(mesh.magSf()*traction)/rho
       + g()
       - dampingCoeff()*fvc::ddt(D);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linGeomTotalDispSolid::linGeomTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    foamPetscSnesHelper
    (
        fileName
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        ),
        D(),
        solidModel::twoD(),
        solidModelDict().lookupOrDefault<Switch>("stopOnPetscError", true),
        bool(solutionAlg() == solutionAlgorithm::PETSC_SNES)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    A_
    (
        IOobject
        (
            "A",
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength/pow(dimTime, 2), vector::zero)
    ),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false)),
    ds_
    (
        IOobject
        (
            "ds",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("ds", (dimForce/dimVolume)/dimVelocity, 1.0)
    )
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    // For consistent restarts, we will calculate the gradient field
    D().correctBoundaryConditions();
    D().storePrevIter();
    mechanical().grad(D(), gradD());

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // Check the gradScheme
    const word gradDScheme
    (
        mesh().gradScheme("grad(" + D().name() +')')
    );

    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        if (gradDScheme != "leastSquaresS4f")
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithm" << abort(FatalError);
        }

        // Set extrapolateValue to true for solidTraction boundaries
        forAll(D().boundaryField(), patchI)
        {
            if
            (
                isA<solidTractionFvPatchVectorField>
                (
                    D().boundaryField()[patchI]
                )
            )
            {
                Info<< "    Setting `extrapolateValue` to `true` on the "
                    << mesh().boundary()[patchI].name() << " patch of the D "
                    << "field" << endl;

                solidTractionFvPatchVectorField& tracPatch =
                    refCast<solidTractionFvPatchVectorField>
                    (
                        D().boundaryFieldRef()[patchI]
                    );

                tracPatch.extrapolateValue() = true;
            }
        }
    }
    else if (solutionAlg() != solutionAlgorithm::EXPLICIT)
    {
        if (gradDScheme == "leastSquaresS4f")
        {
            FatalErrorIn(type() + "::" + type())
                << "The `leastSquaresS4f` gradScheme should only be used for "
                << "`grad(D)` when using the "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " and "
                << solidModel::solutionAlgorithmNames_
                   [
                       solidModel::solutionAlgorithm::PETSC_SNES
                   ]
                << " solution algorithms" << abort(FatalError);
        }

    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void linGeomTotalDispSolid::setDeltaT(Time& runTime)
{
    if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        // Max wave speed in the domain
        const scalar waveSpeed = max
        (
            Foam::sqrt(mechanical().impK()/mechanical().rho())
        ).value();

        // deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
        // In the current discretisation, information can move two cells per
        // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
        // calculating the required stable time-step
        // i.e. deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
        // For safety, we should use a time-step smaller than this e.g. Abaqus uses
        // stableTimeStep/sqrt(2): we will default to this value
        const scalar requiredDeltaT =
            1.0/
            gMax
            (
                DimensionedField<scalar, Foam::surfaceMesh>
                (
                    mesh().surfaceInterpolation::
                    deltaCoeffs().internalField()
                   *waveSpeed
                )
            );

        // Lookup the desired Courant number
        const scalar maxCo =
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

        const scalar newDeltaT = maxCo*requiredDeltaT;

        // Update print info
        physicsModel::printInfo() = bool
        (
            runTime.timeIndex() % infoFrequency() == 0
         || mag(runTime.value() - runTime.endTime().value()) < SMALL
        );

        physicsModel::printInfo() = false;

        if (time().timeIndex() == 1)
        {
            Info<< nl << "Setting deltaT = " << newDeltaT
                << ", maxCo = " << maxCo << endl;
        }

        runTime.setDeltaT(newDeltaT);
    }
}


bool linGeomTotalDispSolid::evolve()
{
    if (solutionAlg() == solutionAlgorithm::PETSC_SNES)
    {
        return evolveSnes();
    }
    // else if (solutionAlg() == solutionAlgorithm::IMPLICIT_COUPLED)
    // {
    //     // Not yet implmented, although coupledUnsLinGeomLinearElasticSolid
    //     // could be combined with PETSc to achieve this.. todo!
    //     return evolveImplicitCoupled();
    // }
    else if (solutionAlg() == solutionAlgorithm::IMPLICIT_SEGREGATED)
    {
        return evolveImplicitSegregated();
    }
    else if (solutionAlg() == solutionAlgorithm::EXPLICIT)
    {
        return evolveExplicit();
    }
    else
    {
        FatalErrorIn("bool vertexCentredLinGeomSolid::evolve()")
            << "Unrecognised solution algorithm. Available options are "
            // << solutionAlgorithmNames_.names() << endl;
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::PETSC_SNES
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::IMPLICIT_SEGREGATED
               ]
            << solidModel::solutionAlgorithmNames_
               [
                   solidModel::solutionAlgorithm::EXPLICIT
               ]
            << endl;
    }

    // Keep compiler happy
    return true;
}


tmp<vectorField> linGeomTotalDispSolid::residualMomentum
(
    const volVectorField& D
)
{
    // Prepare result
    tmp<vectorField> tresidual(new vectorField(D.size(), vector::zero));
    vectorField& residual = tresidual.ref();

    // Enforce the boundary conditions
    const_cast<volVectorField&>(D).correctBoundaryConditions();

    // Update gradient of displacement
    mechanical().grad(D, gradD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    // Unit normal vectors at the faces
    const surfaceVectorField n(mesh().Sf()/mesh().magSf());

    // Traction vectors at the faces
    surfaceVectorField traction(n & fvc::interpolate(sigma()));

    // Add stabilisation to the traction
    // We add this before enforcing the traction condition as the stabilisation
    // is set to zero on traction boundaries
    // To-do: add a stabilisation traction function to momentumStabilisation
    const scalar scaleFactor =
        readScalar(stabilisation().dict().lookup("scaleFactor"));
    const surfaceTensorField gradDf(fvc::interpolate(gradD()));
    traction += scaleFactor*impKf_*(fvc::snGrad(D) - (n & gradDf));

    // Enforce traction boundary conditions
    enforceTractionBoundaries(traction, D, n);

    // The residual vector is defined as
    // F = div(sigma) + rho*g
    //     - rho*d2dt2(D) - dampingCoeff*rho*ddt(D) + stabilisationTerm
    // where, here, we roll the stabilisationTerm into the div(sigma)
    residual =
        fvc::div(mesh().magSf()*traction)
      + rho()
       *(
            g() - fvc::d2dt2(D) - dampingCoeff()*fvc::ddt(D)
        );

    // Make residual extensive as fvc operators are intensive (per unit volume)
    residual *= mesh().V();

    // Add optional fvOptions, e.g. MMS body force
    // Note that "source()" is already multiplied by the volumes
    residual -= fvOptions()(ds_, const_cast<volVectorField&>(D))().source();

    return tresidual;
}


tmp<sparseMatrix> linGeomTotalDispSolid::JacobianMomentum
(
    const volVectorField& D
)
{
    // Count the number of non-zeros for a Laplacian discretisation
    // This equals the sum of one plus the number of internal faces for each,
    // which can be calculated as nCells + 2*nInternalFaces
    // Multiply by the blockSize since we will form the block matrix
    const int blockSize = solidModel::twoD() ? 2 : 3;
    const label numNonZeros =
        blockSize*returnReduce
        (
            mesh().nCells() + 2.0*mesh().nInternalFaces(), sumOp<label>()
        );

    // Calculate a segregated approximation of the Jacobian
    fvVectorMatrix approxJ
    (
        fvm::laplacian(impKf_, D, "laplacian(DD,D)")
      - rho()*fvm::d2dt2(D)
    );

    if (dampingCoeff().value() > SMALL)
    {
        approxJ -= dampingCoeff()*rho()*fvm::ddt(D);
    }

    // Optional: under-relaxation of the linear system
    approxJ.relax();

    // Convert fvMatrix matrix to sparseMatrix

    // Initialise matrix
    tmp<sparseMatrix> tmatrix(new sparseMatrix(numNonZeros));
    sparseMatrix& matrix = tmatrix.ref();

    // Insert the diagonal
    {
        const vectorField diag(approxJ.DD());
        forAll(diag, blockRowI)
        {
            const tensor coeff
            (
                diag[blockRowI][vector::X], 0, 0,
                0, diag[blockRowI][vector::Y], 0,
                0,  0, diag[blockRowI][vector::Z]
            );

            const label globalBlockRowI =
                foamPetscSnesHelper::globalCells().toGlobal(blockRowI);

            matrix(globalBlockRowI, globalBlockRowI) = coeff;
        }
    }

    // Insert the off-diagonal
    {
        const labelUList& own = mesh().owner();
        const labelUList& nei = mesh().neighbour();
        const scalarField& upper = approxJ.upper();
        forAll(own, faceI)
        {
            const tensor coeff(upper[faceI]*I);

            const label blockRowI = own[faceI];
            const label blockColI = nei[faceI];

            const label globalBlockRowI =
                foamPetscSnesHelper::globalCells().toGlobal(blockRowI);
            const label globalBlockColI =
                foamPetscSnesHelper::globalCells().toGlobal(blockColI);

            matrix(globalBlockRowI, globalBlockColI) = coeff;
            matrix(globalBlockColI, globalBlockRowI) = coeff;
        }
    }

    // Collect the global cell indices from neighbours at processor boundaries
    // These are used to insert the off-processor coefficients
    // First, send the data
    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pD = D.boundaryField()[patchI];
        if (pD.type() == "processor")
        {
            // Take a copy of the faceCells (local IDs) and convert them to
            // global IDs
            labelList globalFaceCells(mesh().boundary()[patchI].faceCells());
            foamPetscSnesHelper::globalCells().inplaceToGlobal(globalFaceCells);

            // Send global IDs to the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh().boundary()[patchI]);
            procPatch.send
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );
        }
    }
    // Next, receive the data
    PtrList<labelList> neiProcGlobalIDs(D.boundaryField().size());
    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pD = D.boundaryField()[patchI];
        if (pD.type() == "processor")
        {
            neiProcGlobalIDs.set(patchI, new labelList(pD.size()));
            labelList& globalFaceCells = neiProcGlobalIDs[patchI];

            // Receive global IDs from the neighbour proc
            const processorFvPatch& procPatch =
                refCast<const processorFvPatch>(mesh().boundary()[patchI]);
            procPatch.receive
            (
                Pstream::commsTypes::blocking, globalFaceCells
            );
        }
    }

    // Insert the off-processor coefficients
    forAll(D.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pD = D.boundaryField()[patchI];

        if (pD.type() == "processor")
        {
            const vectorField& intCoeffs = approxJ.internalCoeffs()[patchI];
            const vectorField& neiCoeffs = approxJ.boundaryCoeffs()[patchI];
            const unallocLabelList& faceCells =
                mesh().boundary()[patchI].faceCells();
            const labelList& neiGlobalFaceCells = neiProcGlobalIDs[patchI];

            forAll(pD, faceI)
            {
                const label globalBlockRowI =
                    foamPetscSnesHelper::globalCells().toGlobal
                    (
                        faceCells[faceI]
                    );

                // On-proc diagonal coefficient
                {
                    const tensor coeff
                    (
                        intCoeffs[faceI][vector::X], 0, 0,
                        0, intCoeffs[faceI][vector::Y], 0,
                        0, 0, intCoeffs[faceI][vector::Z]
                    );

                    matrix(globalBlockRowI, globalBlockRowI) += coeff;
                }

                // Off-proc off-diagonal coefficient
                {
                    // Take care: we need to flip the sign
                    const tensor coeff
                    (
                        -neiCoeffs[faceI][vector::X], 0, 0,
                        0, -neiCoeffs[faceI][vector::Y], 0,
                        0, 0, -neiCoeffs[faceI][vector::Z]
                    );

                    const label globalBlockColI = neiGlobalFaceCells[faceI];

                    matrix(globalBlockRowI, globalBlockColI) += coeff;
                }
            }
        }
        else if (pD.coupled()) // coupled but not a processor boundary
        {
            FatalErrorIn
            (
                "tmp<sparseMatrix> linGeomTotalDispSolid::JacobianMomentum"
            )   << "Coupled boundaries (except processors) not implemented"
                << abort(FatalError);
        }
        // else non-coupled boundary contributions have already been added to
        // the diagonal
    }

    return tmatrix;
}


tmp<vectorField> linGeomTotalDispSolid::tractionBoundarySnGrad
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
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField n(patch.nf());

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - impK*pGradD))
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
