/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coupledTotalLagNewtonRaphsonBeam.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "HermiteSpline.H"
#include "spinTensor.H"
#include "pseudoVector.H"

#include "mergePoints.H"
#include "scalarMatrices.H"
#include "denseMatrixHelperFunctions.H"
#include "BlockEigenSolverOF.H"
#include "samplingFluid.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coupledTotalLagNewtonRaphsonBeam::evolve()
{
    beamModel::evolve();

    const dictionary& beamDict = beamProperties();

    const int nCorr
    (
        beamDict.lookupOrDefault<int>("nCorrectors", 1000)
    );

    // Tolerance for equation (force) residuals (R = Ax - B)
    scalar residualTol = beamDict.lookupOrDefault<scalar>("residualTol", -1);

    if (residualTol < 0)
    {
        if (beamDict.found("convergenceTol"))
        {
            residualTol = readScalar(beamDict.lookup("convergenceTol"));
        }
        else
        {
            residualTol = 1e-8;
        }
    }

    // Tolerance - displacement (DW) & rotation (DTheta_) correction variables
    const scalar solutionTol
    (
        beamDict.lookupOrDefault<scalar>("solutionTol", 1e-10)
    );

    // Absolute tolerance
    const scalar absoluteTol
    (
        beamDict.lookupOrDefault<scalar>("absoluteTol", 1e-30)
    );

    // Tolerance to check of the solver has diverged
    const scalar divTol
    (
        beamDict.lookupOrDefault<scalar>("divergenceTol", 1e6)
    );

    const label writeResidualFrequency
    (
        beamDict.lookupOrDefault<scalar>("infoFrequency", 1)
    );

    scalar initialResidualNorm = 1;
    scalar currentResidualNorm = GREAT;
    scalar deltaXNorm = GREAT;
    scalar XNorm = GREAT;

    // Let momentum contributions (e.g. almDrag) do any once-per-evolve work
    // such as external fluid sampling.  The result is frozen for the
    // duration of the inner NR loop below to avoid the 2-cycle that arose
    // when sampling was repeated each NR iteration and the sample point
    // drifted across fluid cell boundaries.
    forAll(momentumContribPtr_, i)
    {
        momentumContribPtr_[i].preEvolve(*this);
    }

    iOuterCorr() = 0;
    do
    {
        if (debug)
        {
            Info<< "iOuterCorr: " << iOuterCorr() << endl;
        }

        scalar tStart = runTime().elapsedCpuTime();

        {
            // Initialise the block system
            Field<scalarSquareMatrix> d
            (
                mesh().nCells(), scalarSquareMatrix(6, 0.0)
            );

            // Grab off-diagonal and set it to zero
            Field<scalarSquareMatrix> u
            (
                mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0)
            );
            // tensor6Field& u = WThetaEqn.upper().asSquare();
            // u = tensor6::zero;

            // Grab off-diagonal and set it to zero
            Field<scalarSquareMatrix> l
            (
                mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0)
            );
            // tensor6Field& l = WThetaEqn.lower().asSquare();
            // l = tensor6::zero;

            // Grab source and set it to zero
            Field<scalarRectangularMatrix> source
            (
                mesh().nCells(), scalarRectangularMatrix(6, 1, 0.0)
            );
            // vector6Field& source = WThetaEqn.source();
            // source = vector6::zero;

            scalarField deltaf = 1.0/mesh().deltaCoeffs().internalField();
            // const scalarField& wf = mesh().weights().internalField();

            const labelList& own = mesh().owner(); // unallocLabelList => labelList (ESI)
            const labelList& nei = mesh().neighbour();

            W_.boundaryFieldRef().updateCoeffs();
            Theta_.boundaryFieldRef().updateCoeffs();

            const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

            // Update the coefficients of W_ and Theta_ equations
            // These coefficients are inputs to the Jacobian matrix
            updateEqnCoefficients();

            // SB: Initial accleration and velocity values at 0th iteration
            // Valid for Newmark-beta integration scheme
            if (d2dt2SchemeName_ == "Newmark" && iOuterCorr() == 0)
            {
                // if (iOuterCorr() == 0)
                // {
                Accl_ = -(1/(runTime().deltaT()*betaN_))*U_.oldTime()
                  - (0.5/betaN_ - 1)*Accl_.oldTime();

                U_ = U_.oldTime()
                  + runTime().deltaT()*((1 - gammaN_)*Accl_.oldTime() + gammaN_*Accl_);

                dotOmega_ =
                  - (1/(runTime().deltaT()*betaN_))*Omega_.oldTime()
                  - (0.5/betaN_ - 1)*dotOmega_.oldTime();

                Omega_ =
                    Omega_.oldTime()
                  + runTime().deltaT()
                  *((1 - gammaN_)*dotOmega_.oldTime() + gammaN_*dotOmega_);
                // }
            }

            //- Assembling the diagonal and off-diagonal contributions
            //- of DW_ and DTheta_ to solve in block-coupled
            assembleMatrixCoefficients(d, l, u, source);

            // Add distributed force
            forAll(source, cellI)
            {
                source[cellI](0,0) -= q()[cellI].x()*L()[cellI];
                source[cellI](1,0) -= q()[cellI].y()*L()[cellI];
                source[cellI](2,0) -= q()[cellI].z()*L()[cellI];
            }

            // Add self-weight of the beam
            forAll(source, cellI)
            {
                // Beam index from the local cellZone: parallel-safe, unlike
                // globalCellIndex() whose localToGlobalCellAddressing_ is
                // never populated (segfaults in parallel runs)
                label bI = mesh().cellZones().whichZone(cellI);
                if (bI < 0)
                {
                    bI = 0;
                }
                // Effective density accounts for buoyancy when beam is submerged.
                // TODO: use rhoFluid(bI) if beams can be in different fluids.
                const scalar rhoEff = rho(bI).value() - rhoFluid().value();
                const scalar beamWeight = rhoEff*L()[cellI]*A(bI).value();

                source[cellI](0,0) -= beamWeight*g().component(0).value();
                source[cellI](1,0) -= beamWeight*g().component(1).value();
                source[cellI](2,0) -= beamWeight*g().component(2).value();
            }

            // Add point forces
            forAll(pointForces(), pfI)
            {
                // Get beam relative coordinates
                const label cellI = pointForces()[pfI].first().first();
                const scalar zeta = pointForces()[pfI].first().second();

                const vector F0 = pointForces()[pfI].second()(runTime().value());

                source[cellI](0,0) -= F0.x();
                source[cellI](1,0) -= F0.y();
                source[cellI](2,0) -= F0.z();

                const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

                const surfaceScalarField& dc = mesh().deltaCoeffs();

                vector DR = vector::zero;

                if (zeta > SMALL)
                {
                    const label faceID = own.find(cellI);
                    if (faceID == -1) // last cell
                    {
                        const labelList& faceCells =
                            mesh().boundary()[endPatchIndex()].faceCells();

                        const label bFaceID = faceCells.find(cellI);

                        DR = zeta*dRdS.boundaryField()[endPatchIndex()][bFaceID]
                            /dc.boundaryField()[endPatchIndex()][bFaceID];
                    }
                    else
                    {
                        DR = 0.5*zeta*dRdS.internalField()[faceID]*deltaf[faceID];
                    }
                }
                else if (zeta < -SMALL)
                {
                    const label faceID = nei.find(cellI);
                    if (faceID == -1) // first cell
                    {
                        const labelList& faceCells =
                            mesh().boundary()[startPatchIndex()].faceCells();

                        const label bFaceID = faceCells.find(cellI);

                        DR = zeta*dRdS.boundaryField()[startPatchIndex()][bFaceID]
                            /dc.boundaryField()[startPatchIndex()][bFaceID];
                    }
                    else
                    {
                        DR = 0.5*zeta*dRdS.internalField()[faceID]*deltaf[faceID];
                    }
                }

                const vector M0 = (spinTensor(DR) & F0);

                source[cellI](3,0) -= M0.x();
                source[cellI](4,0) -= M0.y();
                source[cellI](5,0) -= M0.z();
            }

            // Add distributed moment
            forAll(source, cellI)
            {
                source[cellI](3,0) -= m()[cellI].x()*L()[cellI];
                source[cellI](4,0) -= m()[cellI].y()*L()[cellI];
                source[cellI](5,0) -= m()[cellI].z()*L()[cellI];
            }
            // Add body force due to gravity and buoyancy if fluid is present
            // Note that for buoyant force calculation it is assumed
            // that the entire beam is submerged in the fluid
            // forAll(source, cellI)
            // {
            //     source[cellI](0,0) -=
            //         (
            //             rho().value() - rhoFluid().value()
            //         )*L()[cellI]*A().value()*g().component(0).value();

            //     source[cellI](1,0) -=
            //         (
            //             rho().value() - rhoFluid().value()
            //         )*L()[cellI]*A().value()*g().component(1).value();

            //     source[cellI](2,0) -=
            //         (
            //             rho().value() - rhoFluid().value()
            //         )*L()[cellI]*A().value()*g().component(2).value();
            // }

            // Add inertial components if not steadyState
            // SB Note: The inertia terms have two components
            // 1. Inertia term in the force balance equation
            // 2. Inertia term in the moment balance equation

            // The difference between Newmark and Euler is how
            // linear and angular velocities/accelerations are updated.
            // Note: The updations are done after solving the eqns
            // Also the implicit contributions of inertia terms that
            // go into the diagonal of the matrix are different
            // for Euler and Newmark.

            // Temporal Variables
            // U_ = Linear Velocity
            // Accl_ = Linear Acceleration
            // Omega_ = Angular Velocity
            // dotOmega_ = Angular Acceleration

            if (d2dt2SchemeName_ != "steadyState")
            {
                // Throw error if rho is set to zero for transient case
                for (label bI = 0 ; bI < nBeams() ; ++bI)
                {
                    if (rho(bI).value() == 0.0)
                    {
                        FatalErrorInFunction
                            << "The time integration scheme provided is " << d2dt2SchemeName_
                            << " but density 'rho' is not specified!!\n"
                            << "Specify density as scalar either locally in crossSectionModel dict "
                            << "of constant/beamProperties or as a global variable with dimensions "
                            << "inside coupledTotalLagNewtonRaphsonBeamCoeffs sub-dictionary"
                            << abort(FatalError);
                    }
                }

                // The EXPLICIT inertial contributions to source
                // 1. Initialise explicit inertial force contribution to zero
                // vectorField QRho(W_.size(), vector(0,0,0));
                volVectorField QRho
                (
                    IOobject
                    (
                        "QRho",
                        runTime().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedVector("zero", dimForce, vector::zero)
                );

                // 2. Initialise explicit inertial moment contribution to zero
                volVectorField MRho
                (
                    IOobject
                    (
                        "MRhoCoeff",
                        runTime().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedVector("zero", dimForce*dimLength, vector::zero)
                );

                if (d2dt2SchemeName_ == "Newmark")
                {
                    // 1. Inertial force
                    QRho = ARho_*L()*Accl_;

                    // 2. Inertial angular momentum
                    MRho =
                    (
                        L()
                        *(
                            (Lambda_ & (CIRho_ & dotOmega_))
                            + (Lambda_ & (Omega_ ^ (CIRho_ & Omega_)))
                        )
                    );

                }
                else if (d2dt2SchemeName_ == "Euler")
                {
                    // 1. Inertial force
                    QRho = ARho_*L()*fvc::ddt(U_);

                    // 2. Inertial angular momentum
                    MRho =
                    (
                        L()
                        *(
                          (Lambda_ & (CIRho_ & fvc::ddt(Omega_)))
                            + (Lambda_ & (Omega_ ^ (CIRho_ & Omega_)))
                        )
                    );

                }
                else
                {
                    FatalErrorInFunction
                        << "d2dt2SchemeName undefined: " << d2dt2SchemeName_
                        << exit(FatalError);
                }

                forAll(source, cellI)
                {
                    // 1. Add the explicit inertial force contribution
                    source[cellI](0,0) += QRho[cellI].x();
                    source[cellI](1,0) += QRho[cellI].y();
                    source[cellI](2,0) += QRho[cellI].z();

                    // 2. Add the explicit inertial moment contribution
                    source[cellI](3,0) += MRho[cellI].x();
                    source[cellI](4,0) += MRho[cellI].y();
                    source[cellI](5,0) += MRho[cellI].z();
                }

                // The IMPLICIT inertial contributions to source
                // 1. Initialise implicit inertial force contribution to zero
                // scalarField QRhoCoeff(W_.size(), 0.0);
                volScalarField QRhoCoeff
                (
                    IOobject
                    (
                        "QRhoCoeff",
                        runTime().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedScalar
                    (
                        "zero",
                        L().dimensions()*ARho_.dimensions(),
                        0
                    )
                );

                // 2. Initialise implicit inertial moment contribution to zero
                volTensorField MRhoCoeff
                (
                    IOobject
                    (
                        "MRhoCoeff",
                        runTime().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh(),
                    dimensionedTensor("zero", dimForce*dimLength, tensor::zero)
                );

                // if (newmark_)
                if (d2dt2SchemeName_ == "Newmark")
                {
                    // 1. Implicit inertial force coefficient
                    QRhoCoeff =
                        -L()*ARho_/(sqr(runTime().deltaT().value())*betaN_);

                    // 2. Implicit contribution of inertial torque because of
                    // Newton-Raphson linearisation
                    MRhoCoeff =
                        L()
                        *(
                            spinTensor
                            (
                                (Lambda_ & (Omega_ ^ (CIRho_ & Omega_)))
                                + (Lambda_ & (CIRho_ & dotOmega_))
                            )
                            + (
                                (spinTensor(Lambda_ & (CIRho_ & Omega_)))
                                - (
                                    Lambda_
                                    & (
                                        spinTensor(Omega_)
                                        & (CIRho_ &  Lambda_.T())
                                    )
                                )
                            )*(gammaN_/(betaN_*runTime().deltaT()))
                            - (
                                Lambda_ & (CIRho_ & Lambda_.T())
                            )*(1/(betaN_*sqr(runTime().deltaT())))
                        );
                }
                else if (d2dt2SchemeName_ == "Euler")
                {
                    // 1. Implicit inertial force coefficient
                    QRhoCoeff =
                        -L()*ARho_/sqr(runTime().deltaT().value());

                    // 2. ZT code: Implicit contribution of inertial torque
                    // because of Newton-Raphson linearisation
                    MRhoCoeff =
                        L()
                        *(
                            spinTensor(Lambda_ & (CIRho_ & Omega_))/runTime().deltaT()

                            - (Lambda_ & (CIRho_ & Lambda_.T()))/sqr(runTime().deltaT())

                            - spinTensor(Lambda_ & (CIRho_ & Omega_.oldTime()))/runTime().deltaT() // + sign ?

                            - spinTensor(Lambda_ & (spinTensor(Omega_) & (CIRho_ & Omega_)))

                            - spinTensor(Lambda_ & (CIRho_ & Omega_))/runTime().deltaT() // - sign

                            + (Lambda_ & (spinTensor(Omega_) & (CIRho_ & Lambda_.T())))/runTime().deltaT()
                        );
                }

                // Adding the implicit contribution to the diagonal
                forAll(d, cellI)
                    {
                        // 1. Add implicit inertia force contribution to diag
                        d[cellI](0,0) += QRhoCoeff[cellI];
                        d[cellI](1,1) += QRhoCoeff[cellI];
                        d[cellI](2,2) += QRhoCoeff[cellI];

                        // 2. Add implicit inertia momentum contribution to diag
                        d[cellI](3,3) += MRhoCoeff[cellI].xx();
                        d[cellI](3,4) += MRhoCoeff[cellI].xy();
                        d[cellI](3,5) += MRhoCoeff[cellI].xz();

                        d[cellI](4,3) += MRhoCoeff[cellI].yx();
                        d[cellI](4,4) += MRhoCoeff[cellI].yy();
                        d[cellI](4,5) += MRhoCoeff[cellI].yz();

                        d[cellI](5,3) += MRhoCoeff[cellI].zx();
                        d[cellI](5,4) += MRhoCoeff[cellI].zy();
                        d[cellI](5,5) += MRhoCoeff[cellI].zz();
                    }
            }

            // // Add run-time selectable momentum contributions
            if (momentumContribPtr_.size() > 0)
            {
                forAll(momentumContribPtr_, i)
                {
                    const vectorField linMomSource
                    (
                        momentumContribPtr_[i].linearMomentumSource(*this, U_, Accl_)
                    );

                    const vectorField angMomSource
                    (
                        momentumContribPtr_[i].angularMomentumSource(*this, U_, Accl_)
                    );

                    forAll(source, cellI)
                    {
                        source[cellI](0,0) += linMomSource[cellI][vector::X];
                        source[cellI](1,0) += linMomSource[cellI][vector::Y];
                        source[cellI](2,0) += linMomSource[cellI][vector::Z];

                        source[cellI](3,0) += angMomSource[cellI][vector::X];
                        source[cellI](4,0) += angMomSource[cellI][vector::Y];
                        source[cellI](5,0) += angMomSource[cellI][vector::Z];
                    }

                    // Jacobian diagonal contribution
                    const Field<scalarSquareMatrix> diagCoeff
                    (
                         momentumContribPtr_[i].diagCoeff(*this, U_, Accl_)
                    );

                    d += diagCoeff;
                }
            }
            // Block coupled solver call

            // Create Eigen linear solver
            BlockEigenSolverOF eigenSolver(d, l, u, own, nei);

            // Create solution vector
            Field<scalarRectangularMatrix> solVec
            (
                mesh().nCells(), scalarRectangularMatrix(6, 1, 0.0)
            );

            // Solve the linear system
            // currentResidualNorm is the imbalance vector
            currentResidualNorm = eigenSolver.solve(solVec, source); // peak RAM
            reduce(currentResidualNorm, sumOp<scalar>());

            if (iOuterCorr() == 0)
            {
                initialResidualNorm = currentResidualNorm;
            }

            // Copy the solution from solVec into the DW_ and DTheta_ fields
            vectorField& DWI = DW_;
            vectorField& DThetaI = DTheta_;

            forAll(solVec, cellI)
            {
                DWI[cellI][vector::X] = solVec[cellI](0, 0);
                DWI[cellI][vector::Y] = solVec[cellI](1, 0);
                DWI[cellI][vector::Z] = solVec[cellI](2, 0);

                DThetaI[cellI][vector::X] = solVec[cellI](3, 0);
                DThetaI[cellI][vector::Y] = solVec[cellI](4, 0);
                DThetaI[cellI][vector::Z] = solVec[cellI](5, 0);
            }

            DW_.correctBoundaryConditions();
            DTheta_.correctBoundaryConditions();

            // Update all the solution and output variables for
            // the current (Newton) iteration loop
            updateSolutionVariables();

            // Calculating norm of primary correction variables DW_ & DTheta_
            deltaXNorm = sqrt(sum(magSqr(solVec)));

            // Calculating the norm of W_ and Theta_
            XNorm =
                sqrt
                (
                    sum(magSqr(W_.primitiveFieldRef()))
                  + sum(magSqr(Theta_.primitiveFieldRef()))
                );
            reduce(deltaXNorm, sumOp<scalar>());
            reduce(XNorm, sumOp<scalar>());

            if (iOuterCorr() == 0)
            {
                initialResidualNorm = currentResidualNorm;
            }
        }

        scalar tEnd = runTime().elapsedCpuTime();

        totalSolutionTime_ += tEnd-tStart;

        if (debug)
        {
            Pout<< "Current total solution update time: "
                << totalSolutionTime_ << endl;
        }
    }
    while
    (
        !checkConvergence
        (
            currentResidualNorm,
            initialResidualNorm,
            deltaXNorm,
            XNorm,
            iOuterCorr()++,
            nCorr,
            residualTol,
            absoluteTol,
            solutionTol,
            divTol,
            writeResidualFrequency
        )
    );

    totalIter_ += iOuterCorr();

    if (debug)
    {
        Info<< "Total iterations " << totalIter_ << endl;
    }

    return initialResidualNorm;
}

//- Update the matrix coefficients for the Jacobian matrix
void coupledTotalLagNewtonRaphsonBeam::updateEqnCoefficients()
{
    const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

    // Total rotation matrix
    const surfaceTensorField Lambdaf((Lambdaf_ & refLambdaf_));

    CQW_ = (Lambdaf & (CQ_ & Lambdaf.T()));

    explicitQ_ = (Lambdaf & (CQ_ & (Gamma_)));

    explicitM_ = (Lambdaf & (CM_ & (K_)));

    CQTheta_ =
    (
        (
            Lambdaf & (CQ_ & Lambdaf.T())
        )
      & spinTensor(dRdS)
    )
    - spinTensor(Q_);
    // - spinTensor(explicitQ_);

    CMTheta_ = (Lambdaf & (CM_ & Lambdaf.T()));

    // CMTheta2_ = -spinTensor(Lambdaf & (CM_ & (K_ - KP_)))

    CMTheta2_ = -spinTensor(M_);
    // CMTheta2_ = -spinTensor(explicitM_);

    CMQW_ =
        0.5
       *(
            (
                spinTensor(dRdS)
              & (Lambdaf & (CQ_ & Lambdaf.T()) )
            )
          - spinTensor(Q_)
        //   - spinTensor(explicitQ_)
        )/mesh().deltaCoeffs();
        // + (Lambdaf & (CDMDGamma_ & Lambdaf.T())); // Check for Kirchhoff beam

    CMQTheta_ =
        0.5
       *(
            (
                (
                    spinTensor(dRdS)
                  & (Lambdaf & (CQ_ & Lambdaf.T()))
                )
              & spinTensor(dRdS)
            )
            - (
                spinTensor(dRdS)
              & (
                    spinTensor(Q_)
                    // spinTensor(explicitQ_)
                )
            )
        )/mesh().deltaCoeffs();

    explicitMQ_ =
        0.5
       *(
            spinTensor(dRdS) & explicitQ_
        )/mesh().deltaCoeffs();

    // Correct at boundary
    forAll(CMQW_.boundaryFieldRef(), patchI)
    {
        if (!CMQW_.boundaryFieldRef()[patchI].coupled())
        {
            CMQW_.boundaryFieldRef()[patchI] *= 2;
            CMQTheta_.boundaryFieldRef()[patchI] *= 2;
            explicitMQ_.boundaryFieldRef()[patchI] *= 2;
        }
    }
}

//- Update the solution variables like displacements, rotations,
//- translational and rotational strains, forces and moments
void coupledTotalLagNewtonRaphsonBeam::updateSolutionVariables()
{
    Theta_.primitiveFieldRef() +=
        (Lambda_.internalField().T() & DTheta_.internalField()); // Get back to this

    // Theta_.internalField().replace(0, 0);

    Theta_.correctBoundaryConditions();
    Theta_.storePrevIter();

    //WThetaEqn.retrieveSolution(0, DW_.internalField());
    //DW_.boundaryField().evaluateCoupled();
    W_.primitiveFieldRef() += DW_.internalField();

    W_.correctBoundaryConditions();
    W_.storePrevIter();

    // Store U for optional field relaxation after updating it below
    U_.storePrevIter();

    // Update mean line linear velocity and acceleration fields
    if (d2dt2SchemeName_ == "steadyState")
    {
        // Do Nothing, no update of temporal variables reqd
    }
    else if (d2dt2SchemeName_ == "Newmark")
    {
        // SB: Update mean line acceleration field
        U_ += (1/runTime().deltaT())*(gammaN_/betaN_)*DW_;

        // SB: Update mean line acceleration field
        Accl_ += (1/(sqr(runTime().deltaT())*betaN_))*DW_;

    }
    else if (d2dt2SchemeName_ == "Euler")
    {
        U_ = fvc::ddt(W_);
        Accl_ = fvc::ddt(U_);
    }
    else
    {
        FatalErrorInFunction
            << "Provided d2dt2Scheme is not implemented!"
            << "Valid choices are steadyState, Euler, Newmark"
            << abort(FatalError);
    }

    // Apply U relaxation if configured in fvSolution
    U_.relax();

    const surfaceVectorField DThetaf(fvc::interpolate(DTheta_));

    const surfaceScalarField magDThetaf(mag(DThetaf) + SMALL);
    const surfaceTensorField DThetaHat(spinTensor(DThetaf));

    const dimensionedTensor I("I", dimless, tensor::I);

    // Tangent operator
    const surfaceTensorField DT
    (
        (Foam::sin(magDThetaf)/magDThetaf)*I
        + (
            (1.0-Foam::sin(magDThetaf)/magDThetaf)/sqr(magDThetaf)
        )
        *(DThetaf*DThetaf)
        + (
            (1.0-Foam::cos(magDThetaf))/sqr(magDThetaf)
        )*DThetaHat
    );

    // Update bending strain vector
    K_ +=
    (
        (refLambdaf_.T() & Lambdaf_.T())
        // & fvc::snGrad(DTheta_)
        & (DT.T() & fvc::snGrad(DTheta_))
    );

    // Rodrigues formula
    const surfaceTensorField DLambdaf(rotationMatrix(DThetaf));

    // Update rotation matrix
    Lambdaf_ = (DLambdaf & Lambdaf_);

    // Update cell-centre rotation matrix
    const volTensorField DLambda(rotationMatrix(DTheta_));

    Lambda_ = (DLambda & Lambda_);
    // interpolateRotationMatrix(*this, Lambdaf_, Lambda_);

    // Update mean line angular velocity and acceleration fields
    if (d2dt2SchemeName_ == "steadyState")
    {
        // Do nothing; no update of temporal variables reqd
    }
    else if (d2dt2SchemeName_ == "Newmark")
    {
        // Update angular velocity
        // without tangent space
        Omega_ +=
            (
                gammaN_/(betaN_*runTime().deltaT())
            )*(Lambda_.T() & DTheta_);

        // Update angular acceleration
        dotOmega_ +=
            (
                1/(betaN_*sqr(runTime().deltaT()))
            )*(Lambda_.T() & DTheta_);
    }
    else if (d2dt2SchemeName_ == "Euler") // First order Euler scheme
    {
        Omega_ = axialVector(Lambda_.T() & fvc::ddt(Lambda_));

        // dotOmega_ = fvc::ddt(Omega_);
    }
    else
    {
        FatalErrorInFunction
            << "Provided d2dt2Scheme is not implemented!"
            << "Valid choices are steadyState, Euler, Newmark"
            << abort(FatalError);
    }

    // }

    // Update axial and shear strain vector
    {
        // W_.correctBoundaryConditions();

        surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

        Gamma_ = (refLambdaf_.T() & ((Lambdaf_.T() & dRdS) - dR0Ds_));
    }

    // Q_ = (CQ_ & (Gamma_ - GammaP_));
    // M_ = (CM_ & (K_ - KP_));

    // Q_ = ((Lambdaf_ & refLambdaf_) & (CQ_ & (Gamma_ - GammaP_)));
    // M_ = ((Lambdaf_ & refLambdaf_) & (CM_ & (K_ - KP_)));

    // Calculate Q, where we ignore the orientation check
    {
        explicitQ_.setOriented(true);
        surfaceVectorField implicitQ(CQW_ & fvc::snGrad(DW_));
        implicitQ.setOriented(false);
        implicitQ += (CQTheta_ & fvc::interpolate(DTheta_));
        implicitQ.setOriented(true);
        Q_ = explicitQ_ + implicitQ;
    }

    // Calculate M, where we ignore the orientation check
    {
        explicitM_.setOriented(true);
        surfaceVectorField implicitM(CMTheta_ & fvc::snGrad(DTheta_));
        implicitM.setOriented(false);
        implicitM += (CMTheta2_ & fvc::interpolate(DTheta_));
        implicitM.setOriented(true);
        M_ = explicitM_ + implicitM;
    }
}

//- Function to check convergence of the outer (Newton) loop
bool coupledTotalLagNewtonRaphsonBeam::checkConvergence
(
    const scalar currentResidualNorm,
    const scalar initialResidualNorm,
    const scalar deltaXNorm,
    const scalar xNorm,
    const label iteration,
    const label maxIterations,
    const scalar rtol,
    const scalar atol,
    const scalar stol,
    const scalar divtol,
    const label writeResidualFrequency,
    const bool writeConvergedReason
)
{
    // Precompute tolerances
    const scalar relativeResidualTol = rtol*initialResidualNorm;
    const scalar stepTolerance = stol*xNorm;
    // Log residuals if enabled
    if (writeResidualFrequency > 0)
    {
        if (iteration == 0)
        {
            // Print the header with fixed widths
            Info<< setw(10) << "Iteration"
                << setw(20) << "Residual Norm"
                << setw(20) << "Step Norm" << endl;
        }

        // Print each iteration's data with aligned fields
        if (iteration % writeResidualFrequency == 0)
        {
            Info<< setw(10) << iteration
                << setw(20) << currentResidualNorm
                << setw(20) << deltaXNorm << endl;
        }
    }

    // 1. Check Residual Norm against absolute tolerance
    if (currentResidualNorm <= atol)
    {
        if (writeConvergedReason)
        {
            Info<< setw(10) << iteration
                << setw(20) << ": Converged - Absolute residual tolerance met."
                << endl;
        }
        return true;
    }

    // 2. Check Residual Norm Convergence
    if (currentResidualNorm <= relativeResidualTol)
    {
        if (writeConvergedReason)
        {
            Info<< setw(10) << iteration
                << setw(20) << ": Converged - Relative residual tolerance met."
                << endl;
        }
        return true;
    }

    // 3. Check Step Norm Convergence
    if (deltaXNorm <= stepTolerance)
    {
        if (writeConvergedReason)
        {
            Info<< setw(10) << iteration
                << setw(20) << ": Converged - Step norm relative tolerance met."
                << endl;
        }
        return true;
    }

    // 4. Check Divergence
    if (currentResidualNorm >= divtol*initialResidualNorm && iteration > 1)
    {
        FatalErrorInFunction
            << "Iteration " << iteration
            << setw(20) << ": Diverged - Residual grew excessively."
            << "\nYou can try setting the `divergenceTol` in `constant/beamProperties` "
            << "to a different value. The default value is: " << divtol
            << abort(FatalError);
        return false;
    }

    // 5. Check Maximum Iterations
    if (iteration >= maxIterations)
    {
        FatalErrorInFunction
            << "Iteration " << iteration
            << setw(20) << ": Failed - Maximum iterations reached."
            << abort(FatalError);
        return false;
    }

    // 6. Not Converged Yet
    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
