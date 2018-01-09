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

#include "solidModel.H"
#include "volFields.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "RectangularMatrix.H"
#include "solidTractionFvPatchVectorField.H"
#include "blockSolidTractionFvPatchVectorField.H"
#include "fvcGradf.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidModel, 0);
    defineRunTimeSelectionTable(solidModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidModel::checkWedges() const
{
    const fvMesh& mesh = this->mesh();

    label nWedgePatches = 0;
    vector wedgeDirVec = vector::zero;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());

            // Make sure that solidWedge is used instead of wedge
            if
            (
                DD_.boundaryField()[patchI].type() == "wedge"
             && D_.boundaryField()[patchI].type() == "wedge"
            )
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "solidWedge should be used on displacement solution "
                    << "field wedge patches as non-orthogonal corrections "
                    << "are important!"
                    << abort(FatalError);
            }
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        if (nWedgePatches != 2)
        {
            FatalErrorIn("void Foam::solidModel::checkWedges() const")
                << "For axisymmetric cases, there should be exactly two wedge "
                << "patches!" << abort(FatalError);
        }

        Info<< nl << "Axisymmetric case: disabling the solution in the "
            << "out-of-plane direction" << endl;

        // We will use const_cast to disable the out-of-lane direction
        Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;

                wordList dirs(3);
                dirs[0] = "x";
                dirs[1] = "y";
                dirs[2] = "z";
                Info<< "    out-of-plane direction: " << dirs[cmpt] << nl
                    << endl;
            }
            else
            {
                solD[cmpt] = 1;
            }
        }
    }


    // Check all the face normals are in the same direction on the wedge patches
    // This is to avoid the case where a wedge patch is composed of two
    // disconnected regions with one on the front and one on the back
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // Unit face normals on processor
            const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

            if (nf.size() == 0)
            {
                FatalErrorIn("void Foam::solidModel::checkWedges() const")
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                    << nl << "Every processor should have at least one face on "
                    << "each wedge patch"
                    << abort(FatalError);
            }

            // Check that all the wedge face normals point in the same direction

            vector firstFaceNOnMasterProc = vector::zero;

            if (Pstream::master())
            {
                firstFaceNOnMasterProc = nf[0];
            }

            // Sync in parallel so that all processors have the master vector
            reduce(firstFaceNOnMasterProc, sumOp<vector>());

            forAll(nf, faceI)
            {
                if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
                {
                    FatalErrorIn("void Foam::solidModel::checkWedges() const")
                        << "On wedge patch "
                        << mesh.boundaryMesh()[patchI].name() << " there are at "
                        << "least two faces with unit normals in the opposite "
                        << "directions" << nl
                        << "Please check that the wedge patches are correctly "
                        << "defined"
                        << abort(FatalError);
                }
            }
        }
    }
}


void Foam::solidModel::calcGlobalFaceZones() const
{
    // Find global face zones
    if (globalFaceZonesPtr_)
    {
        FatalErrorIn
        (
            "void solidModel::calcGlobalFaceZones() const"
        )
            << "Global face zones already found"
            << abort(FatalError);
    }

    if (Pstream::parRun())
    {
        SLList<label> globalFaceZonesSet;

        // Lookup globalFaceZones from decomposeParDict
        // For FSI cases, we need to look in a different location for the dict

        word decompDictName = "system/decomposeParDict";

        if
        (
            isDir
            (
                mesh().time().rootPath()/mesh().time().caseName()
                /"../system/solid"
            )
        )
        {
            decompDictName = "../system/solid/decomposeParDict";
        }

        Info<< "Reading decomposeParDict " << decompDictName << endl;

    Info<< "line " << __LINE__ << endl;
        IOdictionary decompDict
        (
            IOobject
            (
                decompDictName,
                mesh().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (decompDict.found("globalFaceZones"))
        {
            wordList globalFaceZoneNames(decompDict.lookup("globalFaceZones"));

            const faceZoneMesh& faceZones = mesh().faceZones();

            forAll(globalFaceZoneNames, nameI)
            {
                const label zoneID =
                    faceZones.findZoneID(globalFaceZoneNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn(type() + "::findGlobalFaceZones")
                        << "Cannot find globalFaceZone:"
                        << " " << globalFaceZoneNames[nameI]
                        << abort(FatalError);
                }

                globalFaceZonesSet.insert(zoneID);
            }

            globalFaceZonesPtr_ = new labelList(globalFaceZonesSet);
        }
        else
        {
            globalFaceZonesPtr_ = new labelList(0);
        }
    }
    else
    {
        globalFaceZonesPtr_ = new labelList(0);
    }
}


void Foam::solidModel::calcGlobalToLocalFaceZonePointMap() const
{
    // Find global face zones
    if (globalToLocalFaceZonePointMapPtr_)
    {
        FatalErrorIn
        (
            "void solidModel::calcGlobalToLocalFaceZonePointMap() const"
        )   << "Global to local face zones point map already exists"
            << abort(FatalError);
    }

    globalToLocalFaceZonePointMapPtr_ =
        new labelListList(globalFaceZones().size());

    labelListList& globalToLocalFaceZonePointMap =
        *globalToLocalFaceZonePointMapPtr_;

    const labelList& globalFaceZones = this->globalFaceZones();

    forAll(globalFaceZones, zoneI)
    {
        const label curZoneID = globalFaceZones[zoneI];

        Info<< "Creating faceMap for globalFaceZones "
            << mesh().faceZones()[curZoneID].name()<< endl;

        labelList curMap(mesh().faceZones()[curZoneID]().nPoints(), -1);

        vectorField fzGlobalPoints =
            mesh().faceZones()[curZoneID]().localPoints();

        // Set all slave points to zero because only the master order is used
        if(!Pstream::master())
        {
            fzGlobalPoints = vector::zero;
        }

        // Pass points to all procs
        reduce(fzGlobalPoints, sumOp<vectorField>());

        // Now every proc has the master's list of FZ points
        // every proc must now find the mapping from their local FZ points to
        // the global FZ points

        const vectorField& fzLocalPoints =
            mesh().faceZones()[curZoneID]().localPoints();

        const edgeList& fzLocalEdges =
            mesh().faceZones()[curZoneID]().edges();

        const labelListList& fzPointEdges =
            mesh().faceZones()[curZoneID]().pointEdges();

        scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

        forAll(minEdgeLength, pI)
        {
            const labelList& curPointEdges = fzPointEdges[pI];

            forAll(curPointEdges, eI)
            {
                const scalar Le =
                    fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);

                if (Le < minEdgeLength[pI])
                {
                    minEdgeLength[pI] = Le;
                }
            }
        }

        forAll(fzGlobalPoints, globalPointI)
        {
            boolList visited(fzLocalPoints.size(), false);

            forAll(fzLocalPoints, procPointI)
            {
                if (!visited[procPointI])
                {
                    visited[procPointI] = true;

                    label nextPoint = procPointI;

                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[nextPoint]
                          - fzGlobalPoints[globalPointI]
                        );

                    if (curDist < 1e-4*minEdgeLength[nextPoint])
                    {
                        curMap[globalPointI] = nextPoint;
                        break;
                    }

                    label found = false;

                    while (nextPoint != -1)
                    {
                        const labelList& nextPointEdges =
                            fzPointEdges[nextPoint];

                        scalar minDist = GREAT;
                        label index = -1;
                        forAll(nextPointEdges, edgeI)
                        {
                            label curNgbPoint =
                                fzLocalEdges[nextPointEdges[edgeI]]
                               .otherVertex(nextPoint);

                            if (!visited[curNgbPoint])
                            {
                                visited[curNgbPoint] = true;

                                scalar curDist =
                                    mag
                                    (
                                        fzLocalPoints[curNgbPoint]
                                      - fzGlobalPoints[globalPointI]
                                    );

                                if (curDist < 1e-4*minEdgeLength[curNgbPoint])
                                {
                                    curMap[globalPointI] = curNgbPoint;
                                    found = true;
                                    break;
                                }
                                else if (curDist < minDist)
                                {
                                    minDist = curDist;
                                    index = curNgbPoint;
                                }
                            }
                        }

                        nextPoint = index;
                    }

                    if (found)
                    {
                        break;
                    }
                }
            }
        }

        forAll(curMap, globalPointI)
        {
            if (curMap[globalPointI] == -1)
            {
                FatalErrorIn
                (
                    "solidModel::calcGlobalToLocalFaceZonePointMap()"
                )   << "local to global face zone point map is not correct"
                    << abort(FatalError);
            }
        }

        globalToLocalFaceZonePointMap[zoneI] = curMap;
    }
}


void Foam::solidModel::updateGlobalFaceZoneNewPoints
(
    const pointField& pointDDI,
    pointField& newPoints
)
{
    const labelList& globalFaceZones = this->globalFaceZones();
    const labelListList& globalToLocalFaceZonePointMap =
        this->globalToLocalFaceZonePointMap();

    forAll(globalFaceZones, zoneI)
    {
        const label curZoneID = globalFaceZones[zoneI];

        const labelList& curMap = globalToLocalFaceZonePointMap[zoneI];

        const labelList& curZoneMeshPoints =
            mesh().faceZones()[curZoneID]().meshPoints();

        vectorField curGlobalZonePointDispl
        (
            curZoneMeshPoints.size(),
            vector::zero
        );

        // Inter-proc points are shared by multiple procs
        // pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

        forAll(curGlobalZonePointDispl, globalPointI)
        {
            label localPoint = curMap[globalPointI];

            if(curZoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = curZoneMeshPoints[localPoint];

                curGlobalZonePointDispl[globalPointI] = pointDDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(curGlobalZonePointDispl, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            // Now average the displacement between all procs
            curGlobalZonePointDispl /= pointNumProcs;
        }

        // The curZonePointsDisplGlobal now contains the correct face zone
        // displacement in a global master processor order, now convert them
        // back into the local proc order

        vectorField curZonePointDispl(curZoneMeshPoints.size(), vector::zero);

        forAll(curGlobalZonePointDispl, globalPointI)
        {
            label localPoint = curMap[globalPointI];

            curZonePointDispl[localPoint] =
                curGlobalZonePointDispl[globalPointI];
        }

        forAll(curZonePointDispl, pointI)
        {
            // Unused points
            if (curZoneMeshPoints[pointI] >= mesh().nPoints())
            {
                newPoints[curZoneMeshPoints[pointI]] +=
                    curZonePointDispl[pointI];
            }
        }
    }
}


void Foam::solidModel::makeThermalModel() const
{
    if (!thermalPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeThermalModel() const")
            << "pointer already set!" << abort(FatalError);
    }

    thermalPtr_.set
    (
        new thermalModel(mesh())
    );
}


void Foam::solidModel::makeMechanicalModel() const
{
    if (!mechanicalPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeMechanicalModel() const")
            << "pointer already set!" << abort(FatalError);
    }

    mechanicalPtr_.set
    (
        new mechanicalModel(mesh(), nonLinGeom(), incremental())
    );
}


void Foam::solidModel::makeRho() const
{
    if (!rhoPtr_.empty())
    {
        FatalErrorIn("void Foam::solidModel::makeRho() const")
            << "pointer already set!" << abort(FatalError);
    }

    rhoPtr_.set
    (
        new volScalarField(mechanical().rho())
    );
}


const Foam::pointVectorField& Foam::solidModel::pointDorPointDD() const
{
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Updated Lagrangian approaches move the mesh at the end of each
        // time-step so we use the increment of displacement field to calculate
        // the current deformed face zone points
        return pointDD();
    }
    else
    {
        // As linearGeometry and total Lagrangian approaches do not move the
        // mesh, we use the total displacement field to calculate the current
        // deformed face zone points
        return pointD();
    }
}


// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

Foam::thermalModel& Foam::solidModel::thermal()
{
    if (thermalPtr_.empty())
    {
        makeThermalModel();
    }

    return thermalPtr_();
}


Foam::mechanicalModel& Foam::solidModel::mechanical()
{
    if (mechanicalPtr_.empty())
    {
        makeMechanicalModel();
    }

    return mechanicalPtr_();
}


const Foam::volScalarField& Foam::solidModel::rho() const
{
    if (rhoPtr_.empty())
    {
        makeRho();
    }

    return rhoPtr_();
}


void Foam::solidModel::relaxField(volVectorField& D, int iCorr)
{
    if (relaxationMethod_ == "fixed")
    {
        // Fixed under-relaxation
        D.relax();
    }
    else if (relaxationMethod_ == "Aitken")
    {
        // See Aitken method at:
        // http://empire-multiphysics.com/projects/empire/wiki/Aitken_Relaxation
        // and
        // A partitioned solution approach for electro-thermo-
        // problems, Patrick Erbts, Stefan Hartmann, Alexander Duster.

        // Store aitkenResidual previous iteration
        aitkenResidual_.storePrevIter();

        // Calculate new aitkenResidual
        aitkenResidual_ = D.prevIter() - D;

        if (iCorr == 0)
        {
            // Fixed under-relaxation is applied on the first iteration
            aitkenAlpha_ = 1.0;

            if (mesh().solutionDict().relaxField(D.name()))
            {
                aitkenAlpha_ =
                    mesh().solutionDict().fieldRelaxationFactor(D.name());
            }
        }
        else
        {
            const volVectorField aitkenResidualDelta =
                aitkenResidual_.prevIter() - aitkenResidual_;

            // Update the relaxation factor field
            aitkenAlpha_ =
                aitkenAlpha_*(aitkenResidual_.prevIter() & aitkenResidualDelta)
               /(
                    magSqr(aitkenResidualDelta)
                  + dimensionedScalar("SMALL", dimLength*dimLength, SMALL)
                );

            // Bound alpha between 0.0 and 2.0
            // This may not be necessary but it seems to help convergence
            aitkenAlpha_ = max(0.0, min(2.0, aitkenAlpha_));
        }

        // Relax the field
        D -= aitkenAlpha_*aitkenResidual_;
    }
    else if (relaxationMethod_ == "QuasiNewton")
    {
        // This method is a modified form of the IQNILS by Degroote et al.

        // J. Degroote, K.-J. Bathe and J. Vierendeels.
        // A fluid solid interaction solver with IQN-ILS coupling algorithm.
        // Performance of a new partitioned procedure versus a monolithic
        // procedure in fluid-solid interaction. Computers & Solids

        if (iCorr == 0 || iCorr % QuasiNewtonRestartFreq_ == 0)
        {
            // Clean up data from old time steps

            if (debug)
            {
                Info<< "Modes before clean-up : " << QuasiNewtonT_.size();
            }

            while (true)
            {
                if (QuasiNewtonT_.size())
                {
                    if
                    (
                        runTime().timeIndex() > QuasiNewtonT_[0]
                     || iCorr % QuasiNewtonRestartFreq_ == 0
                    )
                    {
                        for (label i = 0; i < QuasiNewtonT_.size() - 1; i++)
                        {
                            QuasiNewtonT_[i] = QuasiNewtonT_[i + 1];
                            QuasiNewtonV_[i] = QuasiNewtonV_[i + 1];
                            QuasiNewtonW_[i] = QuasiNewtonW_[i + 1];
                        }

                        QuasiNewtonT_.remove();
                        QuasiNewtonV_.remove();
                        QuasiNewtonW_.remove();
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }

            if (debug)
            {
                Info<< ", modes after clean-up : " << QuasiNewtonT_.size()
                    << endl;
            }
        }
        else if (iCorr == 1 || iCorr % QuasiNewtonRestartFreq_ == 1)
        {
            // Set reference in the first coupling iteration
            unrelaxedDRef_ = D;
            DRef_ = D.prevIter();
        }
        else
        {
            // Store the input vector field, defined as the previous iteration
            // D field (after relaxation) minus the Dp previous iteration field
            // in the first iteration (after relaxation)
            QuasiNewtonV_.append
            (
                (D - D.prevIter()) - (unrelaxedDRef_ - DRef_)
            );

            // V should be (from FSI paper):
            // DeltaR^{k-1} = R^{k-1} - R^k = DPrevIt.PrevIt - DPrevIter
            // DeltaR^{k-2} = R^{k-2} - R^k = D.PI.PI.PI - D.PI
            // ...
            // DeltaR^{0} =  R^0 - R^k = DRef - D.prevIter
            // Or in the general paper:
            // V_i = p_k - p_i    for i = 0, 1, ..., k - 1
            // V_{k-1} = p_k - p_{k-1} = D.PI - D.PI.PI
            // V_{k-2} = p_k - p_{k-2} = D.PI - D.PI.PI.PI
            // ...
            // V_{0} = p_k - p_{0} = D.PI - DRef
            // BUT, the implemented code does this:
            // V_i = p_{i+1} - p_0    for i = 0, 1, ..., k - 1
            // V_{k-1} = p_{k} - p_0 = D.PI - DRef
            // V_{k-2} = p_{k-1} - p_0 = D.PI{k-1} - DRef
            // ...
            // V_{0} = p_{1} - p_{0} = D.PI_1 - DRef
            // This means that we just append  the following line each
            // iteration:
            // V_{k-1} = p_{k} - p_0 = D.PI - DRef
            // It this equivalent?
            // We could try implementing it as described in the paper, but this
            // will require D and D.prevIter and their history

            // Store the output vector field, defined as the current iteration
            // D field (before relaxation) minus the D field in the first
            // iteration (before relaxation)
            QuasiNewtonW_.append(D - unrelaxedDRef_);

            // Store the time index
            QuasiNewtonT_.append(runTime().timeIndex());
        }

        if (QuasiNewtonT_.size() > 1)
        {
            // Consider QuasiNewtonV as a matrix V
            // with as columns the items
            // in the DynamicList and calculate the QR-decomposition of V
            // with modified Gram-Schmidt
            label cols = QuasiNewtonV_.size();
            RectangularMatrix<scalar> R(cols, cols, 0.0);
            RectangularMatrix<scalar> C(cols, 1);
            RectangularMatrix<scalar> Rcolsum(1, cols);
            // philipc: do need for dynamic list for Q
            //DynamicList<vectorField> Q(cols);
            List<vectorField> Q(cols);

            for (label i = 0; i < cols; i++)
            {
                //Q.append(QuasiNewtonV_[cols - 1 - i]);
                Q[i] = QuasiNewtonV_[cols - 1 - i];
            }

            for (label i = 0; i < cols; i++)
            {
                // Normalize column i
                R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
                Q[i] /= R[i][i];

                // Orthogonalize columns to the right of column i
                for (label j = i+1; j < cols; j++)
                {
                    R[i][j] = sum(Q[i] & Q[j]);
                    Q[j] -= R[i][j]*Q[i];
                }

                // Project minus the residual vector on the Q
                C[i][0] =
                    sum
                    (
                        Q[i]
                      & (
                          D.prevIter().internalField()
                        - D.internalField()
                      )
                    );
            }

            // Solve the upper triangular system
            for (label j = 0; j < cols; j++)
            {
                Rcolsum[0][j] = 0.0;
                for (label i = 0; i < (j + 1); i++)
                {
                    Rcolsum[0][j] += cmptMag(R[i][j]);
                }
            }
            scalar epsilon = 1.0E-10*max(Rcolsum);
            for (label i = 0; i < cols; i++)
            {
                if (cmptMag(R[i][i]) > epsilon)
                {
                    for (label j = i + 1; j < cols; j++)
                    {
                        R[i][j] /= R[i][i];
                    }
                    C[i][0] /= R[i][i];
                    R[i][i] = 1.0;
                }
            }
            for (label j = (cols - 1); j >= 0; j--)
            {
                if (cmptMag(R[j][j]) > epsilon)
                {
                    for (label i = 0; i < j; i++)
                    {
                        C[i][0] -= C[j][0]*R[i][j];
                    }
                }
                else
                {
                    C[j][0] = 0.0;
                }
            }

            // Update D
            for (label i = 0; i < cols; i++)
            {
                D.internalField() += QuasiNewtonW_[i]*C[cols - 1 - i][0];
            }

            D.correctBoundaryConditions();
        }
        else
        {
            // Fixed under-relaxation during startup
            D.relax();
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::solidModel::relaxField(volVectorField& D, int iCorr)"
        )   << "relaxationMethod '" << relaxationMethod_ << "' unknown!"
            << " Options are fixed, Aitken or QuasiNewton" << abort(FatalError);
    }
}


bool Foam::solidModel::converged
(
    const int iCorr,
    const scalar solverPerfInitRes,
    const int solverPerfNIters,
    const volVectorField& D
)
{
    // We will check three residuals:
    // - relative displacement residual
    // - linear equation residual
    // - material model residual
    bool converged = false;

    // Calculate displacement residual based on the relative change of D
    scalar denom = gMax(mag(D.internalField() - D.oldTime().internalField()));
    if (denom < SMALL)
    {
        denom = max(gMax(mag(D.internalField())), SMALL);
    }
    const scalar residualDD =
        gMax(mag(D.internalField() - D.prevIter().internalField()))/denom;

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at least 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol_)
    {
        if
        (
            solverPerfInitRes < solutionTol_
         && residualDD < solutionTol_
        )
        {
            Info<< "    Both residuals have converged" << endl;
            converged = true;
        }
        else if (residualDD < alternativeTol_)
        {
            Info<< "    The relative residual has converged" << endl;
            converged = true;
        }
        else if (solverPerfInitRes < alternativeTol_)
        {
            Info<< "    The solver residual has converged" << endl;
            converged = true;
        }
        else
        {
            converged = false;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res, relRes, matRes, iters" << endl;
    }
    else if (iCorr % infoFrequency_ == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfInitRes
            << ", " << residualDD
            << ", " << materialResidual
            << ", " << solverPerfNIters << endl;

        if (residualFilePtr_.valid())
        {
            residualFilePtr_()
                << solverPerfInitRes << " "
                << residualDD << " "
                << materialResidual
                << endl;
        }

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr_ - 1)
    {
        maxIterReached_++;
        Warning
            << "Max iterations reached within momentum loop" << endl;
    }

    return converged;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModel::solidModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime),
    IOdictionary
    (
        IOobject
        (
            "solidProperties",
            // If region == "region0" then read from the main case
            // Otherwise, read from the region/sub-mesh directory e.g.
            // constant/fluid or constant/solid
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),
    solidProperties_(subDict(type + "Coeffs")),
    thermalPtr_(NULL),
    mechanicalPtr_(NULL),
    D_
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    DD_
    (
        IOobject
        (
            "DD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    U_
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength/dimTime, vector::zero)
    ),
    pMesh_(mesh()),
    pointD_
    (
        IOobject
        (
            "pointD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pointDD_
    (
        IOobject
        (
            "pointDD",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    gradD_
    (
        IOobject
        (
            "grad(" + D_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    gradDD_
    (
        IOobject
        (
            "grad(" + DD_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    rhoPtr_(NULL),
    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    solutionTol_
    (
        solidProperties().lookupOrDefault<scalar>("solutionTolerance", 1e-06)
    ),
    alternativeTol_
    (
        solidProperties().lookupOrDefault<scalar>("alternativeTolerance", 1e-07)
    ),
    materialTol_
    (
        solidProperties().lookupOrDefault<scalar>("materialTolerance", 1e-05)
    ),
    infoFrequency_
    (
        solidProperties().lookupOrDefault<int>("infoFrequency", 100)
    ),
    nCorr_(solidProperties().lookupOrDefault<int>("nCorrectors", 10000)),
    maxIterReached_(0),
    residualFilePtr_(NULL),
    enforceLinear_(false),
    relaxationMethod_
    (
        solidProperties_.lookupOrDefault<word>("relaxationMethod", "fixed")
    ),
    aitkenAlpha_
    (
        IOobject
        (
            "aitkenAlpha",
            runTime.constant(),
            meshPtr_(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    aitkenResidual_
    (
        IOobject
        (
            "aitkenResidual",
            runTime.constant(),
            meshPtr_(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    QuasiNewtonRestartFreq_
    (
        solidProperties_.lookupOrDefault<int>("QuasiNewtonRestartFrequency", 25)
    ),
    QuasiNewtonV_(QuasiNewtonRestartFreq_ + 2),
    QuasiNewtonW_(QuasiNewtonRestartFreq_ + 2),
    QuasiNewtonT_(QuasiNewtonRestartFreq_ + 2),
    DRef_
    (
        IOobject
        (
            "DRef",
            runTime.constant(),
            meshPtr_(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    unrelaxedDRef_
    (
        IOobject
        (
            "unrelaxedDRef",
            runTime.constant(),
            meshPtr_(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshPtr_(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL)
{
    // Force old time fields to be stored
    D_.oldTime().oldTime();
    DD_.oldTime().oldTime();
    pointD_.oldTime();
    pointDD_.oldTime();
    gradD_.oldTime();
    gradDD_.oldTime();
    sigma_.oldTime();

    // Check that at least one of D or DD were read from file
    if
    (
       !IOobject
        (
            "D", runTime.timeName(), mesh(), IOobject::MUST_READ
        ).headerOk()
    && !IOobject
        (
            "DD", runTime.timeName(), mesh(), IOobject::MUST_READ
        ).headerOk()
    )
    {
        FatalErrorIn("solidModel::solidModel(...)")
            << "Either D or DD should be specified!" << abort(FatalError);
    }

    // Print out the relaxation factor
    Info<< "    under-relaxation method: " << relaxationMethod_ << endl;
    if (relaxationMethod_ == "QuasiNewton")
    {
        Info<< "        restart frequency: " << QuasiNewtonRestartFreq_ << endl;
    }

    // If requested, create the residual file
    if (solidProperties_.lookupOrDefault<Switch>("residualFile", false))
    {
        if (Pstream::master())
        {
            Info<< "Creating residual.dat" << endl;
            residualFilePtr_.set
            (
                new OFstream(runTime.path()/"residual.dat")
            );
        }
    }

    // If the case is axisymmetric, we will disable solving in the out-of-plane
    // direction
    checkWedges();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModel::~solidModel()
{
    thermalPtr_.clear();
    mechanicalPtr_.clear();
    deleteDemandDrivenData(globalFaceZonesPtr_);
    deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::thermalModel& Foam::solidModel::thermal() const
{
    if (thermalPtr_.empty())
    {
        makeThermalModel();
    }

    return thermalPtr_();
}


const Foam::mechanicalModel& Foam::solidModel::mechanical() const
{
    if (mechanicalPtr_.empty())
    {
        makeMechanicalModel();
    }

    return mechanicalPtr_();
}


void Foam::solidModel::setTraction
(
    const label patchID,
    const vectorField& traction
)
{
    solidModel::setTraction(solutionD().boundaryField()[patchID], traction);
}


void Foam::solidModel::setPressure
(
    const label patchID,
    const scalarField& pressure
)
{
    solidModel::setPressure(solutionD().boundaryField()[patchID], pressure);
}


Foam::vector Foam::solidModel::pointU(const label pointID) const
{
    pointVectorField pointU
    (
        IOobject
        (
            "pointU",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimVelocity, vector::zero)
    );

    mechanical().volToPoint().interpolate(U(), pointU);

    return pointU.internalField()[pointID];
}


Foam::tmp<Foam::vectorField>
Foam::solidModel::faceZonePointDisplacementIncrement
(
    const label zoneID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(),
            vector::zero
        )
    );
    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& pointDDI = pointDD().internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        //- Inter-proc points are shared by multiple procs
        //  pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if (zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] = pointDDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            //- now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] =
                zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        tPointDisplacement() =
            vectorField
            (
                pointDDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZonePointDisplacement
(
    const label zoneID
) const
{
    tmp<vectorField> tPointDisplacement
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints().size(),
            vector::zero
        )
    );
    vectorField& pointDisplacement = tPointDisplacement();

    const vectorField& oldPointDI = pointD().oldTime().internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        // Inter-proc points are shared by multiple procs
        // pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if(zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] =
                    oldPointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            // Now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] = zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        tPointDisplacement() =
            vectorField
            (
                oldPointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    return tPointDisplacement;
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZoneAcceleration
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tAcceleration
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& acceleration = tAcceleration();

    const volVectorField a = fvc::d2dt2(D());

    vectorField patchAcceleration = a.boundaryField()[patchID];

    const label patchStart = mesh().boundaryMesh()[patchID].start();

    forAll(patchAcceleration, I)
    {
        acceleration[mesh().faceZones()[zoneID].whichFace(patchStart + I)] =
            patchAcceleration[I];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(acceleration, sumOp<vectorField>());

    return tAcceleration;
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZoneVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tVelocity
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& velocity = tVelocity();

    vectorField patchVelocity = U().boundaryField()[patchID];

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchVelocity, I)
    {
        velocity[mesh().faceZones()[zoneID].whichFace(patchStart + I)] =
            patchVelocity[I];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(velocity, sumOp<vectorField>());

    return tVelocity;
}


Foam::tmp<Foam::tensorField> Foam::solidModel::faceZoneSurfaceGradientOfVelocity
(
    const label zoneID,
    const label patchID
) const
{
    tmp<tensorField> tVelocityGradient
    (
        new tensorField
        (
            mesh().faceZones()[zoneID]().size(),
            tensor::zero
        )
    );
    tensorField& velocityGradient = tVelocityGradient();

    vectorField pPointU =
        mechanical().volToPoint().interpolate
        (
            mesh().boundaryMesh()[patchID],
            U()
        );

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints +=
        pointDorPointDD().boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    tensorField patchGradU = fvc::fGrad(patch, pPointU);

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart =
            mesh().boundaryMesh()[patchID].start();

        forAll(patchGradU, i)
        {
            velocityGradient
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] = patchGradU[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(velocityGradient, sumOp<tensorField>());
    }
    else
    {
        velocityGradient = patchGradU;
    }

    return tVelocityGradient;
}


Foam::tmp<Foam::vectorField> Foam::solidModel::currentFaceZonePoints
(
    const label zoneID
) const
{
    vectorField pointDisplacement
    (
        mesh().faceZones()[zoneID]().localPoints().size(),
        vector::zero
    );

    // Lookup the point displacement field
    const vectorField& pointDI = pointDorPointDD().internalField();

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone
        const labelList& curPointMap =
            globalToLocalFaceZonePointMap()[globalZoneIndex];

        const labelList& zoneMeshPoints =
            mesh().faceZones()[zoneID]().meshPoints();

        vectorField zonePointsDisplGlobal
        (
            zoneMeshPoints.size(),
            vector::zero
        );

        // Inter-proc points are shared by multiple procs
        // pointNumProc is the number of procs which a point lies on
        scalarField pointNumProcs(zoneMeshPoints.size(), 0);

        forAll(zonePointsDisplGlobal, globalPointI)
        {
            label localPoint = curPointMap[globalPointI];

            if (zoneMeshPoints[localPoint] < mesh().nPoints())
            {
                label procPoint = zoneMeshPoints[localPoint];

                zonePointsDisplGlobal[globalPointI] = pointDI[procPoint];

                pointNumProcs[globalPointI] = 1;
            }
        }

        if (Pstream::parRun())
        {
            reduce(zonePointsDisplGlobal, sumOp<vectorField>());
            reduce(pointNumProcs, sumOp<scalarField>());

            // now average the displacement between all procs
            zonePointsDisplGlobal /= pointNumProcs;
        }

        forAll(pointDisplacement, globalPointI)
        {
            const label localPoint = curPointMap[globalPointI];

            pointDisplacement[localPoint] = zonePointsDisplGlobal[globalPointI];
        }
    }
    else
    {
        pointDisplacement =
            vectorField
            (
                pointDI,
                mesh().faceZones()[zoneID]().meshPoints()
            );
    }

    tmp<vectorField> tCurrentPoints
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().localPoints()
          + pointDisplacement
        )
    );

    return tCurrentPoints;
}


Foam::tmp<Foam::vectorField> Foam::solidModel::faceZoneNormal
(
    const label zoneID,
    const label patchID
) const
{
    tmp<vectorField> tNormals
    (
        new vectorField
        (
            mesh().faceZones()[zoneID]().size(),
            vector::zero
        )
    );
    vectorField& normals = tNormals();

    const faceList& localFaces =
        mesh().boundaryMesh()[patchID].localFaces();

    const pointVectorField& pointD = pointDorPointDD();

    // Current deformed patch points
    vectorField localPoints =
        mesh().boundaryMesh()[patchID].localPoints();
    localPoints += pointD.boundaryField()[patchID].patchInternalField();

    PrimitivePatch<face, List, const pointField&> patch
    (
        localFaces,
        localPoints
    );

    vectorField patchNormals(patch.size(), vector::zero);

    forAll(patchNormals, faceI)
    {
        patchNormals[faceI] =
            localFaces[faceI].normal(localPoints);
    }

    label globalZoneIndex = findIndex(globalFaceZones(), zoneID);

    if (globalZoneIndex != -1)
    {
        // global face zone

        const label patchStart =
            mesh().boundaryMesh()[patchID].start();

        forAll(patchNormals, i)
        {
            normals
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ] = patchNormals[i];
        }

        // Parallel data exchange: collect field on all processors
        reduce(normals, sumOp<vectorField>());
    }
    else
    {
        normals = patchNormals;
    }

    return tNormals;
}


void Foam::solidModel::updateTotalFields()
{
    mechanical().updateTotalFields();
}


void Foam::solidModel::end()
{
    if (maxIterReached_ > 0)
    {
        WarningIn(type() + "::end()")
            << "The maximum momentum correctors were reached in "
            << maxIterReached_ << " time-steps" << nl << endl;
    }
    else
    {
        Info<< "The momentum equation converged in all time-steps"
            << nl << endl;
    }
}


Foam::autoPtr<Foam::solidModel> Foam::solidModel::New
(
    Time& runTime,
    const word& region
)
{
    word solidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the flow is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                runTime.caseConstant()/region,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidModel")
            >> solidModelTypeName;
    }

    Info<< "Selecting solidModel " << solidModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::New(Time&, const word&)"
        )   << "Unknown solidModel type " << solidModelTypeName
            << endl << endl
            << "Valid solidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidModel>(cstrIter()(runTime, region));
}


const Foam::labelList& Foam::solidModel::globalFaceZones() const
{
    if (!globalFaceZonesPtr_)
    {
        calcGlobalFaceZones();
    }

    return *globalFaceZonesPtr_;
}


const Foam::labelListList&
Foam::solidModel::globalToLocalFaceZonePointMap() const
{
    if (!globalToLocalFaceZonePointMapPtr_)
    {
        calcGlobalToLocalFaceZonePointMap();
    }

    return *globalToLocalFaceZonePointMapPtr_;
}


void Foam::solidModel::setTraction
(
    fvPatchVectorField& tractionPatch,
    const vectorField& traction
)
{
    if (tractionPatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(tractionPatch);

        patchD.traction() = traction;
    }
    else if
    (
        tractionPatch.type() == blockSolidTractionFvPatchVectorField::typeName
    )
    {
        blockSolidTractionFvPatchVectorField& patchD =
            refCast<blockSolidTractionFvPatchVectorField>(tractionPatch);

        patchD.traction() = traction;
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
            << "for patch " << tractionPatch.patch().name()
            << " should instead be type "
            << solidTractionFvPatchVectorField::typeName
            << " or "
            << blockSolidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }
}


void Foam::solidModel::setPressure
(
    fvPatchVectorField& pressurePatch,
    const scalarField& pressure
)
{
    if (pressurePatch.type() == solidTractionFvPatchVectorField::typeName)
    {
        solidTractionFvPatchVectorField& patchD =
            refCast<solidTractionFvPatchVectorField>(pressurePatch);

        patchD.pressure() = pressure;
    }
    else if
    (
        pressurePatch.type() == blockSolidTractionFvPatchVectorField::typeName
    )
    {
        blockSolidTractionFvPatchVectorField& patchD =
            refCast<blockSolidTractionFvPatchVectorField>(pressurePatch);

        patchD.pressure() = pressure;
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
            << solidTractionFvPatchVectorField::typeName
            << " or "
            << blockSolidTractionFvPatchVectorField::typeName
            << abort(FatalError);
    }
}


void Foam::solidModel::setTraction
(
    const label patchID,
    const label zoneID,
    const vectorField& faceZoneTraction
)
{
    vectorField patchTraction(mesh().boundary()[patchID].size(), vector::zero);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchTraction, i)
    {
        patchTraction[i] =
            faceZoneTraction
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setTraction(patchID, patchTraction);
}


void Foam::solidModel::setPressure
(
    const label patchID,
    const label zoneID,
    const scalarField& faceZonePressure
)
{
    scalarField patchPressure(mesh().boundary()[patchID].size(), 0.0);

    const label patchStart =
        mesh().boundaryMesh()[patchID].start();

    forAll(patchPressure, i)
    {
        patchPressure[i] =
            faceZonePressure
            [
                mesh().faceZones()[zoneID].whichFace(patchStart + i)
            ];
    }

    setPressure(patchID, patchPressure);
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const volScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


Foam::Switch& Foam::solidModel::checkEnforceLinear(const surfaceScalarField& J)
{
    scalar minJ = min(J).value();
    reduce(minJ, minOp<scalar>());

    scalar maxJ = max(J).value();
    reduce(maxJ, maxOp<scalar>());

    if ((minJ < 0.01) || (maxJ > 100))
    {
        Info<< "Enforcing linear geometry: "
            << "minJ: " << minJ << ", maxJ: " << maxJ << endl;

        // Enable enforce linear to try improve convergence
        enforceLinear() = true;
    }

    return enforceLinear();
}


void Foam::solidModel::writeFields(const Time& runTime)
{
    // Calculaute equivalent (von Mises) stress
    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma())))
    );

    Info<< "Max sigmaEq = " << gMax(sigmaEq) << endl;

    runTime.write();
}


Foam::scalar Foam::solidModel::newDeltaT()
{
    return mechanical().newDeltaT();
}

void Foam::solidModel::moveMesh
(
    const pointField& oldPoints,
    const volVectorField& DD,
    pointVectorField& pointDD
)
{
    Info<< "Moving the mesh to the deformed configuration" << nl << endl;

    //- Move mesh by interpolating displacement field to vertices
    // To be checked: sync boundary and global points across procs to make sure
    // numiercal error does not build up and when end up with the error
    // "face area does not match neighbour..."
    // We could sync points as a pointVectorField just as we sync pointDD

    // Interpolate cell displacements to vertices
    mechanical().interpolate(DD, pointDD);

    // Ensure continuous displacement across processor boundary
    // Something strange is happening here
    pointDD.correctBoundaryConditions();

    vectorField& pointDDI = pointDD.internalField();

    vectorField newPoints = oldPoints;

    // Correct symmetryPlane points

    forAll(mesh().boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());

            const vector i(1, 0, 0);
            const vector j(0, 1, 0);
            const vector k(0, 0, 1);

            if (mag(avgN & i) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].x() = 0;
                }
            }
            else if (mag(avgN & j) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].y() = 0;
                }
            }
            else if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
        else if (isA<emptyPolyPatch>(mesh().boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const vector avgN =
                gAverage(mesh().boundaryMesh()[patchI].pointNormals());
            const vector k(0, 0, 1);

            if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDDI[meshPoints[pI]].z() = 0;
                }
            }
        }
    }

    // Note: allPoints will have more points than pointDD if there are
    // globalFaceZones
    forAll (pointDDI, pointI)
    {
        newPoints[pointI] += pointDDI[pointI];
    }

    // Move unused globalFaceZone points
    updateGlobalFaceZoneNewPoints(pointDDI, newPoints);

    twoDPointCorrector twoDCorrector(mesh());
    twoDCorrector.correctPoints(newPoints);
    twoDCorrector.correctPoints(pointDDI);
    mesh().movePoints(newPoints);
    mesh().V00();
    mesh().moving(false);
    mesh().changing(false);
    mesh().setPhi().writeOpt() = IOobject::NO_WRITE;

    // Tell the mechanical model to move the subMeshes, if they exist
    mechanical().moveSubMeshes();
}


bool Foam::solidModel::read()
{
    if (regIOobject::read())
    {
        solidProperties_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
