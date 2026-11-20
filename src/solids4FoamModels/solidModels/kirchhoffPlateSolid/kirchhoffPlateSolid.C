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

#ifndef OPENFOAM_ORG

#include "kirchhoffPlateSolid.H"

// PETSc must be included before faCFD.H, which brings the Foam namespace
// into scope
#if defined(USE_PETSC) && defined(OPENFOAM_COM)
#include <petscksp.h>
#include <vector>
#endif

#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "linearElasticMechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawManager.H"
#include "sigFpe.H"
#include "leastSquaresFaVectors.H"
#include "fixedGradientFaPatchFields.H"
#include "freeEdgeStencil.H"

//#include "BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kirchhoffPlateSolid, 0);
addToRunTimeSelectionTable(solidModel, kirchhoffPlateSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool kirchhoffPlateSolid::converged
(
    const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
    const SolverPerformance<scalar>& solverPerfM,
    const SolverPerformance<scalar>& solverPerfw,
#else
    const lduSolverPerformance& solverPerfM,
    const lduSolverPerformance& solverPerfw,
#endif
    const areaScalarField& M,
    const areaScalarField& w
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar residualM =
        gMax
        (
            (
                mag(M - M.prevIter())
               /max
                (
                    gMax(mag(M - M.oldTime())()), SMALL
                )
            )()
        );

    const scalar residualw =
        gMax
        (
            (
                mag(w - w.prevIter())
               /max
                (
                    gMax(mag(w - w.oldTime())()), SMALL
                )
            )()
        );

    // This model admits only a linearElastic material, whose material
    // residual is zero, so it is not asked for and not tested
    if (iCorr > 1)
    {
        bool convergedM = false;
        bool convergedw = false;

        if
        (
            (
                solverPerfM.initialResidual() < solutionTol()
             && residualM < solutionTol()
            )
         || solverPerfM.initialResidual() < alternativeTol()
         || residualM < alternativeTol()
        )
        {
            convergedM = true;
        }

        if
        (
            (
                solverPerfw.initialResidual() < solutionTol()
             && residualw < solutionTol()
            )
         || solverPerfw.initialResidual() < alternativeTol()
         || residualw < alternativeTol()
        )
        {
            convergedw = true;
        }


        if (convergedM && convergedw)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (M & w), relRes (M & w), iters (M & w)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfM.initialResidual()
            << ", " << solverPerfw.initialResidual()
            << ", " << residualM
            << ", " << residualw
            << ", " << solverPerfM.nIterations()
            << ", " << solverPerfw.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the M-w loop" << endl;
    }

    return converged;
}


const fvPatch& kirchhoffPlateSolid::areaPatch() const
{
    if (areaPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaPatchID_];
}


const fvPatch& kirchhoffPlateSolid::areaShadowPatch() const
{
    if (areaShadowPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaShadowPatchID_];
}


void kirchhoffPlateSolid::calcAreaPatches() const
{
    // Note: face0PatchID may be -1 if this processor has no faces on the
    // finiteArea patch

    // Check that all areaMesh faces map to the same patch

    const polyMesh& pMesh = mesh();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();
    const labelList& faceLabels = aMesh_.faceLabels();
    const label pMeshNFaces = pMesh.nFaces();

    if (faceLabels.size() > 0)
    {
        const label face0ID = faceLabels[0];

        if (face0ID < pMeshNFaces)
        {
            areaPatchID_ = bm.whichPatch(face0ID);

            // Check all faces map to the same fvMesh patch

            forAll(faceLabels, aFaceI)
            {
                const label faceID = faceLabels[aFaceI];

                // Escape if face is beyond active faces, eg belongs to a face
                // zone
                if (faceID < pMeshNFaces)
                {
                    const label curPatchID = bm.whichPatch(face0ID);

                    if (curPatchID != areaPatchID_)
                    {
                        FatalErrorIn
                        (
                            "void kirchhoffPlateSolid::calcAreaPatches() const"
                        )   << "The finiteArea patch should correspond to a "
                            << "patch on the boundary of the polyMesh!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }


    // We will now check if the polyMesh has the same number of cells as the
    // number of faces on the areaPatch, as we are assuming the polyMesh to be
    // one cell thick
    if (pMesh.nCells() != pMesh.boundaryMesh()[areaPatchID_].size())
    {
        FatalErrorIn
        (
            "void kirchhoffPlateSolid::calcAreaPatches() const"
        )   << "The solid polyMesh should be one cell thick, where there is "
            << "the same number of cells as the number of faces on the "
            << "areaPatch" << endl
            << "areaPatchID: " << areaPatchID_
            << abort(FatalError);
    }


    // To find the areaShadowPatch, we will ...
    const labelUList& faceCells =
        pMesh.boundaryMesh()[areaPatchID_].faceCells();

    if (faceCells.size())
    {
        const cellList& cells = pMesh.cells();

        const label face0ID = bm[areaPatchID_].start();
        const label cell0ID = faceCells[0];
        const vector& face0N = bm[areaPatchID_].faceNormals()[0];
        const labelList& curCellFaces = cells[cell0ID];

        scalar mostNegativeDotProduct = GREAT;

        forAll(curCellFaces, fI)
        {
            const label curFaceID = curCellFaces[fI];

            if (curFaceID != face0ID)
            {
                if (!pMesh.isInternalFace(curFaceID))
                {
                    const label otherPatchID = bm.whichPatch(curFaceID);
                    const label curLocalFaceID =
                        curFaceID - bm[otherPatchID].start();

                    const vector& curFaceN =
                        bm[otherPatchID].faceNormals()[curLocalFaceID];

                    const scalar dotProduct = face0N & curFaceN;

                    if (dotProduct < mostNegativeDotProduct)
                    {
                        mostNegativeDotProduct = dotProduct;
                        areaShadowPatchID_ = otherPatchID;
                    }
                }
            }
        }
    }


    // Check if the areaPatch and areaShadowPatch have the same number of faces
    if
    (
        pMesh.boundaryMesh()[areaShadowPatchID_].size()
     != pMesh.boundaryMesh()[areaPatchID_].size()
    )
    {
        FatalErrorIn
        (
            "void kirchhoffPlateSolid::calcAreaPatches() const"
        )   << "The polyMesh should be one cell thick, where there should be "
            << "two patches opposite each other that have the same number of "
            << "faces" << abort(FatalError);
    }
}


void kirchhoffPlateSolid::solveCoupled()
{
#if defined(USE_PETSC) && defined(OPENFOAM_COM)

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "The coupled Kirchhoff plate solver does not run in parallel"
            << abort(FatalError);
    }

    const faMesh& mesh = aMesh_;
    const label nCells = w_.size();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& S = mesh.S();
    const faBoundaryMesh& bm = mesh.boundary();
    const scalar D = bendingStiffness_.value();
    const scalar nu = nu_.value();

    // Unknown numbering: M and W = D*w of cell i at 2i and 2i + 1, followed
    // by M and W of each boundary edge
    labelList bStart(bm.size(), -1);
    label nBoundaryEdges = 0;
    forAll(bm, patchI)
    {
        if (w_.boundaryField()[patchI].size())
        {
            bStart[patchI] = 2*nCells + 2*nBoundaryEdges;
            nBoundaryEdges += bm[patchI].size();
        }
    }
    const label nUnknowns = 2*nCells + 2*nBoundaryEdges;

    auto mCell = [](const label i) { return 2*i; };
    auto wCell = [](const label i) { return 2*i + 1; };
    auto mEdge = [&bStart](const label patchI, const label e)
    {
        return bStart[patchI] + 2*e;
    };
    auto wEdge = [&bStart](const label patchI, const label e)
    {
        return bStart[patchI] + 2*e + 1;
    };

    // Sparse matrix in triplet form
    DynamicList<label> rows(20*nUnknowns);
    DynamicList<label> cols(20*nUnknowns);
    DynamicList<scalar> vals(20*nUnknowns);
    auto add = [&](const label r, const label c, const scalar v)
    {
        rows.append(r);
        cols.append(c);
        vals.append(v);
    };

    // Right-hand side of the boundary rows, which does not change between
    // outer iterations
    scalarField rhsBoundary(nUnknowns, 0.0);

    // Least-squares gradient stencil of w in each cell, as the segregated
    // algorithm uses the leastSquares gradient for theta = -grad(w)
    const leastSquaresFaVectors& lsv = leastSquaresFaVectors::New(mesh);
    const edgeVectorField& pV = lsv.pVectors();
    const edgeVectorField& nV = lsv.nVectors();

    List<DynamicList<label>> gradIdx(nCells);
    List<DynamicList<vector>> gradVec(nCells);

    forAll(own, e)
    {
        const label o = own[e];
        const label n = nei[e];

        gradIdx[o].append(wCell(n));
        gradVec[o].append(pV[e]);
        gradIdx[o].append(wCell(o));
        gradVec[o].append(-pV[e]);

        gradIdx[n].append(wCell(n));
        gradVec[n].append(-nV[e]);
        gradIdx[n].append(wCell(o));
        gradVec[n].append(nV[e]);
    }

    forAll(bm, patchI)
    {
        if (bStart[patchI] < 0)
        {
            continue;
        }

        const labelUList& ef = bm[patchI].edgeFaces();
        const vectorField& pVb = pV.boundaryField()[patchI];

        forAll(ef, e)
        {
            gradIdx[ef[e]].append(wEdge(patchI, e));
            gradVec[ef[e]].append(pVb[e]);
            gradIdx[ef[e]].append(wCell(ef[e]));
            gradVec[ef[e]].append(-pVb[e]);
        }
    }

    // Add coeff*(n & grad(W))_cellI to the given row
    auto addNormalGrad = [&]
    (
        const label row, const label cellI, const vector& n, const scalar coeff
    )
    {
        forAll(gradIdx[cellI], k)
        {
            add(row, gradIdx[cellI][k], coeff*(n & gradVec[cellI][k]));
        }
    };

    // Add coeff*(n & grad(M))_cellI to the given row: the M unknowns use the
    // same stencil as W, with the index of M one before that of W
    auto addNormalGradM = [&]
    (
        const label row, const label cellI, const vector& n, const scalar coeff
    )
    {
        forAll(gradIdx[cellI], k)
        {
            add(row, gradIdx[cellI][k] - 1, coeff*(n & gradVec[cellI][k]));
        }
    };

    // Interior Laplacian coefficients, taken from the finite-area
    // discretisation so that the schemes of faSchemes are used
    const faScalarMatrix lapM(fam::laplacian(M_));
    const faScalarMatrix lapW(fam::laplacian(bendingStiffness_, w_));
    const scalarField& upperM = lapM.upper();
    const scalarField& upperW = lapW.upper();

    // Non-orthogonal correction of the Laplacians: the corrected scheme adds
    // magLe*(k & grad(phi))_e to the edge flux, with the gradient linearly
    // interpolated from the least-squares cell gradients. It is assembled
    // implicitly through the gradient stencils, so that one solve gives the
    // corrected discretisation on any mesh. The correction vectors vanish
    // on an orthogonal mesh and on the boundary edges
    const edgeVectorField& corrVecs = mesh.correctionVectors();
    const edgeScalarField& weights = mesh.weights();
    const scalarField& magLe = mesh.magLe().primitiveField();

    // M equation: -laplacian(M) + rho*h*d2dt2(w) = p
    // w equation: laplacian(D, w) + M = 0, written for W = D*w
    forAll(own, e)
    {
        const label o = own[e];
        const label n = nei[e];

        add(mCell(o), mCell(o), upperM[e]);
        add(mCell(o), mCell(n), -upperM[e]);
        add(mCell(n), mCell(n), upperM[e]);
        add(mCell(n), mCell(o), -upperM[e]);

        add(wCell(o), wCell(o), -upperW[e]/D);
        add(wCell(o), wCell(n), upperW[e]/D);
        add(wCell(n), wCell(n), -upperW[e]/D);
        add(wCell(n), wCell(o), upperW[e]/D);

        if (magSqr(corrVecs[e]) > SMALL)
        {
            const vector k = corrVecs[e];
            const scalar lo = magLe[e]*weights[e];
            const scalar ln = magLe[e]*(1.0 - weights[e]);

            addNormalGrad(wCell(o), o, k, lo);
            addNormalGrad(wCell(o), n, k, ln);
            addNormalGrad(wCell(n), o, k, -lo);
            addNormalGrad(wCell(n), n, k, -ln);

            addNormalGradM(mCell(o), o, k, -lo);
            addNormalGradM(mCell(o), n, k, -ln);
            addNormalGradM(mCell(n), o, k, lo);
            addNormalGradM(mCell(n), n, k, ln);
        }
    }

    // Inertia with the Euler scheme, consistent with the segregated
    // algorithm: rho*h*(w - 2*w0 + w00)/deltaT^2
    const scalar inertiaCoeff =
        rho_.value()*h_.value()/sqr(runTime().deltaTValue());

    forAll(S, i)
    {
        add(wCell(i), mCell(i), S[i]);

        if (rho_.value() > SMALL)
        {
            add(mCell(i), wCell(i), S[i]*inertiaCoeff/D);
        }
    }

    // Boundary edge fluxes and boundary rows
    forAll(bm, patchI)
    {
        if (bStart[patchI] < 0)
        {
            continue;
        }

        const faPatch& fap = bm[patchI];
        const labelUList& ef = fap.edgeFaces();
        const scalarField& dc = fap.deltaCoeffs();
        const scalarField& Lb = mesh.magLe().boundaryField()[patchI];
        const vectorField n(fap.edgeNormals());
        const vectorField& Ce = fap.edgeCentres();
        const pointField& points = mesh.points();

        const faPatchField<scalar>& Mp = M_.boundaryField()[patchI];
        const faPatchField<scalar>& wp = w_.boundaryField()[patchI];

        // Fluxes through the boundary edges in the cell rows
        forAll(ef, e)
        {
            const label c = ef[e];

            add(mCell(c), mCell(c), Lb[e]*dc[e]);
            add(mCell(c), mEdge(patchI, e), -Lb[e]*dc[e]);

            add(wCell(c), wCell(c), -Lb[e]*dc[e]);
            add(wCell(c), wEdge(patchI, e), Lb[e]*dc[e]);
        }

        // Edge connectivity for the free edge conditions
        const labelListList nbrs(freeEdgeStencil::edgeNeighbours(fap));
        const labelListList ends(freeEdgeStencil::edgeEndPoints(fap));
        const vectorField t(freeEdgeStencil::edgeTangents(fap, nbrs));

        // Add coeff*g_j to the row, where g_j = (W_b - W_P)*deltaCoeffs is
        // D times the normal slope on edge j
        auto addSlope = [&](const label row, const label j, const scalar coeff)
        {
            add(row, wEdge(patchI, j), coeff*dc[j]);
            add(row, wCell(ef[j]), -coeff*dc[j]);
        };

        // Add coeff*q_j to the row, where q_j = D*kappa*dw/ds on edge j from
        // the rotation of the neighbouring edge normals and the neighbouring
        // edge values of W
        auto addCurvatureSlope = [&]
        (
            const label row, const label j, const scalar coeff
        )
        {
            const labelList& nbr = nbrs[j];

            if (nbr.size() == 0)
            {
                return;
            }

            label m = j;
            label p = j;

            forAll(nbr, i)
            {
                if (((Ce[nbr[i]] - Ce[j]) & t[j]) > 0)
                {
                    p = nbr[i];
                }
                else
                {
                    m = nbr[i];
                }
            }

            const scalar h = mag(Ce[p] - Ce[m]);
            const scalar f = coeff*(t[j] & (n[p] - n[m]))/sqr(h);

            add(row, wEdge(patchI, p), f);
            add(row, wEdge(patchI, m), -f);
        };

        // M rows
        forAll(ef, e)
        {
            const label c = ef[e];
            const label row = mEdge(patchI, e);

            if (Mp.type() == "clampedMoment")
            {
                // M_b = -D*d2w/dn2 = -(0 - n & grad(W)_P)*dc, as the slope
                // vanishes at the edge
                add(row, row, 1.0);
                addNormalGrad(row, c, n[e], -dc[e]);
            }
            else if (Mp.type() == "freeEdgeMoment")
            {
                // (M_b - M_P)*dc = (1 - nu)*(d2g/ds2 - dq/ds)
                add(row, row, dc[e]);
                add(row, mCell(c), -dc[e]);

                const labelList& nbr = nbrs[e];

                // Distances behind (m) and ahead (p) of the edge centre and
                // the edges, or virtual corner node, that hold the values
                scalar hm = 0.0;
                scalar hp = 0.0;
                label jm = -1;
                label jp = -1;

                if (nbr.size() == 2)
                {
                    jm = nbr[0];
                    jp = nbr[1];

                    if (((Ce[jp] - Ce[jm]) & t[e]) < 0)
                    {
                        Swap(jm, jp);
                    }

                    hm = mag(Ce[jm] - Ce[e]);
                    hp = mag(Ce[jp] - Ce[e]);
                }
                else if (nbr.size() == 1 && ends[e].size() == 1)
                {
                    const label nbrI = nbr[0];
                    const label pointI = ends[e][0];
                    const scalar h0 = mag(points[pointI] - Ce[e]);

                    // Virtual node at the corner: zero at a clamped or
                    // supported corner (jc = -1), mirror of this edge at a
                    // free corner (jc = e)
                    scalar hc = h0;
                    label jc = -1;

                    if (!freeEdgeStencil::fixesValueAtPoint(w_, fap, pointI))
                    {
                        hc = 2.0*h0;
                        jc = e;
                    }

                    if (((Ce[nbrI] - Ce[e]) & t[e]) > 0)
                    {
                        hm = hc;
                        jm = jc;
                        hp = mag(Ce[nbrI] - Ce[e]);
                        jp = nbrI;
                    }
                    else
                    {
                        hm = mag(Ce[nbrI] - Ce[e]);
                        jm = nbrI;
                        hp = hc;
                        jp = jc;
                    }
                }
                else
                {
                    continue;
                }

                const scalar denom = hm*hp*(hm + hp);
                const scalar K = 1.0 - nu;

                // d2g/ds2 = 2*(hm*g_p - (hm + hp)*g_e + hp*g_m)/denom
                addSlope(row, e, 2.0*K*(hm + hp)/denom);
                if (jp != -1)
                {
                    addSlope(row, jp, -2.0*K*hm/denom);
                }
                if (jm != -1)
                {
                    addSlope(row, jm, -2.0*K*hp/denom);
                }

                // dq/ds = (hm^2*q_p - hp^2*q_m + (hp^2 - hm^2)*q_e)/denom
                addCurvatureSlope(row, e, K*(sqr(hp) - sqr(hm))/denom);
                if (jp != -1)
                {
                    addCurvatureSlope(row, jp, K*sqr(hm)/denom);
                }
                if (jm != -1)
                {
                    addCurvatureSlope(row, jm, -K*sqr(hp)/denom);
                }
            }
            else if (Mp.fixesValue())
            {
                add(row, row, 1.0);
                rhsBoundary[row] = Mp[e];
            }
            else if (isA<fixedGradientFaPatchField<scalar>>(Mp))
            {
                add(row, row, dc[e]);
                add(row, mCell(c), -dc[e]);
                rhsBoundary[row] =
                    refCast<const fixedGradientFaPatchField<scalar>>(Mp)
                   .gradient()[e];
            }
            else if (Mp.type() == "zeroGradient")
            {
                add(row, row, 1.0);
                add(row, mCell(c), -1.0);
            }
            else
            {
                FatalErrorInFunction
                    << "The coupled solver does not support the " << Mp.type()
                    << " condition for M on patch " << fap.name()
                    << abort(FatalError);
            }
        }

        // w rows
        forAll(ef, e)
        {
            const label c = ef[e];
            const label row = wEdge(patchI, e);

            if (wp.type() == "freeEdgeDisplacement")
            {
                // (W_b - W_P)*dc = n & grad(W)_P + nu/(1 - nu)*M_P/dc
                // where at a free-free corner the right-hand side is that of
                // the neighbouring edge
                add(row, row, dc[e]);
                add(row, wCell(c), -dc[e]);

                label j = e;

                if
                (
                    nbrs[e].size() == 1
                 && ends[e].size() == 1
                 && !freeEdgeStencil::fixesValueAtPoint(w_, fap, ends[e][0])
                )
                {
                    j = nbrs[e][0];
                }

                addNormalGrad(row, ef[j], n[j], -1.0);
                add(row, mCell(ef[j]), -nu/((1.0 - nu)*dc[j]));
            }
            else if (wp.fixesValue())
            {
                add(row, row, 1.0);
                rhsBoundary[row] = D*wp[e];
            }
            else if (isA<fixedGradientFaPatchField<scalar>>(wp))
            {
                add(row, row, dc[e]);
                add(row, wCell(c), -dc[e]);
                rhsBoundary[row] =
                    D*refCast<const fixedGradientFaPatchField<scalar>>(wp)
                   .gradient()[e];
            }
            else if (wp.type() == "zeroGradient")
            {
                add(row, row, 1.0);
                add(row, wCell(c), -1.0);
            }
            else
            {
                FatalErrorInFunction
                    << "The coupled solver does not support the " << wp.type()
                    << " condition for w on patch " << fap.name()
                    << abort(FatalError);
            }
        }
    }

    // PETSc matrix
    PetscBool petscInitialised;
    PetscInitialized(&petscInitialised);
    if (!petscInitialised)
    {
        PetscInitialize(NULL, NULL, NULL, NULL);
    }

    // Optional options file in the case directory
    {
        fileName optionsFile
        (
            solidModelDict().lookupOrDefault<fileName>
            (
                "optionsFile", "petscOptions"
            )
        );
        optionsFile.expand();

        if (isFile(optionsFile))
        {
            Info<< "    Reading PETSc options from " << optionsFile << endl;
            PetscOptionsInsertFile
            (
                PETSC_COMM_SELF, NULL, optionsFile.c_str(), PETSC_TRUE
            );
        }
    }

    std::vector<PetscInt> nnz(nUnknowns, 0);
    forAll(rows, k)
    {
        ++nnz[rows[k]];
    }

    Mat A;
    MatCreate(PETSC_COMM_SELF, &A);
    MatSetSizes(A, nUnknowns, nUnknowns, nUnknowns, nUnknowns);
    MatSetType(A, MATSEQAIJ);
    MatSeqAIJSetPreallocation(A, 0, nnz.data());

    forAll(rows, k)
    {
        MatSetValue(A, rows[k], cols[k], vals[k], ADD_VALUES);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    Vec b;
    Vec x;
    VecCreateSeq(PETSC_COMM_SELF, nUnknowns, &b);
    VecDuplicate(b, &x);

    // Direct solver by default, with MUMPS when PETSc has it: the PETSc LU
    // factorisation does not pivot and can be inaccurate for this system.
    // The options file can select another solver or package
    KSP ksp;
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPPREONLY);
    {
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCLU);

        PetscBool hasMumps;
        PetscHasExternalPackage("mumps", &hasMumps);
        if (hasMumps)
        {
            PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        }
        else
        {
            WarningInFunction
                << "PETSc was built without MUMPS: the PETSc LU "
                << "factorisation does not pivot and may be inaccurate; "
                << "check the reported true residual" << endl;
        }
    }
    KSPSetFromOptions(ksp);

    Vec r;
    VecDuplicate(b, &r);

    Info<< "    Coupled system: " << nUnknowns << " unknowns, "
        << rows.size() << " coefficients" << endl;

    // The system is linear and fully implicit, so the second iteration only
    // confirms the solution of the first
    const areaScalarField& w0 = w_.oldTime();
    const areaScalarField& w00 = w0.oldTime();

    for (int iter = 0; iter < nCorr(); ++iter)
    {
        scalarField rhs(rhsBoundary);
        forAll(S, i)
        {
            rhs[mCell(i)] = S[i]*p_[i];
            rhs[wCell(i)] = 0.0;

            if (rho_.value() > SMALL)
            {
                rhs[mCell(i)] += S[i]*inertiaCoeff*(2.0*w0[i] - w00[i]);
            }
        }

        {
            PetscScalar* bArray;
            VecGetArray(b, &bArray);
            forAll(rhs, k)
            {
                bArray[k] = rhs[k];
            }
            VecRestoreArray(b, &bArray);
        }

        // The factorisation packages can raise harmless floating point
        // exceptions, so the OpenFOAM trap is suspended during the solve
        {
            sigFpe::ignore ignoreFpe;
            KSPSolve(ksp, b, x);
        }

        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt nIters;
        KSPGetIterationNumber(ksp, &nIters);

        if (reason < 0)
        {
            FatalErrorInFunction
                << "PETSc KSP diverged, reason " << label(reason)
                << abort(FatalError);
        }

        // True residual of the linear solve
        PetscReal normR;
        PetscReal normB;
        MatMult(A, x, r);
        VecAXPY(r, -1.0, b);
        VecNorm(r, NORM_2, &normR);
        VecNorm(b, NORM_2, &normB);
        const scalar relResidual = normR/(normB + VSMALL);

        if (relResidual > 1e-4)
        {
            FatalErrorInFunction
                << "The linear solve is inaccurate: true relative residual "
                << relResidual << ". Select another factorisation package or "
                << "solver in the PETSc options file" << abort(FatalError);
        }

        // Retrieve the solution and measure the change
        const PetscScalar* xArray;
        VecGetArrayRead(x, &xArray);

        scalar maxDeltaM = 0.0;
        scalar maxDeltaW = 0.0;
        scalar maxM = SMALL;
        scalar maxW = SMALL;

        forAll(S, i)
        {
            const scalar Mi = xArray[mCell(i)];
            const scalar wi = xArray[wCell(i)]/D;
            maxDeltaM = max(maxDeltaM, mag(Mi - M_[i]));
            maxDeltaW = max(maxDeltaW, mag(wi - w_[i]));
            maxM = max(maxM, mag(Mi));
            maxW = max(maxW, mag(wi));
            M_[i] = Mi;
            w_[i] = wi;
        }

        forAll(bm, patchI)
        {
            if (bStart[patchI] < 0)
            {
                continue;
            }

            const faPatch& fap = bm[patchI];
            const labelUList& ef = fap.edgeFaces();
            const scalarField& dc = fap.deltaCoeffs();

            scalarField Mb(fap.size());
            scalarField wb(fap.size());
            forAll(ef, e)
            {
                Mb[e] = xArray[mEdge(patchI, e)];
                wb[e] = xArray[wEdge(patchI, e)]/D;
            }

            faPatchField<scalar>& Mp = M_.boundaryFieldRef()[patchI];
            faPatchField<scalar>& wp = w_.boundaryFieldRef()[patchI];

            Mp == Mb;
            wp == wb;

            // Store the slopes so that fixedGradient patches evaluate to
            // the coupled solution
            if (isA<fixedGradientFaPatchField<scalar>>(Mp))
            {
                refCast<fixedGradientFaPatchField<scalar>>(Mp).gradient() =
                    (Mb - Mp.patchInternalField())*dc;
            }
            if (isA<fixedGradientFaPatchField<scalar>>(wp))
            {
                refCast<fixedGradientFaPatchField<scalar>>(wp).gradient() =
                    (wb - wp.patchInternalField())*dc;
            }
        }

        VecRestoreArrayRead(x, &xArray);

        // Update the angle of rotation and its gradient
        theta_ = -fac::grad(w_);
        gradTheta_ = fac::grad(theta_);

        const scalar relDeltaM = maxDeltaM/maxM;
        const scalar relDeltaW = maxDeltaW/maxW;

        Info<< "    Coupled iteration " << iter + 1
            << ", KSP iterations: " << label(nIters)
            << ", true residual: " << relResidual
            << ", max relative change M: " << relDeltaM
            << ", w: " << relDeltaW << endl;

        if (relDeltaM < solutionTol() && relDeltaW < solutionTol())
        {
            break;
        }
    }

    KSPDestroy(&ksp);
    VecDestroy(&r);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);

#else

    FatalErrorInFunction
        << "The coupled Kirchhoff plate solver requires OpenFOAM.com and "
        << "solids4foam compiled with PETSc (PETSC_DIR set)"
        << abort(FatalError);

#endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kirchhoffPlateSolid::kirchhoffPlateSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    aMesh_(mesh()),
    w_
    (
        IOobject
        (
            "w",
            runTime.timeName(),
#if (OPENFOAM > 2312)
            aMesh_,
#else
            mesh(),
#endif
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    wVf_
    (
        IOobject
        (
            "wVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime.timeName(),
#if (OPENFOAM > 2312)
            aMesh_,
#else
            mesh(),
#endif
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    MVf_
    (
        IOobject
        (
            "MVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure/dimArea, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
#if (OPENFOAM > 2312)
            aMesh_,
#else
            mesh(),
#endif
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    pVf_
    (
        IOobject
        (
            "pVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", p_.dimensions(), 0.0)
    ),
    theta_
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
#if (OPENFOAM > 2312)
            aMesh_,
#else
            mesh(),
#endif
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    thetaVf_
    (
        IOobject
        (
            "thetaVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    gradTheta_(fac::grad(theta_)),
    rho_("zero", dimDensity, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    h_(solidModelDict().lookup("plateThickness")),
    bendingStiffness_("zero", dimPressure*dimVolume, 0.0),
    areaPatchID_(-1),
    areaShadowPatchID_(-1),
    coupled_(solidModelDict().lookupOrDefault<Switch>("coupled", false))
{
    // A Kirchhoff plate has one thickness and one bending stiffness, so one
    // isotropic linear elastic material is what its formulation is written
    // for
    const mechanicalConstitutiveLaw& law = mechanicalManager().singleLaw();

    if (!isA<linearElasticMechanicalConstitutiveLaw>(law))
    {
        FatalErrorInFunction
            << type() << " can only be used with the linearElastic "
            << "mechanical constitutive law" << abort(FatalError);
    }

    const linearElasticMechanicalConstitutiveLaw& mech =
        refCast<const linearElasticMechanicalConstitutiveLaw>(law);

    rho_ = mech.rho();
    E_ = mech.E();
    nu_ = mech.nu();
    bendingStiffness_ = E_*pow(h_, 3)/(12*(1 - pow(nu_, 2)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool kirchhoffPlateSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh_);

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        SolverPerformance<scalar> solverPerfM;
        SolverPerformance<scalar> solverPerfw;
        SolverPerformance<scalar>::debug = 0;
#else
        lduSolverPerformance solverPerfM;
        lduSolverPerformance solverPerfw;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the Kirchhoff plate equation for w and M" << endl;

        if (coupled_)
        {
            // Solve the M and w equations as one linear system
            solveCoupled();
        }
        else
        {
            // Segregated w and M equation loop
            do
            {
                // Algorithm
                // Solve M equation
                // Solve w equation
                // where
                // M is the moment sum
                // w is the transvere (out of plane) displacement
                // The M equation is:
                //     rho*h*fac::d2dt2(w) = fam::laplacian(M) + p
                // and the w equation is:
                //     fam::laplacian(D, w) + M
                // where
                // rho is the density
                // h is the plate thickness
                // p is the net transverse pressure
                // D is the bending stiffness
                // THESE ARE BOTH SCALARS: IDEA: USE BLOCK COUPLED!

                // Store fields for under-relaxation and residual calculation
                M_.storePrevIter();

                // Solve M equation
                // d2dt2 not implemented so calculate it manually using ddt
                // Also, "==" complains so we will move all terms to left
                faScalarMatrix MEqn
                (
                    rho_*h_
                   *(
                       fac::ddt(w_) - fac::ddt(w_.oldTime())
                    )/runTime().deltaT()
                  - fam::laplacian(M_) - p_
                );

                // Relax the linear system
                MEqn.relax();

                // Solve the linear system
                solverPerfM = MEqn.solve();

                // Relax the field
                M_.relax();

                // Store fields for under-relaxation and residual calculation
                w_.storePrevIter();

                // Solve w equation
                faScalarMatrix wEqn
                (
                    fam::laplacian(bendingStiffness_, w_) + M_
                );

                // Relax the linear system
                wEqn.relax();

                // Solve the linear system
                solverPerfw = wEqn.solve();

                // Relax the field
                w_.relax();

                // Update the angle of rotation
                theta_ = -fac::grad(w_);

                // Update the gradient of rotation field, used for the
                // non-orthogonal correction in clamped boundary conditions
                gradTheta_ = fac::grad(theta_);

                // Work in progress
                // Use block-coupled framework to implicitly couple the
                // equations
                //if (coupled_)
                // {
                //     // Prepare block system
                //     BlockLduSystem<vector2, vector2> blockM(aMesh_);

                //     // Grab block diagonal and set it to zero
                //     Field<tensor2>& d = blockM.diag().asSquare();
                //     d = tensor2::zero;

                //     // Grab linear off-diagonal and set it to zero
                //     Field<vector2>& l = blockM.lower().asLinear();
                //     Field<vector2>& u = blockM.upper().asLinear();
                //     u = vector2::zero;
                //     l = vector2::zero;

                //     // Insert MEqn and wEqn (without coupling terms)
                //     to-do

                //     // Insert coupling terms implicitly
                //     to-do

                //     // Solve the linear system
                //     to-do

                //     // Retrieve solution
                //     to-do

                //     // Correct the boundary conditions for M and w
                //     M.correctBoundaryConditions();
                //     w.correctBoundaryConditions();
                // }
            }
            while
            (
                !converged(iCorr, solverPerfM, solverPerfw, M_, w_)
             && ++iCorr < nCorr()
            );
        }

        // Map area fields to vol fields
        mapAreaFieldToSingleLayerVolumeField(M_, MVf_);
        mapAreaFieldToSingleLayerVolumeField(w_, wVf_);
        mapAreaFieldToSingleLayerVolumeField(theta_, thetaVf_);
        mapAreaFieldToSingleLayerVolumeField(p_, pVf_);
        {
            const areaVectorField Ds(w_*aMesh_.faceAreaNormals());
            mapAreaFieldToSingleLayerVolumeField(Ds, D());
        }

        // Interpolate cell displacements to vertices
        volToPoint().interpolate(D(), pointD());
        correctPointDisplacement(pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    // Disable writing of the stress fields, as it is not calculated
    sigma().writeOpt() = IOobject::NO_WRITE;

    return true;
}


tmp<vectorField> kirchhoffPlateSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    notImplemented(type() + "::tractionBoundarySnGrad(...)");

    // Keep compiler happy
    return tmp<vectorField>();
}


void kirchhoffPlateSolid::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    // Map global field to patch field
    const vectorField patchTraction
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTraction)
    );

    // Take normal component of the traction field
    // Note: p is the net pressure on the plate (from both sides)
#ifdef OPENFOAM_NOT_EXTEND
    p_.primitiveFieldRef() =
        aMesh_.faceAreaNormals().primitiveField() & patchTraction;
#else
    p_.internalField() =
        aMesh_.faceAreaNormals().internalField() & patchTraction;
#endif
}


void kirchhoffPlateSolid::writeFields(const Time& runTime)
{
    // Do not call solidModel::writeFields() as we do not want to write the
    // stress and strain fields

    physicsModel::writeFields(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // OPENFOAM_ORG

// ************************************************************************* //
