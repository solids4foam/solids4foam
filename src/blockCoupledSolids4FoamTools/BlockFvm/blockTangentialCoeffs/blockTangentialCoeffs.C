/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "blockTangentialCoeffs.H"

//#include "fv.H"
// To-do: figure out the minimum header files actually needed
#include "correctedSnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::fv::blockFvmInsertCoeffsTang
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& vf,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM,
    const newLeastSquaresVolPointInterpolation& volToPointInterp,
    const int op
)
{
    if (debug)
    {
        Info<< nl << "Adding tangential derivatve terms to matrix"
            << nl << "    start time: " << solidMesh().time().elapsedClockTime()
            << " s" << endl;
    }

    // Const reference to fvMesh
    const fvMesh& mesh = vf.mesh();

    // Grab block diagonal
    Field<tensor>& d = blockM.diag().asSquare();

    // Grab linear off-diagonal
    Field<tensor>& l = blockM.lower().asSquare();
    Field<tensor>& u = blockM.upper().asSquare();

    // Method
    // for all cells
    //     for all faces of cellI
    //         for all edges of faceI
    //             for all:
    //                 interpCells
    //                 interpBndFaces
    //                 interpCyclicFaces
    //                 interpProcFaces
    //                 glInterpNgbProcBndFaceCentres
    //                 glInterpNgbProcCellCentres
    //                 interpNgbProcCellCentres
    //                 interpNgbProcBndFaceCentres
    //             of edge.start() and edge.end()
    //                 add constributions to d, u and l
    //             end
    //         end
    //     end
    // end

    const surfaceVectorField n = mesh.Sf()/mesh.magSf();
    const unallocLabelList& fvOwner = mesh.owner();
    const scalarField& mufI = muf.internalField();
    const scalarField& lambdafI = lambdaf.internalField();
    const vectorField& points = mesh.points();
    const vectorField& nI = n.internalField();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelListList& pointCells = mesh.pointCells();
    const vectorField CI = mesh.C().internalField();

    const bool nonOrthogonalMesh = !mesh.orthogonal();
    const vectorField* corrVecIPtr = NULL;
    if (nonOrthogonalMesh)
    {
        corrVecIPtr = &mesh.correctionVectors().internalField();
    }

    List<vectorField> boundaryDelta(mesh.boundary().size());
    List<scalarField> boundaryDeltaCoeffs(mesh.boundary().size());
    forAll(boundaryDelta, patchI)
    {
        boundaryDelta[patchI] = mesh.boundary()[patchI].delta();
        boundaryDeltaCoeffs[patchI] = mesh.boundary()[patchI].deltaCoeffs();
    }

    // Least squares vol-to-point interpolation data
    // There is something not right when using least squares weights
    const vectorField& origins = volToPointInterp.origins();
    const FieldField<Field, scalar>& weights = volToPointInterp.weights();
    const PtrList<scalarRectangularMatrix>& invMatrices =
        volToPointInterp.invLsMatrices();
    const labelListList& pointBndFaces = volToPointInterp.pointBndFaces();
    const labelListList& pointCyclicFaces =
        volToPointInterp.pointCyclicFaces();
    const labelListList& pointProcFaces = volToPointInterp.pointProcFaces();
    Map<Field<vector> > gPtNgbProcBndFaceFieldData;
    volToPointInterp.globalPointNgbProcBndFaceFieldData
    (
        vf, gPtNgbProcBndFaceFieldData
    );
    Map<Field<vector> > gPtNgbProcCellFieldData;
    volToPointInterp.globalPointNgbProcCellFieldData
    (
        vf, gPtNgbProcCellFieldData
    );
    const Map< List<labelPair> >& pointProcCells =
        volToPointInterp.pointProcCells();
    const List< List<labelPair> >& pointProcBndFaces =
        volToPointInterp.pointProcBndFaces();
    const List<Tuple2<vector, tensor> >& mirrorPlaneTrans =
        volToPointInterp.mirrorPlaneTransformation();

    if (gMax(weights) > (1.0 + SMALL) || gMin(weights) < (1.0 - SMALL))
    {
        FatalErrorIn("insertTangCoeffs()")
            << "Weights currently assumed to be simple average" << nl
            << "Weights should be set to 1.0 in "
            << "newLeastSquaresVolPointInterpolation" << abort(FatalError);
    }

    // Global coefficients
    solidMesh.makeGlobalCoeffs(volToPointInterp);
    List< Map<tensorField> >& pointProcFacesCoeffs =
        solidMesh.pointProcFacesCoeffs();
    List< Map<tensorField> >& pointProcBndFacesCoeffs =
        solidMesh.pointProcBndFacesCoeffs();
    List< Map<tensorField> >& pointProcCellsCoeffs =
        solidMesh.pointProcCellsCoeffs();
    List< Map<tensorField> >& gPtNgbProcBndFaceCoeffs =
        solidMesh.gPtNgbProcBndFaceCoeffs();
    List< Map<tensorField> >& gPtNgbProcCellCoeffs =
        solidMesh.gPtNgbProcCellCoeffs();

    // Fixed points: value and direction
    const bool enforcePointConstraints = false;
    const Map<vector>& pointFixedComp = solidMesh.pointFixedComponent(vf);
    const Map<symmTensor>& pointFixedDir = solidMesh.pointFixedDirection(vf);


    // Perform discretisation

    forAll(cells, cellI)
    {
        if (debug > 1)
        {
            Info<< nl << "cell " << cellI << endl;
        }

        const labelList& curCellFaces = cells[cellI];

        forAll(curCellFaces, faceI)
        {
            const label curFaceID = curCellFaces[faceI];

            const bool intFace = mesh.isInternalFace(curFaceID);

            bool emptyFace = false;

            if (!mesh.isInternalFace(curFaceID))
            {
                const label patchID = mesh.boundaryMesh().whichPatch(curFaceID);

                if
                (
                    mesh.boundaryMesh()[patchID].type()
                == emptyPolyPatch::typeName
                )
                {
                    emptyFace = true;
                }
            }

            if (!emptyFace)
            {
                if (debug > 1)
                {
                    Info<< nl << "face " << curFaceID << " n " << nI[curFaceID]
                        << endl;
                }
                const face& curFace = faces[curFaceID];
                const edgeList curFaceEdges = curFace.edges();

                bool flipFace = false;
                scalar faceMu = 0.0;
                scalar faceLambda = 0.0;
                const vector* faceNPtr = NULL;
                const vector* faceKPtr = NULL;

                if (intFace)
                {
                    const label own = fvOwner[curFaceID];

                    faceMu = mufI[curFaceID];
                    faceLambda = lambdafI[curFaceID];

                    //const vector& faceN = nI[curFaceID];
                    faceNPtr = &(nI[curFaceID]);

                    // Check if cellI is the owner
                    if (cellI != own)
                    {
                        flipFace = true;
                    }

                    // Non-orthogonal correction vector
                    if (nonOrthogonalMesh)
                    {
                        faceKPtr = &((*corrVecIPtr)[curFaceID]);
                    }
                }
                else
                {
                    label patchID = mesh.boundaryMesh().whichPatch(curFaceID);
                    label start = mesh.boundaryMesh()[patchID].start();

                    const label patchFaceID = curFaceID - start;
                    faceMu = muf.boundaryField()[patchID][patchFaceID];
                    faceLambda =
                        lambdaf.boundaryField()[patchID][patchFaceID];

                    faceNPtr = &(n.boundaryField()[patchID][patchFaceID]);

                    if (nonOrthogonalMesh)
                    {
                        const vector& faceN = *faceNPtr;
                        const vector& faceDelta =
                            boundaryDelta[patchID][patchFaceID];
                        const scalar faceDeltaCoeff =
                            boundaryDeltaCoeffs[patchID][patchFaceID];
                        faceKPtr =
                            new vector
                            (
                                pos(faceN & faceDelta)
                               *(faceN - faceDelta*faceDeltaCoeff)
                            );
                    }
                }
                const vector& faceN = *faceNPtr;

                if (faceKPtr == NULL)
                {
                    faceKPtr = new vector(0, 0, 0);
                }
                const vector& faceK = *faceKPtr;

                // Insert coefficients for this face
                blockFvmInsertCoeffsTangForFace
                (
                    d,
                    u,
                    l,
                    blockB[cellI],
                    pointProcFacesCoeffs[cellI],
                    pointProcBndFacesCoeffs[cellI],
                    pointProcCellsCoeffs[cellI],
                    gPtNgbProcBndFaceCoeffs[cellI],
                    gPtNgbProcCellCoeffs[cellI],
                    cellI,
                    faces[curFaceID],
                    flipFace,
                    nonOrthogonalMesh,
                    faceMu,
                    faceLambda,
                    faceN,
                    faceK,
                    points,
                    origins,
                    invMatrices,
                    weights,
                    enforcePointConstraints,
                    pointFixedComp,
                    pointFixedDir,
                    pointCells,
                    mirrorPlaneTrans,
                    solidMesh,
                    pointBndFaces,
                    pointCyclicFaces,
                    pointProcFaces,
                    gPtNgbProcBndFaceFieldData,
                    gPtNgbProcCellFieldData,
                    pointProcCells,
                    pointProcBndFaces,
                    op
                );
            } // if not an empty face
        } // forAll faces of current cell

        if (debug > 1)
        {
            Info<< nl << "    total d[" << cellI << "] " << d[cellI] << nl
                << endl;
        }
    } // forAll cells

    if (fv::debug)
    {
        Info<< "    end time: " << solidMesh().time().elapsedClockTime() << " s"
            << nl << endl;
    }
}


void Foam::fv::blockFvmInsertCoeffsTangForFace
(
    tensorField& d,
    tensorField& u,
    tensorField& l,
    vector& blockB,
    Map<tensorField>& pointProcFacesCoeffs,
    Map<tensorField>& pointProcBndFacesCoeffs,
    Map<tensorField>& pointProcCellsCoeffs,
    Map<tensorField>& gPtNgbProcBndFaceCoeffs,
    Map<tensorField>& gPtNgbProcCellCoeffs,
    const label cellI,
    const face& curFace,
    const bool flipFace,
    const bool nonOrthogonalMesh,
    const scalar faceMu,
    const scalar faceLambda,
    const vector& faceN,
    const vector& faceK,
    const vectorField& points,
    const vectorField& origins,
    const PtrList<scalarRectangularMatrix>& invMatrices,
    const FieldField<Field, scalar>& weights,
    const bool enforcePointConstraints,
    const Map<vector>& pointFixedComp,
    const Map<symmTensor>& pointFixedDir,
    const labelListList& pointCells,
    const List< Tuple2<vector, tensor> >& mirrorPlaneTrans,
    const solidPolyMesh& solidMesh,
    const labelListList& pointBndFaces,
    const labelListList& pointCyclicFaces,
    const labelListList& pointProcFaces,
    const Map<Field<vector> >& gPtNgbProcBndFaceFieldData,
    const Map<Field<vector> >& gPtNgbProcCellFieldData,
    const Map< List<labelPair> >& pointProcCells,
    const List< List<labelPair> >& pointProcBndFaces,
    const int op
)
{
    const edgeList curFaceEdges = curFace.edges();

    forAll(curFaceEdges, edgeI)
    {
        if (debug > 1)
        {
            Info<< nl << "edge " << curFaceEdges[edgeI].centre(points) << endl;
        }

        const edge& curEdge = curFaceEdges[edgeI];

        // Calculate Le*faceN which features in all the coefficients
        // Projected edge vector
        vector e = curEdge.vec(points);
        e -= faceN*(faceN & e);

        // Edge length vector
        vector Le = e^faceN;
        Le *= curFace.edgeDirection(curEdge);

        // LeFaceN term is the tensor component of the coeff
        // Change sign if cellI is not the owner as faceN is
        // opposite
        tensor LeFaceN = Le*faceN;
        if (flipFace)
        {
            LeFaceN *= -1;
        }

        // Non-orthogonal correction component
        scalar kDotLe = 0.0;
        if (nonOrthogonalMesh)
        {
            kDotLe = (faceK & Le);

            if (flipFace)
            {
                kDotLe *= -1;
            }
        }

        //const label curEdgeID = faceEdges[curFaceID][edgeI];
        const label startPointID = curEdge.start();
        const label endPointID = curEdge.end();

        // We add coeffs to:
        //     interpCells
        //     interpBndFaces
        //     interpCyclicFaces
        //     interpProcFaces
        //     glInterpNgbProcBndFaceCentres
        //     glInterpNgbProcCellCentres
        //     interpNgbProcCellCentres
        //     interpNgbProcBndFaceCentres

        // Put start and end point ID in a list so we will perform
        // the same actions to both start and end vertex
        labelList sePointIDs(2, -1);
        sePointIDs[0] = startPointID;
        sePointIDs[1] = endPointID;

        // We will add a coeff contribution to all neighbour
        // computational nodes

        forAll(sePointIDs, pointI)
        {
            // First add coeffs to pointCells on current proc,
            // similar to method for internal edges.

            // Start/end point ID
            const label sePointID = sePointIDs[pointI];

            // ID for neighbour interpolation point
            label neiID = 0;

            // Vector between current point and the average
            // position of the neighbours
            vector dr = points[sePointID] - origins[sePointID];

            // Least square inverse matrix terms
            const scalarRectangularMatrix& curInvMatrix =
                invMatrices[sePointID];

            // Weights for the current point
            const scalarField& w = weights[sePointID];

            // Calculate sum of weights but we do not include mirror
            // points
            //const scalar sumW = sum(sqr(w));
            const scalar sumW = sum(w);

            // Check if the point has any fixed components
            bool pointHasFixedComp = false;
            symmTensor sePointFixedDir = symmTensor::zero;
            vector sePointFixedComp = vector::zero;
            if (enforcePointConstraints)
            {
                pointHasFixedComp = pointFixedComp.found(sePointID);

                if (pointHasFixedComp)
                {
                    sePointFixedDir = pointFixedDir[sePointID];
                    sePointFixedComp = pointFixedComp[sePointID];
                }
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~         Add coeffs for local point cells         ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


            const labelList& sePointCells = pointCells[sePointID];

            for (label i = 0; i < sePointCells.size(); i++)
            {
                // Least squares weights
                const vector m =
                    vector
                    (
                        curInvMatrix[0][neiID],
                        curInvMatrix[1][neiID],
                        curInvMatrix[2][neiID]
                    );

                // Calculate the coefficient
                tensor coeff = tensor::zero;
                blockFvmCalculateCoeff
                (
                    coeff,
                    blockB,
                    w[neiID]/sumW,
                    m,
                    dr,
                    faceMu,
                    faceLambda,
                    faceN,
                    LeFaceN,
                    kDotLe,
                    nonOrthogonalMesh,
                    mirrorPlaneTrans[sePointID],
                    pointHasFixedComp,
                    sePointFixedComp,
                    sePointFixedDir,
                    op
                );

                // Add coeff contribution to cellI from
                // pointCell
                const label sePointCellI = sePointCells[i];

                if (sePointCellI == cellI)
                {
                    // Contributions to diag

                    if (debug > 1)
                    {
                        Info<< "    d[" << cellI << "] p " << coeff << endl;
                    }
                    d[cellI] += coeff;
                }
                else
                {
                    // Contributions to upper/lower

                    // Find which implicit bond connects cellI
                    // to pointCellI
                    const label curImplicitBondID =
                        solidMesh.findCellCellImplicitBond(cellI, sePointCellI);

                    if (curImplicitBondID == -1)
                    {
                        FatalErrorIn
                        (
                            "pointGaussLsDivSigmaScheme::insertCoeffBc()"
                        )   << "curImplicitBond not found"
                            << abort(FatalError);
                    }

                    if (cellI > sePointCellI)
                    {
                        if (debug > 1)
                        {
                            Info<< "    l[" << curImplicitBondID << "] p "
                                << coeff << endl;
                        }
                        l[curImplicitBondID] += coeff;
                    }
                    else
                    {
                        if (debug > 1)
                        {
                            Info<< "    u[" << curImplicitBondID << "] p "
                                << coeff << endl;
                        }
                        u[curImplicitBondID] += coeff;
                    }
                }

                neiID++;
            } // forAll pointCells


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~        Add coeffs for local boundary faces       ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            const labelList& sePointBndFaces = pointBndFaces[sePointID];

            for (label i = 0; i < sePointBndFaces.size(); i++)
            {
                // Least squares weights
                const vector m =
                    vector
                    (
                        curInvMatrix[0][neiID],
                        curInvMatrix[1][neiID],
                        curInvMatrix[2][neiID]
                    );

                // Calculate the coefficient
                tensor coeff = tensor::zero;
                blockFvmCalculateCoeff
                (
                    coeff,
                    blockB,
                    w[neiID]/sumW,
                    m,
                    dr,
                    faceMu,
                    faceLambda,
                    faceN,
                    LeFaceN,
                    kDotLe,
                    nonOrthogonalMesh,
                    mirrorPlaneTrans[sePointID],
                    pointHasFixedComp,
                    sePointFixedComp,
                    sePointFixedDir,
                    op
                );

                // Add coeff contribution to the upper of cellI from
                // boundary the face

                const label sePointBndFaceI = sePointBndFaces[i];

                // Find which implicit bond connects cellI to
                // pointCellI
                const label curImplicitBondID =
                    solidMesh.findCellFaceImplicitBond(cellI, sePointBndFaceI);

                if (curImplicitBondID == -1)
                {
                    FatalErrorIn
                    (
                        "pointGaussLsDivSigmaScheme::insertCoeffBc()"
                    )   << "bndFace: curImplicitBond not found"
                        << abort(FatalError);
                }

                // Contributions to upper
                if (debug > 1)
                {
                    Info<< "    l[" << curImplicitBondID << "] p "
                        << coeff << endl;
                }
                u[curImplicitBondID] += coeff;

                neiID++;
            } // for point boundary faces


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~        Add coeffs for cyclic boundary faces      ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            const labelList& sePointCyclicFaces =
                pointCyclicFaces[sePointID];

            if (sePointCyclicFaces.size())
            {
                FatalErrorIn
                (
                    "void insertCoeffsTang\n"
                    "(\n"
                    "    const solidPolyMesh& solidMesh,\n"
                    "    const surfaceScalarField& muf,\n"
                    "    const surfaceScalarField& lambdaf,\n"
                    "    GeometricField<vector, fvPatchField, "
                    "volMesh>& vf,\n"
                    "    Field<vector>& blockB,\n"
                    "    BlockLduMatrix<vector>& blockM\n"
                    ")\n"
                )   << "cyclic faces not implemented"
                    << abort(FatalError);
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~        Add coeffs for proc boundary faces        ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            const labelList& sePointProcFaces =
                pointProcFaces[sePointID];

            for (label i = 0; i < sePointProcFaces.size(); i++)
            {
                // Least squares weights
                const vector m =
                    vector
                    (
                        curInvMatrix[0][neiID],
                        curInvMatrix[1][neiID],
                        curInvMatrix[2][neiID]
                    );

                // Calculate the coefficient
                tensor coeff = tensor::zero;
                blockFvmCalculateCoeff
                (
                    coeff,
                    blockB,
                    w[neiID]/sumW,
                    m,
                    dr,
                    faceMu,
                    faceLambda,
                    faceN,
                    LeFaceN,
                    kDotLe,
                    nonOrthogonalMesh,
                    mirrorPlaneTrans[sePointID],
                    pointHasFixedComp,
                    sePointFixedComp,
                    sePointFixedDir,
                    op
                );

                // Add coeff contribution to globalCoeff
                pointProcFacesCoeffs[sePointID][i] += coeff;

                neiID++;
            } // forAll proc boundary faces


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~        Add coeffs for global boundary faces      ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            // These are boundary faces on another proc which are
            // adjoining a global point

            if (gPtNgbProcBndFaceFieldData.found(sePointID))
            {
                const Field<vector> glNgbProcBndFaces =
                    gPtNgbProcBndFaceFieldData[sePointID];

                for (label i = 0; i < glNgbProcBndFaces.size(); i++)
                {
                    // Least squares weights
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate the coefficient
                    tensor coeff = tensor::zero;
                    blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB,
                        w[neiID]/sumW,
                        m,
                        dr,
                        faceMu,
                        faceLambda,
                        faceN,
                        LeFaceN,
                        kDotLe,
                        nonOrthogonalMesh,
                        mirrorPlaneTrans[sePointID],
                        pointHasFixedComp,
                        sePointFixedComp,
                        sePointFixedDir,
                        op
                    );

                    // Add coeff contribution to globalCoeff
                    gPtNgbProcBndFaceCoeffs[sePointID][i] +=
                        coeff;

                    neiID++;
                } // forAll proc boundary faces
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~           Add coeffs for global cells            ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            // These are cells on another proc which are adjoining a
            // global point

            if (gPtNgbProcCellFieldData.found(sePointID))
            {
                // const List<labelPair>& glNgbProcCells =
                const Field<vector> glNgbProcCells =
                    gPtNgbProcCellFieldData[sePointID];

                for (label i = 0; i < glNgbProcCells.size(); i++)
                {
                    // Least squares weights
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate the coefficient
                    tensor coeff = tensor::zero;
                    blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB,
                        w[neiID]/sumW,
                        m,
                        dr,
                        faceMu,
                        faceLambda,
                        faceN,
                        LeFaceN,
                        kDotLe,
                        nonOrthogonalMesh,
                        mirrorPlaneTrans[sePointID],
                        pointHasFixedComp,
                        sePointFixedComp,
                        sePointFixedDir,
                        op
                    );

                    // Add coeff contribution to globalCoeff
                    gPtNgbProcCellCoeffs[sePointID][i] += coeff;

                    neiID++;
                } // forAll proc boundary faces
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~             Add coeffs for proc cells            ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            if (pointProcCells.found(sePointID))
            {
                const List<labelPair>& sePointProcCells =
                    pointProcCells[sePointID];

                for (label i = 0; i < sePointProcCells.size(); i++)
                {
                    // Least squares weights
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate the coefficient
                    tensor coeff = tensor::zero;
                    blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB,
                        w[neiID]/sumW,
                        m,
                        dr,
                        faceMu,
                        faceLambda,
                        faceN,
                        LeFaceN,
                        kDotLe,
                        nonOrthogonalMesh,
                        mirrorPlaneTrans[sePointID],
                        pointHasFixedComp,
                        sePointFixedComp,
                        sePointFixedDir,
                        op
                    );

                    // Add coeff contribution to globalCoeff
                    pointProcCellsCoeffs[sePointID][i] += coeff;

                    neiID++;
                } // forAll proc boundary faces
            }


            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
            //~   Add coeffs for boundary faces on other procs   ~//
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

            if (pointProcBndFaces[sePointID].size())
            {
                const List<labelPair>& sePointProcBndFaces =
                    pointProcBndFaces[sePointID];

                for (label i = 0; i < sePointProcBndFaces.size(); i++)
                {
                    // Least squares weights
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate the coefficient
                    tensor coeff = tensor::zero;
                    blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB,
                        w[neiID]/sumW,
                        m,
                        dr,
                        faceMu,
                        faceLambda,
                        faceN,
                        LeFaceN,
                        kDotLe,
                        nonOrthogonalMesh,
                        mirrorPlaneTrans[sePointID],
                        pointHasFixedComp,
                        sePointFixedComp,
                        sePointFixedDir,
                        op
                    );

                    // Add coeff contribution to globalCoeff
                    pointProcBndFacesCoeffs[sePointID][i] += coeff;

                    neiID++;
                } // forAll proc boundary faces
            }
        } // for the edge start and end points
    } // forAll edges of current face
}


void Foam::fv::blockFvmCalculateCoeff
(
    tensor& coeff,
    vector& blockB,
    const scalar wCell,
    const vector& m,
    const vector& dr,
    const scalar faceMu,
    const scalar faceLambda,
    const vector& faceN,
    const tensor& LeFaceN,
    const scalar kDotLe,
    const bool nonOrthogonalMesh,
    const Tuple2<vector, tensor>& mirrorPlaneTrans,
    const bool pointHasFixedComp,
    const vector& pointFixedComp,
    const symmTensor& pointFixedDir,
    const int op
)
{
    // Op decides which operator is discretised:
    // 0 == Laplacian
    // 1 == Laplacian Transpose
    // 2 == Laplacian Trace
    // 3 == Laplacian + Transpose + Trace

    // Calculate the weight to interpolate from the current point to the
    // neighbour cell/face
    // wCell comes from simple averaging and the (m & dr) is a correction based
    // on a least squares fit
    const scalar weight = wCell + (m & dr);

    if (op == 0 || op == 3)
    {
        // Laplacian
        if (nonOrthogonalMesh)
        {
            coeff += tensor(0.5*I*faceMu*kDotLe)*weight;
        }
    }

    if (op == 1 || op == 3)
    {
        // Laplacian Transpose
        coeff += 0.5*(faceMu*LeFaceN)*weight;

        if (nonOrthogonalMesh)
        {
            coeff += tensor(0.5*sqr(faceN)*faceMu*kDotLe)*weight;
        }
    }

    if (op == 2 || op == 3)
    {
        // Laplacian Trace
        coeff += 0.5*(faceLambda*LeFaceN.T())*weight;

        if (nonOrthogonalMesh)
        {
            coeff += tensor(0.5*sqr(faceN)*faceLambda*kDotLe)*weight;
        }
    }

    if (op != 0 && op != 1 && op != 2 && op != 3)
    {
        FatalErrorIn("blockFvmCalculateCoeff()")
            << "op must be 0, 1, 2 or 3" << abort(FatalError);
    }

    // Add contribution from mirrored cell/face, if any
    if (mag(mirrorPlaneTrans.first()) > SMALL)
    {
        const tensor& T = mirrorPlaneTrans.second();

        coeff += transform(T, coeff);
    }

    // Check if the point has a fixed component
    if (pointHasFixedComp)
    {
        FatalErrorIn("void Foam::fv::blockFvmCalculateCoeff(...)")
            << "Enforcing fixed points is disabled!" << abort(FatalError);

        // Add explicitly to the source
        blockB -= pointFixedComp & coeff;

        // Remove coeff in fixed direction
        coeff = ((I - pointFixedDir) & coeff);
    }
}


// ************************************************************************* //
