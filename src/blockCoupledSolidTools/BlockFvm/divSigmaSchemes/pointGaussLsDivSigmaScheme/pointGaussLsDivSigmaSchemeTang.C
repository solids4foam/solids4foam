/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pointGaussLsDivSigmaScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //


void pointGaussLsDivSigmaScheme::insertCoeffsTang
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& vf,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
)
{
    if (debug)
    {
        Info<< nl << "Adding tangential derivatve terms to matrix" << endl;
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
    // WarningIn("insertTangCoeff()")
    //     << "Internal non-orthogonal correction set to zero"
    //     << endl;
    // const bool nonOrthogonalMesh = false;
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
    const vectorField& origins = volToPointInterp().origins();
    const FieldField<Field, scalar>& weights = volToPointInterp().weights();
    const PtrList<scalarRectangularMatrix>& invMatrices =
        volToPointInterp().invLsMatrices();
    const labelListList& pointBndFaces = volToPointInterp().pointBndFaces();
    const labelListList& pointCyclicFaces =
        volToPointInterp().pointCyclicFaces();
    const labelListList& pointProcFaces = volToPointInterp().pointProcFaces();
    Map<Field<vector> > gPtNgbProcBndFaceFieldData;
    volToPointInterp().globalPointNgbProcBndFaceFieldData
    (
        vf, gPtNgbProcBndFaceFieldData
    );
    Map<Field<vector> > gPtNgbProcCellFieldData;
    volToPointInterp().globalPointNgbProcCellFieldData
    (
        vf, gPtNgbProcCellFieldData
    );
    const Map< List<labelPair> >& pointProcCells =
        volToPointInterp().pointProcCells();
    const List< List<labelPair> >& pointProcBndFaces =
        volToPointInterp().pointProcBndFaces();
    const List<Tuple2<vector, tensor> >& mirrorPlaneTrans =
        volToPointInterp().mirrorPlaneTransformation();

    if (gMax(weights) > (1.0 + SMALL) || gMin(weights) < (1.0 - SMALL))
    {
        FatalErrorIn("insertTangCoeffs()")
            << "Weights currently assumed to be simple average" << nl
            << "Weights should be set to 1.0 in "
            << "newLeastSquaresVolPointInterpolation" << abort(FatalError);
    }

    // Global coefficients
    solidMesh.makeGlobalCoeffs(volToPointInterp());
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
    const Map<vector>& pointFixedComp = solidMesh.pointFixedComponent(vf);
    const Map<symmTensor>& pointFixedDir = solidMesh.pointFixedDirection(vf);


    // Perform discretisation

    if (debug)
    {
        Info<< "    Adding tangential coeffs for each cell" << endl;
    }

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
                label patchID = mesh.boundaryMesh().whichPatch(curFaceID);

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
                        //faceKPtr = new vector((I - sqr(faceN)) & faceDelta);
                        faceKPtr =
                            new vector
                            (
                                pos(faceN & faceDelta)
                               *(faceN - faceDelta*faceDeltaCoeff)
                            );
                    }
                }
                const vector& faceN = *faceNPtr;

                // Store (2*mu + lambda) as it is used for each edge
                const scalar faceTwoMuLambda = 2*faceMu + faceLambda;

                forAll(curFaceEdges, edgeI)
                {
                    if (debug > 1)
                    {
                        Info<< nl << "edge "
                            << curFaceEdges[edgeI].centre(mesh.points())
                            << endl;
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
                    scalar* kDotLePtr = NULL;
                    if (nonOrthogonalMesh)
                    {
                        kDotLePtr = new scalar((*faceKPtr) & Le);

                        if (flipFace)
                        {
                            (*kDotLePtr) *= -1;
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
                        //vector dr = vector::zero;

                        // Least square inverse matrix terms
                        const scalarRectangularMatrix& curInvMatrix =
                            invMatrices[sePointID];

                        // Weights for the current point
                        const scalarField& w = weights[sePointID];

                        // Calculate sumof weights but we do not include mirror
                        // points
                        //const scalar sumW = sum(sqr(w));
                        const scalar sumW = sum(w);

                        // Sum least squares matrix weights
                        // Hmmnn sum(m) is zero or very small...
                        // vector sumM = vector::zero;
                        // const label nNei = w.size();
                        // for (label coeffI = 0; coeffI < 3; coeffI++)
                        // {
                        //     for (label neiI = 0; neiI < nNei; neiI++)
                        //     {
                        //         sumM[coeffI] += curInvMatrix[coeffI][neiI];
                        //     }
                        // }

                        // Check if the point has any fixed components
                        const bool pointHasFixedComp =
                            pointFixedComp.found(sePointID);
                        symmTensor sePointFixedDir = symmTensor::zero;
                        if (pointHasFixedComp)
                        {
                            sePointFixedDir = pointFixedDir[sePointID];
                        }


                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                        //~         Add coeffs for local point cells         ~//
                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


                        const labelList& sePointCells = pointCells[sePointID];

                        for (label i = 0; i < sePointCells.size(); i++)
                        {
                            scalar wCell = w[neiID]/sumW;

                            // Coefficient is composed of the contribution from
                            // pure averaging (inverse distance) and a
                            // correction from least squares

                            const vector m
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                            // Contribution from pure inverse distance average
                            // (wCell) and a correction from least squares due
                            // to a nonuniform distribution of points (m & dr)
                            // So least squares correction is zero if the vertex
                            // is located at the average position of the
                            // neighbours
                            // Note: LeFaceN.T()
                            //     == (LeFaceN.T() & (I - sqr(faceN)))
                            // because n & LeFaceN == 0
                            tensor coeff =
                                //0.5*(wCell + ((m - sumM) & dr))
                                0.5*(wCell + (m & dr))
                               *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                            // if (sePointID == 1210)
                            // {
                            //     Info<< "wCell " << wCell << nl
                            //         << "m & dr " << (m & dr) << nl
                            //         << "sumM & dr " << (sumM & dr) << endl;
                            // }

                            // Add non-orthogonal correction term
                            if (nonOrthogonalMesh)
                            {
                                coeff +=
                                    tensor
                                    (
                                        0.5*(wCell + (m & dr))*
                                        (
                                            sqr(faceN)*
                                            faceTwoMuLambda*(*kDotLePtr)
                                            + (I - sqr(faceN))*
                                            faceMu*(*kDotLePtr)
                                        )
                                    );
                            }

                            // Add contribution from mirrored cell, if any
                            if
                            (
                                mag(mirrorPlaneTrans[sePointID].first()) > SMALL
                            )
                            {
                                const tensor& T =
                                    mirrorPlaneTrans[sePointID].second();

                                coeff += transform(T, coeff);
                            }

                            // Check if the point has a fixed component
                            if (pointHasFixedComp)
                            {
                                // Remove coeff in fixed direction
                                coeff = (coeff & (I - sePointFixedDir));

                                // Add explicitly to the source: todo below
                                //blockB += coeff & pointFixedComp(sePointID);
                            }

                            // Add coeff contribution to cellI from
                            // pointCell
                            const label sePointCellI = sePointCells[i];

                            if (sePointCellI == cellI)
                            {
                                // Contributions to diag

                                if (debug > 1)
                                {
                                    Info<< "    d[" << cellI << "] p "
                                        << coeff << endl;
                                }
                                d[cellI] += coeff;
                            }
                            else
                            {
                                // Contributions to upper/lower

                                // Find which implicit bond connects cellI
                                // to pointCellI
                                const label curImplicitBondID =
                                    solidMesh.findCellCellImplicitBond
                                    (
                                        cellI,
                                        sePointCellI
                                    );

                                if (curImplicitBondID == -1)
                                {
                                    FatalErrorIn
                                    (
                                        "pointGaussLsDivSigmaScheme::"
                                        "insertCoeffBc()"
                                    )   << "curImplicitBond not found"
                                        << abort(FatalError);
                                }

                                if (cellI > sePointCellI)
                                {
                                    if (debug > 1)
                                    {
                                        Info<< "    l[" << curImplicitBondID
                                            << "] p "
                                            << coeff << endl;
                                    }
                                    l[curImplicitBondID] += coeff;
                                }
                                else
                                {
                                    if (debug > 1)
                                    {
                                        Info<< "    u[" << curImplicitBondID
                                            << "] p "
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

                        const labelList& sePointBndFaces =
                            pointBndFaces[sePointID];

                        for(label i = 0; i < sePointBndFaces.size(); i++)
                        {
                            scalar wCell = w[neiID]/sumW;

                            // Coefficient is composed of the contribution from
                            // pure averaging (inverse distance) and a
                            // correction from least squares

                            const vector m
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                            // Contribution from pure inverse distance average
                            // (wCell) and a correction from least squares due
                            // to a nonuniform distribution of points (m & dr)
                            // So least squares correction is zero if the vertex
                            // is located at the average position of the
                            // neighbours
                            // Note: LeFaceN.T()
                            //     == (LeFaceN.T() & (I - sqr(faceN)))
                            // because n & LeFaceN == 0
                            tensor coeff =
                                0.5*(wCell + (m & dr))
                               *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                            // Add non-orthogonal correction term
                            if (nonOrthogonalMesh)
                            {
                                coeff +=
                                    tensor
                                    (
                                        0.5*(wCell + (m & dr))*
                                        (
                                            sqr(faceN)*
                                            faceTwoMuLambda*(*kDotLePtr)
                                          + (I - sqr(faceN))*faceMu*(*kDotLePtr)
                                        )
                                    );
                            }

                            // Add contribution from mirrored cell, if any
                            if
                            (
                                mag(mirrorPlaneTrans[sePointID].first()) > SMALL
                            )
                            {
                                const tensor& T =
                                    mirrorPlaneTrans[sePointID].second();

                                coeff += transform(T, coeff);
                            }

                            // Check if the point has a fixed component
                            if (pointHasFixedComp)
                            {
                                // Remove coeff in fixed direction
                                coeff = (coeff & (I - sePointFixedDir));

                                // Add explicitly to the source: todo below
                                //blockB += coeff & pointFixedComp(sePointID);
                            }

                            // Add coeff contribution to the upper of cellI from
                            // boundary the face

                            const label sePointBndFaceI = sePointBndFaces[i];

                            // Find which implicit bond connects cellI to
                            // pointCellI
                            const label curImplicitBondID =
                                solidMesh.findCellFaceImplicitBond
                                (
                                    cellI,
                                    sePointBndFaceI
                                );

                            if (curImplicitBondID == -1)
                            {
                                FatalErrorIn
                                (
                                    "pointGaussLsDivSigmaScheme::"
                                    "insertCoeffBc()"
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
                            scalar wCell = w[neiID]/sumW;

                            // Coefficient is composed of the contribution from
                            // pure averaging (inverse distance) and a
                            // correction from least squares

                            vector m
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                            // Contribution from pure inverse distance average
                            // (wCell) and a correction from least squares due
                            // to a nonuniform distribution of points (m & dr)
                            // So least squares correction is zero if the vertex
                            // is located at the average position of the
                            // neighbours
                            // Note: LeFaceN.T()
                            //     == (LeFaceN.T() & (I - sqr(faceN)))
                            // because n & LeFaceN == 0
                            tensor coeff =
                                0.5*(wCell + (m & dr))
                               *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                            // Add non-orthogonal correction term
                            if (nonOrthogonalMesh)
                            {
                                coeff +=
                                    tensor
                                    (
                                        0.5*(wCell + (m & dr))*
                                        (
                                            sqr(faceN)*
                                            faceTwoMuLambda*(*kDotLePtr)
                                            + (I - sqr(faceN))*
                                            faceMu*(*kDotLePtr)
                                        )
                                    );
                            }

                            // Add contribution from mirrored cell, if any
                            if
                            (
                                mag(mirrorPlaneTrans[sePointID].first()) > SMALL
                            )
                            {
                                const tensor& T =
                                    mirrorPlaneTrans[sePointID].second();

                                coeff += transform(T, coeff);
                            }

                            // Add coeff contribution to globalCoeff
                            pointProcFacesCoeffs[cellI][sePointID][i] += coeff;

                            neiID++;
                        } // forAll proc boundary faces


                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                        //~        Add coeffs for global boundary faces      ~//
                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                        // These are boundary faces on another proc which are
                        // adjoining a global point

                        if (gPtNgbProcBndFaceFieldData.found(sePointID))
                        {
                            // const List<labelPair>& glNgbProcBndFaces =
                            const Field<vector> glNgbProcBndFaces =
                                gPtNgbProcBndFaceFieldData[sePointID];

                            for (label i = 0; i < glNgbProcBndFaces.size(); i++)
                            {
                                scalar wCell = w[neiID]/sumW;

                                // Contribution from pure inverse distance
                                // average (wCell) and a correction from least
                                // squares due to a nonuniform distribution of
                                // points (m & dr) so least squares correction
                                // is zero if the vertex is located at the
                                // average position of the neighbours

                                vector m
                                (
                                    curInvMatrix[0][neiID],
                                    curInvMatrix[1][neiID],
                                    curInvMatrix[2][neiID]
                                );

                                tensor coeff =
                                    0.5*(wCell + (m & dr))
                                    *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                                // Add non-orthogonal correction term
                                if (nonOrthogonalMesh)
                                {
                                    coeff +=
                                        tensor
                                        (
                                            0.5*(wCell + (m & dr))*
                                            (
                                                sqr(faceN)*
                                                faceTwoMuLambda*(*kDotLePtr)
                                                + (I - sqr(faceN))*
                                                faceMu*(*kDotLePtr)
                                            )
                                        );
                                }

                                // Add contribution from mirrored cell, if any
                                if
                                (
                                    mag(mirrorPlaneTrans[sePointID].first())
                                    > SMALL
                                )
                                {
                                    const tensor& T =
                                        mirrorPlaneTrans[sePointID].second();

                                    coeff += transform(T, coeff);
                                }

                                // Add coeff contribution to globalCoeff
                                gPtNgbProcBndFaceCoeffs[cellI][sePointID][i] +=
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
                                scalar wCell = w[neiID]/sumW;

                                // Contribution from pure inverse distance
                                // average (wCell) and a correction from least
                                // squares due to a nonuniform distribution of
                                // points (m & dr) so least squares correction
                                // is zero if the vertex is located at the
                                // average position of the neighbours

                                vector m
                                (
                                    curInvMatrix[0][neiID],
                                    curInvMatrix[1][neiID],
                                    curInvMatrix[2][neiID]
                                );

                                tensor coeff =
                                    0.5*(wCell + (m & dr))
                                    *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                                // Add non-orthogonal correction term
                                if (nonOrthogonalMesh)
                                {
                                    coeff +=
                                        tensor
                                        (
                                            0.5*(wCell + (m & dr))*
                                            (
                                                sqr(faceN)*
                                                faceTwoMuLambda*(*kDotLePtr)
                                                + (I - sqr(faceN))*
                                                faceMu*(*kDotLePtr)
                                            )
                                        );
                                }

                                // Add contribution from mirrored cell, if any
                                if
                                (
                                    mag(mirrorPlaneTrans[sePointID].first())
                                    > SMALL
                                )
                                {
                                    const tensor& T =
                                        mirrorPlaneTrans[sePointID].second();

                                    coeff += transform(T, coeff);
                                }

                                // Add coeff contribution to globalCoeff
                                gPtNgbProcCellCoeffs[cellI][sePointID][i] +=
                                    coeff;

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
                                scalar wCell = w[neiID]/sumW;

                                // Contribution from pure inverse distance
                                // average (wCell) and a correction from least
                                // squares due to a nonuniform distribution of
                                // points (m & dr) so least squares correction
                                // is zero if the vertex is located at the
                                // average position of the neighbours

                                vector m
                                (
                                    curInvMatrix[0][neiID],
                                    curInvMatrix[1][neiID],
                                    curInvMatrix[2][neiID]
                                );

                                tensor coeff =
                                    0.5*(wCell + (m & dr))
                                    *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                                // Add non-orthogonal correction term
                                if (nonOrthogonalMesh)
                                {
                                    coeff +=
                                        tensor
                                        (
                                            0.5*(wCell + (m & dr))*
                                            (
                                                sqr(faceN)*
                                                faceTwoMuLambda*(*kDotLePtr)
                                                + (I - sqr(faceN))*
                                                faceMu*(*kDotLePtr)
                                            )
                                        );
                                }

                                // Add contribution from mirrored cell, if any
                                if
                                (
                                    mag(mirrorPlaneTrans[sePointID].first())
                                    > SMALL
                                )
                                {
                                    const tensor& T =
                                        mirrorPlaneTrans[sePointID].second();

                                    coeff += transform(T, coeff);
                                }

                                // Add coeff contribution to globalCoeff
                                pointProcCellsCoeffs[cellI][sePointID][i] +=
                                    coeff;

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

                            for
                            (
                                label i = 0; i < sePointProcBndFaces.size(); i++
                            )
                            {
                                scalar wCell = w[neiID]/sumW;

                                // Contribution from pure inverse distance
                                // average (wCell) and a correction from least
                                // squares due to a nonuniform distribution of
                                // points (m & dr) so least squares correction
                                // is zero if the vertex is located at the
                                // average position of the neighbours

                                vector m
                                (
                                    curInvMatrix[0][neiID],
                                    curInvMatrix[1][neiID],
                                    curInvMatrix[2][neiID]
                                );

                                tensor coeff =
                                    0.5*(wCell + (m & dr))
                                    *(faceMu*LeFaceN + faceLambda*LeFaceN.T());

                                // Add non-orthogonal correction term
                                if (nonOrthogonalMesh)
                                {
                                    coeff +=
                                        tensor
                                        (
                                            0.5*(wCell + (m & dr))*
                                            (
                                                sqr(faceN)*
                                                faceTwoMuLambda*(*kDotLePtr)
                                                + (I - sqr(faceN))*
                                                faceMu*(*kDotLePtr)
                                            )
                                        );
                                }

                                // Add contribution from mirrored cell, if any
                                if
                                (
                                    mag(mirrorPlaneTrans[sePointID].first())
                                    > SMALL
                                )
                                {
                                    const tensor& T =
                                        mirrorPlaneTrans[sePointID].second();

                                    coeff += transform(T, coeff);
                                }

                                // Add coeff contribution to globalCoeff
                                pointProcBndFacesCoeffs[cellI][sePointID][i] +=
                                    coeff;

                                neiID++;
                            } // forAll proc boundary faces
                        }


                    } // for the edge start and end points
                } // forAll edges of current face
            } // if internal face - now all faces
        } // forAll faces of current cell

        if (debug > 1)
        {
            Info<< nl << "    total d[" << cellI << "] " << d[cellI] << nl
                << endl;
        }
    } // forAll cells
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
