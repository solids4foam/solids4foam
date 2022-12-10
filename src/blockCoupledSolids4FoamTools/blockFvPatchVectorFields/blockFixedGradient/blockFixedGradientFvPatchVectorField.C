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

#include "blockFixedGradientFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blockFixedGradientFvPatchVectorField::blockFixedGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    blockFvPatchVectorField(),
    gradient_(p.size(), vector::zero)
{}


blockFixedGradientFvPatchVectorField::blockFixedGradientFvPatchVectorField
(
    const blockFixedGradientFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    blockFvPatchVectorField(),
    gradient_(ptf.gradient_, mapper)
{}


blockFixedGradientFvPatchVectorField::blockFixedGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    blockFvPatchVectorField(),
    gradient_("gradient", dict, p.size())
{}


blockFixedGradientFvPatchVectorField::blockFixedGradientFvPatchVectorField
(
    const blockFixedGradientFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    blockFvPatchVectorField(),
    gradient_(ptf.gradient_)

{}


blockFixedGradientFvPatchVectorField::blockFixedGradientFvPatchVectorField
(
    const blockFixedGradientFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    blockFvPatchVectorField(),
    gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<vector> > blockFixedGradientFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    FatalErrorIn("gradientBoundaryCoeffs()")
        << "This should not be called" << abort(FatalError);

    // Keep the compiler happy
    return *this;
}


void blockFixedGradientFvPatchVectorField::insertBlockCoeffs
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
) const
{
    // Const reference to polyPatch and the fvMesh
    const polyPatch& ppatch = patch().patch();
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Grab block diagonal
    Field<tensor>& d = blockM.diag().asSquare();

    // Grab linear off-diagonal
    Field<tensor>& l = blockM.lower().asSquare();
    Field<tensor>& u = blockM.upper().asSquare();

    // Mesh fields
    const faceList& faces = mesh.faces();
    const vectorField& points = mesh.points();
    const labelListList& pointCells = mesh.pointCells();

    const bool nonOrthogonalMesh = !mesh.orthogonal();
    //const bool nonOrthogonalMesh = false;
    //Warning
    //    << "    non-orthogonal correction disabled on boundaries" << endl;
    //Warning
    //    << "    least squares dr is zero!" << endl;

    // Least squares vol-to-point interpolation origins, weights and matrices

    const newLeastSquaresVolPointInterpolation& volToPointInterp =
        solidMesh.volToPointInterp();

    const vectorField& origins = volToPointInterp.origins();
    const scalarField& refL = volToPointInterp.refL();
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
        U, gPtNgbProcBndFaceFieldData
    );
    Map<Field<vector> > gPtNgbProcCellFieldData;
    volToPointInterp.globalPointNgbProcCellFieldData
    (
        U, gPtNgbProcCellFieldData
    );
    const Map< List<labelPair> >& pointProcCells =
        volToPointInterp.pointProcCells();
    const List< List<labelPair> >& pointProcBndFaces =
        volToPointInterp.pointProcBndFaces();
    const List<Tuple2<vector, tensor> >& mirrorPlaneTrans =
        volToPointInterp.mirrorPlaneTransformation();

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
    const Map<vector>& pointFixedComp = solidMesh.pointFixedComponent(U);
    const Map<symmTensor>& pointFixedDir = solidMesh.pointFixedDirection(U);

    const fvPatch& fpatch = patch();

    const unallocLabelList& pFaceCells = fpatch.faceCells();
    const vectorField pN = fpatch.nf();
    const scalarField& pMagSf = fpatch.magSf();
    const scalarField& pMu = muf.boundaryField()[fpatch.index()];
    //const scalarField& pLambda = lambdaf.boundaryField()[fpatch.index()];
    const scalarField& pDeltaCoeffs = fpatch.deltaCoeffs();
    const vectorField pDelta = fpatch.delta();

    // Index offset for addressing the diagonal of the boundary faces
    const label start = ppatch.start();

    forAll(fpatch, faceI)
    {
        const label curFaceID = start + faceI;
        //const label varI = offset + faceI;
        const label varI = solidMesh.findOldVariableID(curFaceID);

        const label own = pFaceCells[faceI];

        const vector& faceN = pN[faceI];
        const scalar faceMu = pMu[faceI];
        //const scalar faceLambda = pLambda[faceI];
        const scalar faceMagSf = pMagSf[faceI];
        const scalar faceDeltaCoeff = pDeltaCoeffs[faceI];

        // Normal derivative terms
        // Non-orthogonal correction added with tangential derivative
        // terms
        // Todo: we probably don't need to included faceMu, but it shouldn't
        // make a difference

        const tensor coeff = I*faceMu*faceMagSf*faceDeltaCoeff;

        // Diag contribution for the boundary face
        // We flip the sign on all coefficients so as the diagonal is
        // negative
        // i.e. snGrad is (U_b - U_p)*deltaCoeff
        // so we use -(U_b - U_p)*deltaCoeff

        d[varI] -= coeff;

        // Source contribution
        blockB[varI] -= gradient()[faceI]*faceMagSf;

        // Contribution to the lower from owner cell

        // Find which implicit bond connects cellI to
        // pointCellI

        const label curImplicitBondID =
            solidMesh.findCellFaceImplicitBond(own, curFaceID);

        if (curImplicitBondID == -1)
        {
            FatalErrorIn
                (
                    "pointGaussLsDivSigmaScheme::"
                    "insertCoeffBc(): solidTraction"
                )   << "curImplicitBond not found"
                    << abort(FatalError);
        }

        l[curImplicitBondID] += coeff;


        // Tangential derivative terms

        const face& curFace = faces[curFaceID];
        const edgeList curFaceEdges = curFace.edges();

        // Non-orthogonal correction vector
        const vector* faceKPtr = NULL;
        if (nonOrthogonalMesh)
        {
            //faceKPtr = &((*corrVecIPtr)[curFaceID]);
            const vector& faceDelta = pDelta[faceI];
            //faceKPtr = new vector((I - sqr(faceN)) & faceDelta);
            //faceKPtr = new vector(vector::zero);
            faceKPtr =
                new vector
                (
                    pos(faceN & faceDelta)*(faceN - faceDelta*faceDeltaCoeff)
                );
        }

        forAll(curFaceEdges, edgeI)
        {
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
            //const tensor LeFaceN = Le*faceN;

            // Non-orthogonal correction component
            scalar* kDotLePtr = NULL;
            if (nonOrthogonalMesh)
            {
                kDotLePtr = new scalar((*faceKPtr) & Le);
            }

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
            labelList sePointIDs(2, label(-1));
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

                // ID for neighbour
                label neiID = 0;

                // Vector between current point and the average
                // position of the neighbours
                const vector dr =
                    (points[sePointID] - origins[sePointID])/refL[sePointID];

                // Least square inverse matrix terms
                const scalarRectangularMatrix& curInvMatrix =
                    invMatrices[sePointID];

                // Weights for the current point
                const scalarField& w = weights[sePointID];

                // Calculate sum of weights but we do not include mirror
                // points
                //const scalar sumW = sum(sqr(w));
                const scalar sumW = sum(w);

                // Sum least squares matrix weights
                // vector sumM = vector::zero;
                // const label nNei = w.size();
                // for (label coeffI = 0; coeffI < 3; coeffI++)
                // {
                //     for (label neiI = 0; neiI < nNei; neiI++)
                //     {
                //         sumM[coeffI] += curInvMatrix[coeffI][neiI];
                //     }
                // }

                // Mirror point, if any
                const Tuple2<vector, tensor>& sePointMirr =
                    mirrorPlaneTrans[sePointID];
                bool pointHasMirror = false;
                if (mag(sePointMirr.first()) > SMALL)
                {
                    pointHasMirror = true;
                }

                // Check if the point has any fixed components
                const bool pointHasFixedComp =
                    pointFixedComp.found(sePointID);
                ///const bool pointHasFixedComp = false;
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
                    // correction least squares

                    const vector m
                    (
                        curInvMatrix[0][neiID],
                        curInvMatrix[1][neiID],
                        curInvMatrix[2][neiID]
                    );

                    // Contribution from pure inverse distance average
                    // (wCell) and a correction from least squares due
                    // a nonuniform distribution of points (m & dr)
                    // So least squares correction is zero if the vertex
                    // is located at the average position of the
                    // neighbours
                    tensor coeff = tensor::zero;

                    // Add non-orthogonal correction term
                    if (nonOrthogonalMesh)
                    {
                        coeff +=
                            //coeff -=
                            tensor
                            (
                                0.5*(wCell + (m & dr))*
                                (
                                    I*faceMu*(*kDotLePtr)
                                )
                            );
                    }

                    // Add contribution from mirrored cell, if any
                    if (pointHasMirror)
                    {
                        coeff += transform(sePointMirr.second(), coeff);
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

                    // All pointCell contributions will be to lower as
                    // boundary faces come after internal cells in the
                    // diag

                    // Find which implicit bond connects cellI to
                    // pointCellI
                    const label curImplicitBondID =
                        solidMesh.findCellFaceImplicitBond
                        (
                            sePointCellI,
                            curFaceID
                        );

                    if (curImplicitBondID == -1)
                    {
                        FatalErrorIn
                        (
                            "pointGaussLsDivSigmaScheme::"
                            "insertCoeffBC() solidTraction"
                        )   << "curImplicitBond not found"
                            << abort(FatalError);
                    }

                    if (debug > 1)
                    {
                        Info<< "    l[" << curImplicitBondID << "] p "
                            << coeff << endl;
                    }

                    // Note: the sign is opposite to the internal faces
                    // discretisation in inserCoeffsTang
                    l[curImplicitBondID] -= coeff;

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
                    // correction least squares

                    const vector m
                    (
                        curInvMatrix[0][neiID],
                        curInvMatrix[1][neiID],
                        curInvMatrix[2][neiID]
                    );

                    // Contribution from pure inverse distance average
                    // (wCell) and a correction from least squares due
                    // a nonuniform distribution of points (m & dr)
                    // So least squares correction is zero if the vertex
                    // is located at the average position of the
                    // neighbours
                    tensor coeff = tensor::zero;

                    // Add non-orthogonal correction term
                    if (nonOrthogonalMesh)
                    {
                        coeff +=
                            //coeff -=
                            tensor
                            (
                                0.5*(wCell + (m & dr))*
                                (
                                    I*faceMu*(*kDotLePtr)
                                )
                            );
                    }

                    // Add contribution from mirrored cell, if any
                    if (pointHasMirror)
                    {
                        coeff += transform(sePointMirr.second(), coeff);
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

                    if (sePointBndFaceI == curFaceID)
                    {
                        if (debug > 1)
                        {
                            Info<< "    d[" << varI
                                << "] p " << coeff << endl;
                        }
                        d[varI] -= coeff;
                    }
                    else
                    {
                        // Find which implicit bond connects the
                        // current face to the neighbour boundary
                        // face
                        const label curImplicitBondID =
                            solidMesh.findFaceFaceImplicitBond
                            (
                                curFaceID,
                                sePointBndFaceI
                            );

                        if (curImplicitBondID == -1)
                        {
                            FatalErrorIn
                            (
                                "pointGaussLsDivSigmaScheme::"
                                "insertCoeffBc() solidTraction"
                            )   << "bndFace: curImplicitBond not found"
                                << abort(FatalError);
                        }

                        // Contributions to upper/lower
                        if (curFaceID > sePointBndFaceI)
                        {
                            if (debug > 1)
                            {
                                Info<< "    l[" << curImplicitBondID
                                    << "] p "
                                    << coeff << endl;
                            }
                            l[curImplicitBondID] -= coeff;
                            //l[curImplicitBondID] += coeff;
                        }
                        else
                        {
                            if (debug > 1)
                            {
                                Info<< "    u[" << curImplicitBondID
                                    << "] p "
                                    << coeff << endl;
                            }
                            u[curImplicitBondID] -= coeff;
                            //u[curImplicitBondID] += coeff;
                        }
                    }

                    neiID++;
                } // for boundary faces


                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                //~        Add coeffs for cyclic boundary faces      ~//
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                const labelList& sePointCyclicFaces =
                    pointCyclicFaces[sePointID];

                if (sePointCyclicFaces.size())
                {
                    FatalErrorIn
                        (
                            "void insertCoeffsBc\n"
                            "(\n"
                            "    const solidPolyMesh& solidMesh,\n"
                            "    const surfaceScalarField& muf,\n"
                            "    const surfaceScalarField& lambdaf,\n"
                            "    GeometricField<vector, fvPatchField, "
                            "volMesh>& U,\n"
                            "    Field<vector>& blockB,\n"
                            "    BlockLduMatrix<vector>& blockM\n"
                            ")\n"
                        )   << "cyclic faces not implemented"
                            << abort(FatalError);
                }


                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                //~        Add coeffs for proc boundary faces        ~//
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                // These are the cells across the proc boundary which
                // are accessible using standard proc-proc data transfer

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
                    tensor coeff = tensor::zero;

                    // Add non-orthogonal correction term
                    if (nonOrthogonalMesh)
                    {
                        coeff +=
                            tensor
                            (
                                0.5*(wCell + (m & dr))*
                                (
                                    I*faceMu*(*kDotLePtr)
                                )
                            );
                    }

                    // Add contribution from mirrored cell, if any
                    if (pointHasMirror)
                    {
                        coeff += transform(sePointMirr.second(), coeff);
                    }

                    // Add coeff contribution to globalCoeff
                    pointProcFacesCoeffs[varI][sePointID][i] -= coeff;

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

                        tensor coeff = tensor::zero;

                        // Add non-orthogonal correction term
                        if (nonOrthogonalMesh)
                        {
                            coeff +=
                                tensor
                                (
                                    0.5*(wCell + (m & dr))*
                                    (
                                        I*faceMu*(*kDotLePtr)
                                    )
                                );
                        }

                        // Add contribution from mirrored cell, if any
                        if (pointHasMirror)
                        {
                            coeff +=
                                transform(sePointMirr.second(), coeff);
                        }

                        // Add coeff contribution to globalCoeff
                        gPtNgbProcBndFaceCoeffs[varI][sePointID][i] -=
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

                        tensor coeff = tensor::zero;

                        // Add non-orthogonal correction term
                        if (nonOrthogonalMesh)
                        {
                            coeff +=
                                tensor
                                (
                                    0.5*(wCell + (m & dr))*
                                    (
                                        I*faceMu*(*kDotLePtr)
                                    )
                                );
                        }

                        // Add contribution from mirrored cell, if any
                        if (pointHasMirror)
                        {
                            coeff +=
                                transform(sePointMirr.second(), coeff);
                        }

                        // Add coeff contribution to globalCoeff
                        gPtNgbProcCellCoeffs[varI][sePointID][i] -=
                            coeff;

                        neiID++;
                    } // forAll proc boundary faces
                }


                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                //~             Add coeffs for proc cells            ~//
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                // These are cells on the neighbour proc that are
                // directly accessible through the shadow face

                // This is gathering the data: we will do this in
                // the global patch
                //Field<Type> interpNgbProcCellData(0);

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

                        tensor coeff = tensor::zero;

                        // Add non-orthogonal correction term
                        if (nonOrthogonalMesh)
                        {
                            coeff +=
                                tensor
                                (
                                    0.5*(wCell + (m & dr))*
                                    (
                                        I*faceMu*(*kDotLePtr)
                                    )
                                );
                        }

                        // Add contribution from mirrored cell, if any
                        if (pointHasMirror)
                        {
                            coeff +=
                                transform(sePointMirr.second(), coeff);
                        }

                        // Add coeff contribution to globalCoeff
                        pointProcCellsCoeffs[varI][sePointID][i] -=
                            coeff;

                        neiID++;
                    } // forAll proc boundary faces
                }


                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                //~         Add coeffs for proc boundary faces       ~//
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

                        tensor coeff = tensor::zero;

                        // Add non-orthogonal correction term
                        if (nonOrthogonalMesh)
                        {
                            coeff +=
                                tensor
                                (
                                    0.5*(wCell + (m & dr))*
                                    (
                                        I*faceMu*(*kDotLePtr)
                                    )
                                );
                        }

                        // Add contribution from mirrored cell, if any
                        if (pointHasMirror)
                        {
                            coeff +=
                                transform(sePointMirr.second(), coeff);
                        }

                        // Add coeff contribution to globalCoeff
                        pointProcBndFacesCoeffs[varI][sePointID][i] -=
                            coeff;

                        neiID++;
                    } // forAll proc boundary faces
                }
            } // for the edge start and end points
        } // forAll edges of current face
    } // forAll faces of the patch
}


void blockFixedGradientFvPatchVectorField::write(Ostream& os) const
{
    gradient_.writeEntry("gradient", os);
    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    blockFixedGradientFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
