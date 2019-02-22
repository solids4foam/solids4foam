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

#include "blockFixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "blockTangentialCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blockFixedDisplacementZeroShearFvPatchVectorField::
blockFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    blockFvPatchVectorField(),
    dispSeries_()
{}


blockFixedDisplacementZeroShearFvPatchVectorField::
blockFixedDisplacementZeroShearFvPatchVectorField
(
    const blockFixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    blockFvPatchVectorField(),
    dispSeries_(ptf.dispSeries_)
{}


blockFixedDisplacementZeroShearFvPatchVectorField::
blockFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    blockFvPatchVectorField(),
    dispSeries_()
{
    fvPatchVectorField::operator=(patchInternalField());

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }

    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));
    }
}


blockFixedDisplacementZeroShearFvPatchVectorField::
blockFixedDisplacementZeroShearFvPatchVectorField
(
    const blockFixedDisplacementZeroShearFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    blockFvPatchVectorField(),
    dispSeries_(ptf.dispSeries_)
{}


blockFixedDisplacementZeroShearFvPatchVectorField::
blockFixedDisplacementZeroShearFvPatchVectorField
(
    const blockFixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    blockFvPatchVectorField(),
    dispSeries_(ptf.dispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> >
blockFixedDisplacementZeroShearFvPatchVectorField::snGrad() const
{
    // snGrad without non-orthogonal correction
    // return (*this - patchInternalField())*this->patch().deltaCoeffs();

    // Lookup previous boundary gradient
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    // Calculate correction vector
    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = ((I - sqr(n)) & delta);

    return
    (
        //*this - patchInternalField()
        //*this - (patchInternalField() + (k & gradField.patchInternalField()))
        //*this - (patchInternalField() + (k & gradField))
        (*this - (k & gradField)) - patchInternalField()
    )*this->patch().deltaCoeffs();
}

tmp<Field<vector> > blockFixedDisplacementZeroShearFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    FatalErrorIn("gradientBoundaryCoeffs()")
        << "This should not be called" << abort(FatalError);

    // Keep the compiler happy
    return *this;
}


void blockFixedDisplacementZeroShearFvPatchVectorField::insertBlockCoeffs
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
) const
{
    // Update displacement if time-varying
    vectorField disp = *this;
    if (dispSeries_.size())
    {
        disp = dispSeries_(this->db().time().timeOutputValue());
    }

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

    // Least squares vol-to-point interpolation origins, weights and matrices
    const newLeastSquaresVolPointInterpolation& volToPointInterp =
        solidMesh.volToPointInterp();

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
    const bool enforcePointConstraints = false;
    const Map<vector>& pointFixedComp = solidMesh.pointFixedComponent(U);
    const Map<symmTensor>& pointFixedDir = solidMesh.pointFixedDirection(U);

    // This is effecitively a symmetry plane with a non-zero normal
    // displacement i.e. a moving symmetry plane

    const fvPatch& fpatch = patch();

    const unallocLabelList& pFaceCells = fpatch.faceCells();
    const vectorField pN = fpatch.nf();
    const scalarField& pMagSf = fpatch.magSf();
    const scalarField& pMu = muf.boundaryField()[fpatch.index()];
    const scalarField& pLambda = lambdaf.boundaryField()[fpatch.index()];
    const scalarField& pDeltaCoeffs = fpatch.deltaCoeffs();
    const vectorField pDelta = fpatch.delta();

    // Index offset for addressing the diagonal of the boundary faces
    const label start = ppatch.start();

    // We will calculate the average of the current diagonal to scale the
    // coefficients for the fixedValue boundary conditions
    // If we don't do this then the convergence can be worse and also the
    // results can be strange
    // Note that this scale factor is just approximate and its exact value
    // does not matter
    const tensor averageDiag = gAverage(d);
    scalar diagSign =
        (averageDiag.xx() + averageDiag.yy() + averageDiag.zz());
    diagSign /= mag(diagSign);
    const scalar scaleFac = diagSign*(1.0/sqrt(3.0))*mag(averageDiag);

    if (mag(scaleFac) < SMALL)
    {
        FatalErrorIn
        (
            "void blockFixedDisplacementZeroShearFvPatchVectorField"
            "::insertBlockCoeffs\n"
            "(\n"
            "    const solidPolyMesh& solidMesh,\n"
            "    const surfaceScalarField& muf,\n"
            "    const surfaceScalarField& lambdaf,\n"
            "    const GeometricField<vector, fvPatchField, volMesh>& U,\n"
            "    Field<vector>& blockB,\n"
            "    BlockLduMatrix<vector>& blockM\n"
            ") const"
        )   << "The average diagonal coefficient is zero! The internal faces "
            << "should be discretised before inserting the boundary condition "
            << "equations" << abort(FatalError);
    }

    forAll(fpatch, faceI)
    {
        const label curFaceID = start + faceI;
        const label varI = solidMesh.findOldVariableID(curFaceID);

        const label own = pFaceCells[faceI];

        const vector& faceN = pN[faceI];
        const scalar faceMu = pMu[faceI];
        const scalar faceLambda = pLambda[faceI];
        //const scalar faceTwoMuLambda = 2*faceMu + faceLambda;
        const scalar faceMagSf = pMagSf[faceI];
        const scalar faceDeltaCoeff = pDeltaCoeffs[faceI];

        // Normal derivative terms

        const tensor coeff =
            //sqr(faceN)*faceTwoMuLambda*faceMagSf*faceDeltaCoeff +
            (I - sqr(faceN))*faceMu*faceMagSf*faceDeltaCoeff;

        // Diag contribution for the boundary face
        d[varI] -= coeff;

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
                "insertCoeffBc(): blockFixedDisplacementZeroShear"
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
            const tensor LeFaceN = Le*faceN;

            // Non-orthogonal correction component
            scalar kDotLe = 0.0;
            if (nonOrthogonalMesh)
            {
                kDotLe = (*faceKPtr) & Le;
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
                const vector dr = points[sePointID] - origins[sePointID];

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
                    // Coefficient is composed of the contribution from
                    // pure averaging (inverse distance) and a
                    // correction least squares
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate coefficient
                    tensor coeff = tensor::zero;
                    fv::blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB[varI],
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
                        3
                    );

                    // Remove normal component
                    coeff = (I - sqr(faceN)) & coeff;

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
                            "insertCoeffBC() "
                            "blockFixedDisplacementZeroShear"
                        )   << "curImplicitBond not found"
                            << abort(FatalError);
                    }

                    if (debug > 1)
                    {
                        Info<< "    l[" << curImplicitBondID << "] p "
                            << coeff << endl;
                    }

                    l[curImplicitBondID] -= coeff;

                    neiID++;
                } // forAll pointCells


                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
                //~        Add coeffs for local boundary faces       ~//
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

                const labelList& sePointBndFaces = pointBndFaces[sePointID];

                for (label i = 0; i < sePointBndFaces.size(); i++)
                {
                    // Coefficient is composed of the contribution from
                    // pure averaging (inverse distance) and a
                    // correction least squares
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate coefficient
                    tensor coeff = tensor::zero;
                    fv::blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB[varI],
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
                        3
                    );

                    // Remove normal component
                    coeff = (I - sqr(faceN)) & coeff;

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
                        // current face to the neighbour boundary face
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
                                "insertCoeffBc()"
                                " blockFixedDisplacementZeroShear"
                            )   << "bndFace: curImplicitBond not found"
                                << abort(FatalError);
                        }

                        // Contributions to upper/lower
                        if (curFaceID > sePointBndFaceI)
                        {
                            if (debug > 1)
                            {
                                Info<< "    l[" << curImplicitBondID << "] p "
                                    << coeff << endl;
                            }
                            l[curImplicitBondID] -= coeff;
                        }
                        else
                        {
                            if (debug > 1)
                            {
                                Info<< "    u[" << curImplicitBondID << "] p "
                                    << coeff << endl;
                            }
                            u[curImplicitBondID] -= coeff;
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

                const labelList& sePointProcFaces = pointProcFaces[sePointID];

                for (label i = 0; i < sePointProcFaces.size(); i++)
                {
                    // Coefficient is composed of the contribution from
                    // pure averaging (inverse distance) and a
                    // correction least squares
                    const vector m =
                        vector
                        (
                            curInvMatrix[0][neiID],
                            curInvMatrix[1][neiID],
                            curInvMatrix[2][neiID]
                        );

                    // Calculate coefficient
                    tensor coeff = tensor::zero;
                    fv::blockFvmCalculateCoeff
                    (
                        coeff,
                        blockB[varI],
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
                        3
                    );

                    // Remove normal component
                    coeff = (I - sqr(faceN)) & coeff;

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
                        // Coefficient is composed of the contribution from
                        // pure averaging (inverse distance) and a
                        // correction least squares
                        const vector m =
                            vector
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                        // Calculate coefficient
                        tensor coeff = tensor::zero;
                        fv::blockFvmCalculateCoeff
                        (
                            coeff,
                            blockB[varI],
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
                            3
                        );

                        // Remove normal component
                        coeff = (I - sqr(faceN)) & coeff;

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
                        // Coefficient is composed of the contribution from
                        // pure averaging (inverse distance) and a
                        // correction least squares
                        const vector m =
                            vector
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                        // Calculate coefficient
                        tensor coeff = tensor::zero;
                        fv::blockFvmCalculateCoeff
                        (
                            coeff,
                            blockB[varI],
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
                            3
                        );

                        // Remove normal component
                        coeff = (I - sqr(faceN)) & coeff;

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
                        // Coefficient is composed of the contribution from
                        // pure averaging (inverse distance) and a
                        // correction least squares
                        const vector m =
                            vector
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                        // Calculate coefficient
                        tensor coeff = tensor::zero;
                        fv::blockFvmCalculateCoeff
                        (
                            coeff,
                            blockB[varI],
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
                            3
                        );

                        // Remove normal component
                        coeff = (I - sqr(faceN)) & coeff;

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

                    for (label i = 0; i < sePointProcBndFaces.size(); i++)
                    {
                        // Coefficient is composed of the contribution from
                        // pure averaging (inverse distance) and a
                        // correction least squares
                        const vector m =
                            vector
                            (
                                curInvMatrix[0][neiID],
                                curInvMatrix[1][neiID],
                                curInvMatrix[2][neiID]
                            );

                        // Calculate coefficient
                        tensor coeff = tensor::zero;
                        fv::blockFvmCalculateCoeff
                        (
                            coeff,
                            blockB[varI],
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
                            3
                        );

                        // Remove normal component
                        coeff = (I - sqr(faceN)) & coeff;

                        // Add coeff contribution to globalCoeff
                        pointProcBndFacesCoeffs[varI][sePointID][i] -=
                            coeff;

                        neiID++;
                    } // forAll proc boundary faces
                }


            } // for the edge start and end points
        } // forAll edges of current face

        // Now we enforce the normal displacement (fixedValue) condition

        // Diagonal contribution for the boundary face
        d[varI] = scaleFac*sqr(faceN) + ((I - sqr(faceN)) & d[varI]);

        // Source contribution
        // This is the only difference with a symmetry plane: we force a
        // non-zero normal dislacement
        blockB[varI] += scaleFac*(sqr(faceN) & disp[faceI]);
    }
}


void blockFixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    blockFixedDisplacementZeroShearFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
