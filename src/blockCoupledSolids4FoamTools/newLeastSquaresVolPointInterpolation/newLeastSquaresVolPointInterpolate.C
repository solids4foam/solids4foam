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

#include "newLeastSquaresVolPointInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicFvPatch.H"
#include "processorPolyPatch.H"
#include "transform.H"
#ifdef FOAMEXTEND
    #include "cyclicGgiPolyPatch.H"
    #include "cyclicGgiFvPatchFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void newLeastSquaresVolPointInterpolation::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
    const Field<Type>& vfI = vf.primitiveField();
    Field<Type>& pfI = pf.primitiveFieldRef();
#else
    const Field<Type>& vfI = vf.internalField();
    Field<Type>& pfI = pf.internalField();
#endif

    const vectorField& points = mesh().points();

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();

    label nCoeffs = 3;
    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
    const labelListList& ptCyclicBndFaces = pointCyclicBndFaces();
#ifdef FOAMEXTEND
    const labelListList& ptCyclicGgiFaces = pointCyclicGgiFaces();
    const labelListList& ptCyclicGgiBndFaces = pointCyclicGgiBndFaces();
#endif
    const labelListList& ptProcFaces = pointProcFaces();

    Map<Field<Type> > gPtNgbProcBndFaceFieldData;
    globalPointNgbProcBndFaceFieldData(vf, gPtNgbProcBndFaceFieldData);

    Map<Field<Type> > gPtNgbProcCellFieldData;
    globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

#ifdef FOAMEXTEND
    // Cyclic GGI neighbour fields
    FieldField<Field, Type> cyclicGgiNgbFields(mesh().boundary().size());
    forAll(vf.boundaryField(), patchI)
    {
        if (isA<cyclicGgiFvPatchField<Type> >(vf.boundaryField()[patchI]))
        {
            cyclicGgiNgbFields.set
            (
                patchI,
                vf.boundaryField()[patchI].patchNeighbourField()
            );
        }
    }
#endif

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();
    const scalarField& L = refL();

    FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);
    FieldField<Field, Type> procBndFaceVf = procBndFacesFieldData(vf);

    forAll(pfI, pointI)
    {
        const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

        const scalarField& W = w[pointI];

        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
        const labelList& interpCyclicBndFaces = ptCyclicBndFaces[pointI];
#ifdef FOAMEXTEND
        const labelList& interpCyclicGgiFaces = ptCyclicGgiFaces[pointI];
        const labelList& interpCyclicGgiBndFaces = ptCyclicGgiBndFaces[pointI];
#endif
        const labelList& interpProcFaces = ptProcFaces[pointI];

        Field<Type> glInterpNgbProcBndFaceData(0);
        if (gPtNgbProcBndFaceFieldData.found(pointI))
        {
            glInterpNgbProcBndFaceData = gPtNgbProcBndFaceFieldData[pointI];
        }

        Field<Type> glInterpNgbProcCellData(0);
        if (gPtNgbProcCellFieldData.found(pointI))
        {
            glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
        }

        Field<Type> interpNgbProcCellData(0);
        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellData.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellData[cI] =
                    procCellVfI
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        Field<Type> interpNgbProcBndFaceData(0);
        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceData.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceData[fI] =
                    procBndFaceVf
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
          + interpCyclicBndFaces.size()
#ifdef FOAMEXTEND
          + interpCyclicGgiFaces.size()
          + interpCyclicGgiBndFaces.size()
#endif
          + interpProcFaces.size()
          + glInterpNgbProcBndFaceData.size()
          + glInterpNgbProcCellData.size()
          + interpNgbProcCellData.size()
          + interpNgbProcBndFaceData.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        Type avg = pTraits<Type>::zero;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID] = vfI[interpCells[i]];
            avg += sqr(W[pointID])*vfI[interpCells[i]];
            pointID++;
        }

        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            const cyclicFvPatch& cycPatch =
                refCast<const cyclicFvPatch>(mesh().boundary()[patchID]);

            label sizeby2 = faceCells.size()/2;

#ifdef OPENFOAM_ORG
            if (!isA<processorPolyPatch>(cycPatch) || pTraits<Type>::rank == 0)
#else
            if (!(cycPatch.parallel() || pTraits<Type>::rank == 0))
#endif
            {
                if (localFaceID < sizeby2)
                {
                    source[pointID] =
#ifdef OPENFOAM_ORG
                        cycPatch.transform().transformPosition
                        (
                            vfI[faceCells[localFaceID + sizeby2]]
                        );
#else
                        transform
                        (
                            cycPatch.forwardT()[0],
                            vfI[faceCells[localFaceID + sizeby2]]
                        );
#endif
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
#ifdef OPENFOAM_ORG
                        cycPatch.transform().invTransformPosition
                        (
                            vfI[faceCells[localFaceID + sizeby2]]
                        );
#else
                        transform
                        (
                            cycPatch.reverseT()[0],
                            vfI[faceCells[localFaceID - sizeby2]]
                        );
#endif
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
            else
            {
                if (localFaceID < sizeby2)
                {
                    source[pointID] = vfI[faceCells[localFaceID + sizeby2]];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] = vfI[faceCells[localFaceID - sizeby2]];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
        }

        for (label i=0; i<interpCyclicBndFaces.size(); i++)
        {
            label cycFaceID = interpCyclicFaces[0];
            label cycPatchID = mesh().boundaryMesh().whichPatch(cycFaceID);
            label cycStart = mesh().boundaryMesh()[cycPatchID].start();
            label cycLocalFaceID = cycFaceID - cycStart;

            const cyclicPolyPatch& cycPolyPatch =
                refCast<const cyclicPolyPatch>
                (
                    mesh().boundaryMesh()[cycPatchID]
                );

            label faceID = interpCyclicBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);
            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            label sizeby2 = cycPolyPatch.size()/2;

#ifdef OPENFOAM_ORG
            if (!isA<processorPolyPatch>(cycPolyPatch))
#else
            if (!(cycPolyPatch.parallel()))
#endif
            {
                if (cycLocalFaceID < sizeby2)
                {
                    source[pointID] =
#ifdef OPENFOAM_ORG
                        cycPolyPatch.transform().transformPosition
                        (
                            vf.boundaryField()[patchID][localFaceID]
                        );
#else
                        transform
                        (
                            cycPolyPatch.forwardT()[0],
                            vf.boundaryField()[patchID][localFaceID]
                        );
#endif
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
#ifdef OPENFOAM_ORG
                        cycPolyPatch.transform().invTransformPosition
                        (
                            vf.boundaryField()[patchID][localFaceID]
                        );
#else
                        transform
                        (
                            cycPolyPatch.reverseT()[0],
                            vf.boundaryField()[patchID][localFaceID]
                        );
#endif
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
            else
            {
                if (cycLocalFaceID < sizeby2)
                {
                    source[pointID] =
                        vf.boundaryField()[patchID][localFaceID];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
                        vf.boundaryField()[patchID][localFaceID];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
        }

#ifdef FOAMEXTEND
        for (label i=0; i<interpCyclicGgiFaces.size(); i++)
        {
            label faceID = interpCyclicGgiFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = cyclicGgiNgbFields[patchID][localFaceID];
            avg += sqr(W[pointID])*cyclicGgiNgbFields[patchID][localFaceID];
            pointID++;
        }

        if (interpCyclicGgiBndFaces.size())
        {
            label cycGgiFaceID = interpCyclicGgiFaces[0];
            label cycGgiPatchID =
                mesh().boundaryMesh().whichPatch(cycGgiFaceID);
            label cycGgiStart = mesh().boundaryMesh()[cycGgiPatchID].start();
            label cycGgiLocalFaceID = cycGgiFaceID - cycGgiStart;

            const cyclicGgiPolyPatch& cycGgiPatch =
                refCast<const cyclicGgiPolyPatch>
                (
                    mesh().boundaryMesh()[cycGgiPatchID]
                );

            const labelListList& addr
                (
                    cycGgiPatch.master() ?
                    cycGgiPatch.patchToPatch().masterAddr()
                  : cycGgiPatch.patchToPatch().slaveAddr()
                );

            const labelList& zoneAddr =
                cycGgiPatch.zoneAddressing();

            label zoneLocalFaceID = zoneAddr[cycGgiLocalFaceID];

            label ngbZoneLocalFaceID = addr[zoneLocalFaceID][0];

            label ngbLocalFaceID =
                findIndex
                (
                    cycGgiPatch.shadow().zoneAddressing(),
                    ngbZoneLocalFaceID
                );
            if (ngbLocalFaceID == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makePointFaces() const"
                )
                    << "Can not find cyclic ggi patch ngb face id"
                        << abort(FatalError);
            }

            const face& ngbFace =
                cycGgiPatch.shadow().localFaces()[ngbLocalFaceID];

            const vectorField& p =
                cycGgiPatch.shadow().localPoints();
            const labelList& meshPoints =
                cycGgiPatch.shadow().meshPoints();

            scalar S = ngbFace.mag(p);
            scalar L = ::sqrt(S);

            label ngbPointI = -1;
            forAll(ngbFace, pI)
            {
                vector transCurPoint = points[pointI];

                if (!cycGgiPatch.parallel())
                {
                    transCurPoint =
                        transform
                        (
                            cycGgiPatch.reverseT()[0],
                            transCurPoint
                        );
                }

                if (cycGgiPatch.separated())
                {
                    transCurPoint += cycGgiPatch.separationOffset();
                }

                scalar dist =
                    mag
                    (
                        p[ngbFace[pI]]
                      - transCurPoint
                    );

//                 scalar dist =
//                     mag
//                     (
//                         p[ngbFace[pI]]
//                       - transform
//                         (
//                             cycGgiPatch.reverseT()[0],
//                             points[pointI]
//                         )
//                     );

                if (dist < 1e-2*L)
                {
                    ngbPointI = meshPoints[ngbFace[pI]];
                    break;
                }
            }
            if (ngbPointI == -1)
            {
                FatalErrorIn
                (
                    "leastSquaresVolPointInterpolation"
                    "::makePointFaces() const"
                )
                    << "Can not find cyclic ggi patch ngb point"
                        << abort(FatalError);
            }

            for (label i=0; i<interpCyclicGgiBndFaces.size(); i++)
            {
                label faceID = interpCyclicGgiBndFaces[i];
                label patchID =
                    mesh().boundaryMesh().whichPatch(faceID);
                label start = mesh().boundaryMesh()[patchID].start();
                label localFaceID = faceID - start;

                if (cycGgiPatch.parallel())
                {
                    source[pointID] =
                        vf.boundaryField()[patchID][localFaceID];
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
                else
                {
                    source[pointID] =
                        transform
                        (
                            cycGgiPatch.forwardT()[0],
                            vf.boundaryField()[patchID][localFaceID]
                        );
                    avg += sqr(W[pointID])*source[pointID];
                    pointID++;
                }
            }
        }
#endif

        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = glInterpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*glInterpNgbProcBndFaceData[i];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcCellData.size(); i++)
        {
            source[pointID] = glInterpNgbProcCellData[i];
            avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcCellData.size(); i++)
        {
            source[pointID] = interpNgbProcCellData[i];
            avg += sqr(W[pointID])*interpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = interpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
            pointID++;
        }

        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const tensor& T = mirrorPlaneTransformation()[pointI].second();

            label oldSize = source.size();

            source.setSize(2*oldSize);

            for (label i=oldSize; i<source.size(); i++)
            {
                source[i] = transform(T, source[i-oldSize]);
            }

            avg += transform(T, avg);
        }

//         avg /= source.size() + SMALL;
        avg /= sum(sqr(W));

        source -= avg;

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (points[pointI] - o[pointI])/L[pointI];

        pfI[pointI] =
            avg
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    // ZT: Correct only patches which are not coupled
    // (for example, processor patches should already be sinchronized)
    // pf.correctBoundaryConditions();
    forAll(pf.boundaryField(), patchI)
    {
        if (!pf.boundaryField()[patchI].coupled())
        {
            pf.boundaryField()[patchI].initEvaluate();
        }
    }
    forAll(pf.boundaryField(), patchI)
    {
        if (!pf.boundaryField()[patchI].coupled())
        {
            pf.boundaryField()[patchI].evaluate();
        }
    }

    // Correct axis point values
    forAll(mesh().boundaryMesh(), patchI)
    {
        if
        (
            mesh().boundaryMesh()[patchI].type()
         == wedgePolyPatch::typeName
        )
        {
            const labelList& meshPoints =
                mesh().boundaryMesh()[patchI].meshPoints();

            const wedgePolyPatch& wedge =
                refCast<const wedgePolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            vector n =
                transform
                (
                    wedge.faceT(),
                    wedge.centreNormal()
                );
            n /= mag(n);

            forAll(meshPoints, pointI)
            {
                label curMeshPoint = meshPoints[pointI];

                if (pointAxisEdges().found(curMeshPoint))
                {
                    // Keep only component along axis
                    pfI[curMeshPoint] =
                        transform
                        (
                            sqr(wedge.axis()),
                            pfI[curMeshPoint]
                        );
                }
                else
                {
                    pfI[curMeshPoint] =
                        transform
                        (
                            I-sqr(n),
                            pfI[curMeshPoint]
                        );
                }
            }
        }
    }
}


template<class Type>
tmp<Field<Type> > newLeastSquaresVolPointInterpolation::interpolate
(
    const polyPatch& patch,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
    const Field<Type>& vfI = vf.primitiveField();
#else
    const Field<Type>& vfI = vf.internalField();
#endif

    tmp<Field<Type> > tppf
    (
        new Field<Type>
        (
            patch.nPoints(),
            pTraits<Type>::zero
        )
    );
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
    Field<Type>& ppf = tppf.ref();
#else
    Field<Type>& ppf = tppf();
#endif

    const labelList& meshPoints = patch.meshPoints();

    const vectorField& points = mesh().points();

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();

    label nCoeffs = 3;
    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
#ifdef FOAMEXTEND
    const labelListList& ptCyclicGgiFaces = pointCyclicGgiFaces();
#endif
    const labelListList& ptProcFaces = pointProcFaces();

    Map<Field<Type> > gPtNgbProcBndFaceFieldData;
    globalPointNgbProcBndFaceFieldData(vf, gPtNgbProcBndFaceFieldData);

    Map<Field<Type> > gPtNgbProcCellFieldData;
    globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();
    const scalarField& L = refL();

    FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);
    FieldField<Field, Type> procBndFaceVf = procBndFacesFieldData(vf);

    forAll(ppf, pI)
    {
        label pointI = meshPoints[pI];

        const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

        const scalarField& W = w[pointI];

        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
#ifdef FOAMEXTEND
        const labelList& interpCyclicGgiFaces = ptCyclicGgiFaces[pointI];
#endif
        const labelList& interpProcFaces = ptProcFaces[pointI];

        Field<Type> glInterpNgbProcBndFaceData(0);
        if (gPtNgbProcBndFaceFieldData.found(pointI))
        {
            glInterpNgbProcBndFaceData = gPtNgbProcBndFaceFieldData[pointI];
        }

        Field<Type> glInterpNgbProcCellData(0);
        if (gPtNgbProcCellFieldData.found(pointI))
        {
            glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
        }

        Field<Type> interpNgbProcCellData(0);
        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellData.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellData[cI] =
                    procCellVfI
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        Field<Type> interpNgbProcBndFaceData(0);
        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceData.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceData[fI] =
                    procBndFaceVf
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
#ifdef FOAMEXTEND
          + interpCyclicGgiFaces.size()
#endif
          + interpProcFaces.size()
          + glInterpNgbProcBndFaceData.size()
          + glInterpNgbProcCellData.size()
          + interpNgbProcCellData.size()
          + interpNgbProcBndFaceData.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        Type avg = pTraits<Type>::zero;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID] = vfI[interpCells[i]];
            avg += sqr(W[pointID])*vfI[interpCells[i]];
            pointID++;
        }

        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            label sizeby2 = faceCells.size()/2;

            if (localFaceID < sizeby2)
            {
                source[pointID] = vfI[faceCells[localFaceID + sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID + sizeby2]];
                pointID++;
            }
            else
            {
                source[pointID] = vfI[faceCells[localFaceID - sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID - sizeby2]];
                pointID++;
            }
        }

        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = glInterpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*glInterpNgbProcBndFaceData[i];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcCellData.size(); i++)
        {
            source[pointID] = glInterpNgbProcCellData[i];
            avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcCellData.size(); i++)
        {
            source[pointID] = interpNgbProcCellData[i];
            avg += sqr(W[pointID])*interpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = interpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
            pointID++;
        }

        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const tensor& T = mirrorPlaneTransformation()[pointI].second();

            label oldSize = source.size();

            source.setSize(2*oldSize);

            for (label i=oldSize; i<source.size(); i++)
            {
                source[i] = transform(T, source[i-oldSize]);
            }

            avg += transform(T, avg);
        }

//         avg /= source.size() + SMALL;
        avg /= sum(sqr(W));

        source -= avg;

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (points[pointI] - o[pointI])/L[pointI];

        ppf[pI] =
            avg
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    return tppf;
}


// template<class Type>
// tmp<Field<Type> > newLeastSquaresVolPointInterpolation::interpolate
// (
//     const polyPatch& patch,
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// ) const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::interpolate("
//             << "const GeometricField<Type, fvPatchField, volMesh>&, "
//             << "GeometricField<Type, pointPatchField, pointMesh>&) : "
//             << "interpolating field from cells to points"
//             << endl;
//     }

//     Info << "patch cell to point interpolation" << endl;

// #if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
//     const Field<Type>& vfI = vf.primitiveField();
// #else
//     const Field<Type>& vfI = vf.internalField();
// #endif

//     tmp<Field<Type> > tppf
//     (
//         new Field<Type>
//         (
//             patch.nPoints(),
//             pTraits<Type>::zero
//         )
//     );
// #if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
//     Field<Type>& ppf = tppf.ref();
// #else
//     Field<Type>& ppf = tppf();
// #endif

//     const labelList& meshPoints = patch.meshPoints();

//     const vectorField& points = mesh().points();

//     label nCoeffs = 3;
//     const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
//     const labelListList& ptCells = mesh().pointCells();
//     const labelListList& ptBndFaces = pointBndFaces();
//     const labelListList& ptProcFaces = pointProcFaces();

//     const FieldField<Field, scalar>& w = weights();
//     const vectorField& o = origins();

//     Map<Field<Type> > ptNgbProcBndFaceFieldData;
//     pointNgbProcBndFaceFieldData(vf, ptNgbProcBndFaceFieldData);

//     Map<Field<Type> > gPtNgbProcCellFieldData;
//     globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

//     const Map<List<labelPair> >& ptProcCells = pointProcCells();

//     FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);

//     forAll(ppf, pI)
//     {
//         label pointI = meshPoints[pI];

//         const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

//         const scalarField& W = w[pointI];

//         const labelList& interpCells = ptCells[pointI];
//         const labelList& interpBndFaces = ptBndFaces[pointI];
//         const labelList& interpProcFaces = ptProcFaces[pointI];

//         Field<Type> interpNgbProcBndFaceData(0);
//         if (ptNgbProcBndFaceFieldData.found(pointI))
//         {
//             interpNgbProcBndFaceData = ptNgbProcBndFaceFieldData[pointI];
//         }

//         Field<Type> glInterpNgbProcCellData(0);
//         if (gPtNgbProcCellFieldData.found(pointI))
//         {
//             glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
//         }

//         Field<Type> interpNgbProcCellData(0);
//         if (ptProcCells.found(pointI))
//         {
//             const List<labelPair>& pc = ptProcCells[pointI];

//             interpNgbProcCellData.setSize(pc.size());

//             forAll(pc, cI)
//             {
//                 interpNgbProcCellData[cI] =
//                     procCellVfI
//                     [
//                         pc[cI].first()
//                     ]
//                     [
//                         pc[cI].second()
//                     ];
//             }
//         }

//         Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
//         Field<Type> source
//         (
//             interpCells.size()
//           + interpBndFaces.size()
//           + interpProcFaces.size()
//           + interpNgbProcBndFaceData.size()
//           + glInterpNgbProcCellData.size()
//           + interpNgbProcCellData.size(),
//             pTraits<Type>::zero
//         );

//         label pointID = 0;

//         Type avg = pTraits<Type>::zero;

//         for (label i=0; i<interpCells.size(); i++)
//         {
//             source[pointID] = vfI[interpCells[i]];
//             avg += sqr(W[pointID])*vfI[interpCells[i]];
//             pointID++;
//         }

//         for (label i=0; i<interpBndFaces.size(); i++)
//         {
//             label faceID = interpBndFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpProcFaces.size(); i++)
//         {
//             label faceID = interpProcFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
//         {
//             source[pointID] = interpNgbProcBndFaceData[i];
//             avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
//             pointID++;
//         }

//         for (label i=0; i<glInterpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = glInterpNgbProcCellData[i];
//             avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = interpNgbProcCellData[i];
//             avg += sqr(W[pointID])*interpNgbProcCellData[i];
//             pointID++;
//         }

//         if (mirrorPlaneTransformation().found(pointI))
//         {
//             const tensor& T = mirrorPlaneTransformation()[pointI].second();

//             label oldSize = source.size();

//             source.setSize(2*oldSize);

//             for (label i=oldSize; i<source.size(); i++)
//             {
//                 source[i] = transform(T, source[i-oldSize]);
//             }

//             avg += transform(T, avg);
//         }

// //         avg /= source.size() + SMALL;
//         avg /= sum(sqr(W));

//         source -= avg;

//         for (label i=0; i<nCoeffs; i++)
//         {
//             for (label j=0; j<source.size(); j++)
//             {
//                 coeffs[i] += curInvMatrix[i][j]*source[j];
//             }
//         }

//         vector dr = points[pointI] - o[pointI];

//         ppf[pI] =
//             avg
//           + coeffs[0]*dr.x()
//           + coeffs[1]*dr.y()
//           + coeffs[2]*dr.z();
//     }

//     return tppf;
// }


template<class Type>
Type newLeastSquaresVolPointInterpolation::interpolate
(
    const label pointIndex,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::interpolate("
            << "const label, "
            << "const GeometricField<Type, fvPatchField, volMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
    const Field<Type>& vfI = vf.primitiveField();
#else
    const Field<Type>& vfI = vf.internalField();
#endif

    Type pf = pTraits<Type>::zero;

    const vectorField& points = mesh().points();

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();

    label nCoeffs = 3;
    const labelListList& ptCells = mesh().pointCells();
    const labelListList& ptBndFaces = pointBndFaces();
    const labelListList& ptCyclicFaces = pointCyclicFaces();
#ifdef FOAMEXTEND
    const labelListList& ptCyclicGgiFaces = pointCyclicGgiFaces();
#endif
    const labelListList& ptProcFaces = pointProcFaces();

    Map<Field<Type> > gPtNgbProcBndFaceFieldData;
    globalPointNgbProcBndFaceFieldData(vf, gPtNgbProcBndFaceFieldData);

    Map<Field<Type> > gPtNgbProcCellFieldData;
    globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

    const Map<List<labelPair> >& ptProcCells = pointProcCells();

    const List<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();
//     const Map<List<labelPair> >& ptProcBndFaces = pointProcBndFaces();

    const FieldField<Field, scalar>& w = weights();
    const vectorField& o = origins();
    const scalarField& L = refL();

    FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);
    FieldField<Field, Type> procBndFaceVf = procBndFacesFieldData(vf);

    // Calc point field value
    {
        label pointI = pointIndex;

        const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

        const scalarField& W = w[pointI];

        const labelList& interpCells = ptCells[pointI];
        const labelList& interpBndFaces = ptBndFaces[pointI];
        const labelList& interpCyclicFaces = ptCyclicFaces[pointI];
#ifdef FOAMEXTEND
        const labelList& interpCyclicGgiFaces = ptCyclicGgiFaces[pointI];
#endif
        const labelList& interpProcFaces = ptProcFaces[pointI];

        Field<Type> glInterpNgbProcBndFaceData(0);
        if (gPtNgbProcBndFaceFieldData.found(pointI))
        {
            glInterpNgbProcBndFaceData = gPtNgbProcBndFaceFieldData[pointI];
        }

        Field<Type> glInterpNgbProcCellData(0);
        if (gPtNgbProcCellFieldData.found(pointI))
        {
            glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
        }

        Field<Type> interpNgbProcCellData(0);
        if (ptProcCells.found(pointI))
        {
            const List<labelPair>& pc = ptProcCells[pointI];

            interpNgbProcCellData.setSize(pc.size());

            forAll(pc, cI)
            {
                interpNgbProcCellData[cI] =
                    procCellVfI
                    [
                        pc[cI].first()
                    ]
                    [
                        pc[cI].second()
                    ];
            }
        }

        Field<Type> interpNgbProcBndFaceData(0);
        if (ptProcBndFaces[pointI].size())
//         if (ptProcBndFaces.found(pointI))
        {
            const List<labelPair>& pf = ptProcBndFaces[pointI];

            interpNgbProcBndFaceData.setSize(pf.size());

            forAll(pf, fI)
            {
                interpNgbProcBndFaceData[fI] =
                    procBndFaceVf
                    [
                        pf[fI].first()
                    ]
                    [
                        pf[fI].second()
                    ];
            }
        }

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size()
          + interpBndFaces.size()
          + interpCyclicFaces.size()
#ifdef FOAMEXTEND
          + interpCyclicGgiFaces.size()
#endif
          + interpProcFaces.size()
          + glInterpNgbProcBndFaceData.size()
          + glInterpNgbProcCellData.size()
          + interpNgbProcCellData.size()
          + interpNgbProcBndFaceData.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        Type avg = pTraits<Type>::zero;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID] = vfI[interpCells[i]];
            avg += sqr(W[pointID])*vfI[interpCells[i]];
            pointID++;
        }

        for (label i=0; i<interpBndFaces.size(); i++)
        {
            label faceID = interpBndFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<interpCyclicFaces.size(); i++)
        {
            label faceID = interpCyclicFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            const unallocLabelList& faceCells =
                mesh().boundary()[patchID].faceCells();

            label sizeby2 = faceCells.size()/2;

            if (localFaceID < sizeby2)
            {
                source[pointID] = vfI[faceCells[localFaceID + sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID + sizeby2]];
                pointID++;
            }
            else
            {
                source[pointID] = vfI[faceCells[localFaceID - sizeby2]];
                avg += sqr(W[pointID])*vfI[faceCells[localFaceID - sizeby2]];
                pointID++;
            }
        }

        for (label i=0; i<interpProcFaces.size(); i++)
        {
            label faceID = interpProcFaces[i];
            label patchID =
                mesh().boundaryMesh().whichPatch(faceID);

            label start = mesh().boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID] = vf.boundaryField()[patchID][localFaceID];
            avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = glInterpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*glInterpNgbProcBndFaceData[i];
            pointID++;
        }

        for (label i=0; i<glInterpNgbProcCellData.size(); i++)
        {
            source[pointID] = glInterpNgbProcCellData[i];
            avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcCellData.size(); i++)
        {
            source[pointID] = interpNgbProcCellData[i];
            avg += sqr(W[pointID])*interpNgbProcCellData[i];
            pointID++;
        }

        for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
        {
            source[pointID] = interpNgbProcBndFaceData[i];
            avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
            pointID++;
        }

        if (mag(mirrorPlaneTransformation()[pointI].first()) > SMALL)
//         if (mirrorPlaneTransformation().found(pointI))
        {
            const tensor& T = mirrorPlaneTransformation()[pointI].second();

            label oldSize = source.size();

            source.setSize(2*oldSize);

            for (label i=oldSize; i<source.size(); i++)
            {
                source[i] = transform(T, source[i-oldSize]);
            }

            avg += transform(T, avg);
        }

//         avg /= source.size() + SMALL;
        avg /= sum(sqr(W));

        source -= avg;

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector dr = (points[pointI] - o[pointI])/L[pointI];

        pf =
            avg
          + coeffs[0]*dr.x()
          + coeffs[1]*dr.y()
          + coeffs[2]*dr.z();
    }

    return pf;
}


// template<class Type>
// Type newLeastSquaresVolPointInterpolation::interpolate
// (
//     const label pointIndex,
//     const GeometricField<Type, fvPatchField, volMesh>& vf
// ) const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::interpolate("
//             << "const label, "
//             << "const GeometricField<Type, fvPatchField, volMesh>&) : "
//             << "interpolating field from cells to points"
//             << endl;
//     }

// #if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
//     const Field<Type>& vfI = vf.primitiveField();
// #else
//     const Field<Type>& vfI = vf.internalField();
// #endif
//     Type pf = pTraits<Type>::zero;


//     const vectorField& points = mesh().points();

//     label nCoeffs = 3;
//     const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
//     const labelListList& ptCells = mesh().pointCells();
//     const labelListList& ptBndFaces = pointBndFaces();
//     const labelListList& ptProcFaces = pointProcFaces();

//     const FieldField<Field, scalar>& w = weights();
//     const vectorField& o = origins();

//     Map<Field<Type> > ptNgbProcBndFaceFieldData;
//     pointNgbProcBndFaceFieldData(vf, ptNgbProcBndFaceFieldData);

//     Map<Field<Type> > gPtNgbProcCellFieldData;
//     globalPointNgbProcCellFieldData(vf, gPtNgbProcCellFieldData);

//     const Map<List<labelPair> >& ptProcCells = pointProcCells();

//     FieldField<Field, Type> procCellVfI = procCellsFieldData(vfI);

//     // Calc point field value
//     {
//         label pointI = pointIndex;

//         const scalarRectangularMatrix& curInvMatrix = invMatrices[pointI];

//         const scalarField& W = w[pointI];

//         const labelList& interpCells = ptCells[pointI];
//         const labelList& interpBndFaces = ptBndFaces[pointI];
//         const labelList& interpProcFaces = ptProcFaces[pointI];

//         Field<Type> interpNgbProcBndFaceData(0);
//         if (ptNgbProcBndFaceFieldData.found(pointI))
//         {
//             interpNgbProcBndFaceData = ptNgbProcBndFaceFieldData[pointI];
//         }

//         Field<Type> glInterpNgbProcCellData(0);
//         if (gPtNgbProcCellFieldData.found(pointI))
//         {
//             glInterpNgbProcCellData = gPtNgbProcCellFieldData[pointI];
//         }

//         Field<Type> interpNgbProcCellData(0);
//         if (ptProcCells.found(pointI))
//         {
//             const List<labelPair>& pc = ptProcCells[pointI];

//             interpNgbProcCellData.setSize(pc.size());

//             forAll(pc, cI)
//             {
//                 interpNgbProcCellData[cI] =
//                     procCellVfI
//                     [
//                         pc[cI].first()
//                     ]
//                     [
//                         pc[cI].second()
//                     ];
//             }
//         }

//         Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
//         Field<Type> source
//         (
//             interpCells.size()
//           + interpBndFaces.size()
//           + interpProcFaces.size()
//           + interpNgbProcBndFaceData.size()
//           + glInterpNgbProcCellData.size()
//           + interpNgbProcCellData.size(),
//             pTraits<Type>::zero
//         );

//         label pointID = 0;

//         Type avg = pTraits<Type>::zero;

//         for (label i=0; i<interpCells.size(); i++)
//         {
//             source[pointID] = vfI[interpCells[i]];
//             avg += sqr(W[pointID])*vfI[interpCells[i]];
//             pointID++;
//         }

//         for (label i=0; i<interpBndFaces.size(); i++)
//         {
//             label faceID = interpBndFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpProcFaces.size(); i++)
//         {
//             label faceID = interpProcFaces[i];
//             label patchID =
//                 mesh().boundaryMesh().whichPatch(faceID);

//             label start = mesh().boundaryMesh()[patchID].start();
//             label localFaceID = faceID - start;

//             source[pointID] = vf.boundaryField()[patchID][localFaceID];
//             avg += sqr(W[pointID])*vf.boundaryField()[patchID][localFaceID];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcBndFaceData.size(); i++)
//         {
//             source[pointID] = interpNgbProcBndFaceData[i];
//             avg += sqr(W[pointID])*interpNgbProcBndFaceData[i];
//             pointID++;
//         }

//         for (label i=0; i<glInterpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = glInterpNgbProcCellData[i];
//             avg += sqr(W[pointID])*glInterpNgbProcCellData[i];
//             pointID++;
//         }

//         for (label i=0; i<interpNgbProcCellData.size(); i++)
//         {
//             source[pointID] = interpNgbProcCellData[i];
//             avg += sqr(W[pointID])*interpNgbProcCellData[i];
//             pointID++;
//         }

//         if (mirrorPlaneTransformation().found(pointI))
//         {
//             const tensor& T = mirrorPlaneTransformation()[pointI].second();

//             label oldSize = source.size();

//             source.setSize(2*oldSize);

//             for (label i=oldSize; i<source.size(); i++)
//             {
//                 source[i] = transform(T, source[i-oldSize]);
//             }

//             avg += transform(T, avg);
//         }

// //         avg /= source.size() + SMALL;
//         avg /= sum(sqr(W));

//         source -= avg;

//         for (label i=0; i<nCoeffs; i++)
//         {
//             for (label j=0; j<source.size(); j++)
//             {
//                 coeffs[i] += curInvMatrix[i][j]*source[j];
//             }
//         }

//         vector dr = points[pointI] - o[pointI];

//         pf =
//             avg
//           + coeffs[0]*dr.x()
//           + coeffs[1]*dr.y()
//           + coeffs[2]*dr.z();
//     }

//     return pf;
// }


// template<class Type>
// void newLeastSquaresVolPointInterpolation::pointNgbProcBndFaceFieldData
// (
//     const GeometricField<Type, fvPatchField, volMesh>& vf,
//     Map<Field<Type> >& fieldData
// ) const
// {
//     if (debug)
//     {
//         Info<< "newLeastSquaresVolPointInterpolation::"
//             << "pointNgbProcBndFaceFieldData("
//             << "const GeometricField<Type, fvPatchField, volMesh>&) : "
//             << "extracting bnd face data from ngb processors"
//             << endl;
//     }

//     const labelListList& ptBndFaces = pointBndFaces();

//     if (Pstream::parRun())
//     {
//         forAll(mesh().boundaryMesh(), patchI)
//         {
//             if
//             (
//                 mesh().boundaryMesh()[patchI].type()
//              == processorPolyPatch::typeName
//             )
//             {
//                 const processorPolyPatch& procPatch =
//                     refCast<const processorPolyPatch>
//                     (
//                         mesh().boundaryMesh()[patchI]
//                     );

//                 const labelList& bndPoints = procPatch.boundaryPoints();
//                 const labelList& mshPoints = procPatch.meshPoints();

//                 FieldField<Field, Type> bndPointsFaceData
//                 (
//                     bndPoints.size()
//                 );

//                 forAll(bndPoints, pointI)
//                 {
//                     label curPoint = bndPoints[pointI];
//                     label curMeshPoint = mshPoints[curPoint];

//                     const labelList& curPointBndFaces =
//                         ptBndFaces[curMeshPoint];

//                     Field<Type> curPointFaceData(curPointBndFaces.size());

//                     forAll(curPointFaceData, faceI)
//                     {
//                         label faceID = curPointBndFaces[faceI];
//                         label patchID =
//                             mesh().boundaryMesh().whichPatch(faceID);

//                         label start = mesh().boundaryMesh()[patchID].start();
//                         label localFaceID = faceID - start;

//                         curPointFaceData[faceI] =
//                             vf.boundaryField()[patchID][localFaceID];
//                     }

//                     bndPointsFaceData.set
//                     (
//                         pointI,
//                         new Field<Type>(curPointFaceData)
//                     );
//                 }

//                 // Parallel data exchange
//                 {
//                     OPstream toNeighbProc
//                     (
//                         Pstream::blocking,
//                         procPatch.neighbProcNo()
//                         // size of field
//                     );

//                     toNeighbProc << bndPoints << bndPointsFaceData;
//                 }

//                 FieldField<Field, Type> ngbBndPointsFaceData
//                 (
//                     bndPoints.size()
//                 );

//                 labelList ngbBndPoints(bndPoints.size());

//                 {
//                     IPstream fromNeighbProc
//                     (
//                         Pstream::blocking,
//                         procPatch.neighbProcNo()
//                         // size of field
//                     );

//                     fromNeighbProc >> ngbBndPoints >> ngbBndPointsFaceData;
//                 }

//                 const labelList& glPoints =
//                     mesh().globalData().sharedPointLabels();

//                 forAll(bndPoints, pointI)
//                 {
//                     label curPoint = bndPoints[pointI];
//                     label curMeshPoint = mshPoints[curPoint];

//                     label gpIndex = findIndex(glPoints, curMeshPoint);

//                     // non-global points
//                     if (gpIndex == -1)
//                     {
//                         label curNgbPoint = procPatch.neighbPoints()[curPoint];

//                         label curNgbBndPoint =
//                             findIndex(ngbBndPoints, curNgbPoint);

//                         fieldData.insert
//                         (
//                             curMeshPoint,
//                             ngbBndPointsFaceData[curNgbBndPoint]
//                         );
//                     }
//                 }
//             }
//         }

//         // Global boundary points
//         if (mesh().globalData().nGlobalPoints())
//         {
//             const labelList& spLabels =
//                 mesh().globalData().sharedPointLabels();

//             const labelList& spAddressing =
//                 mesh().globalData().sharedPointAddr();

//             for (label k=0; k<mesh().globalData().nGlobalPoints(); k++)
//             {
//                 List<List<Type> > procBndFaceData(Pstream::nProcs());

//                 label curSpIndex = findIndex(spAddressing, k);

//                 if (curSpIndex != -1)
//                 {
//                     label curMeshPoint = spLabels[curSpIndex];

//                     const labelList& curBndFaces = ptBndFaces[curMeshPoint];

//                     procBndFaceData[Pstream::myProcNo()] =
//                         List<Type>(curBndFaces.size());

//                     forAll (curBndFaces, faceI)
//                     {
//                         label faceID = curBndFaces[faceI];
//                         label patchID =
//                             mesh().boundaryMesh().whichPatch(faceID);

//                         label start = mesh().boundaryMesh()[patchID].start();
//                         label localFaceID = faceID - start;

//                         procBndFaceData[Pstream::myProcNo()][faceI] =
//                             vf.boundaryField()[patchID][localFaceID];
//                     }
//                 }
//                 else
//                 {
//                     procBndFaceData[Pstream::myProcNo()] = List<Type>(0);
//                 }

//                 Pstream::gatherList(procBndFaceData);
//                 Pstream::scatterList(procBndFaceData);

//                 if (curSpIndex != -1)
//                 {
//                     label curMeshPoint = spLabels[curSpIndex];

//                     label nAllFaces = 0;
//                     forAll(procBndFaceData, procI)
//                     {
//                         if (procI != Pstream::myProcNo())
//                         {
//                             nAllFaces += procBndFaceData[procI].size();
//                         }
//                     }

//                     Field<Type> allFaceData(nAllFaces, pTraits<Type>::zero);

//                     label counter = 0;
//                     forAll(procBndFaceData, procI)
//                     {
//                         if (procI != Pstream::myProcNo())
//                         {
//                             forAll(procBndFaceData[procI], faceI)
//                             {
//                                 allFaceData[counter++] =
//                                     procBndFaceData[procI][faceI];
//                             }
//                         }
//                     }

//                     fieldData.insert
//                     (
//                         curMeshPoint,
//                         allFaceData
//                     );
//                 }
//             }
//         }
//     }
// }

template<class Type>
void newLeastSquaresVolPointInterpolation::globalPointNgbProcBndFaceFieldData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Map<Field<Type> >& fieldData
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "globalPointNgbProcBndFaceFieldData("
            << "const GeometricField<Type, fvPatchField, volMesh>&) : "
            << "extracting bnd face data from ngb processors"
            << endl;
    }

//     const labelListList& ptCells = mesh().pointCells();

//     const Field<Type>& vfI = vf.internalField();

    if (Pstream::parRun())
    {
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
        //notImplemented("globalPointNgbProcBndFaceFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "globalPointNgbProcBndFaceFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
#else
        const labelListList& ptBndFaces = pointBndFaces();

        // Global points
        if (mesh().globalData().nGlobalPoints())
        {
            const labelList& spLabels =
                mesh().globalData().sharedPointLabels();

            const labelList& spAddressing =
                mesh().globalData().sharedPointAddr();

            for (label k=0; k<mesh().globalData().nGlobalPoints(); k++)
            {
                List<List<Type> > procBndFaceData(Pstream::nProcs());

                label curSpIndex = findIndex(spAddressing, k);

                if (curSpIndex != -1)
                {
                    label curMeshPoint = spLabels[curSpIndex];

                    const labelList& curBndFaces = ptBndFaces[curMeshPoint];

                    procBndFaceData[Pstream::myProcNo()] =
                        List<Type>(curBndFaces.size());

                    forAll (curBndFaces, faceI)
                    {
                        label faceID = curBndFaces[faceI];
                        label patchID =
                            mesh().boundaryMesh().whichPatch(faceID);
                        label start = mesh().boundaryMesh()[patchID].start();
                        label localFaceID = faceID - start;

                        procBndFaceData[Pstream::myProcNo()][faceI] =
                            vf.boundaryField()[patchID][localFaceID];
                    }
                }
                else
                {
                    procBndFaceData[Pstream::myProcNo()] = List<Type>(0);
                }

                Pstream::gatherList(procBndFaceData);
                Pstream::scatterList(procBndFaceData);

                if (curSpIndex != -1)
                {
                    label curMeshPoint = spLabels[curSpIndex];

                    label nAllBndFaces = 0;
                    forAll(procBndFaceData, procI)
                    {
                        if (procI != Pstream::myProcNo())
                        {
                            nAllBndFaces += procBndFaceData[procI].size();
                        }
                    }

                    Field<Type> allBndFaceData
                    (
                        nAllBndFaces,
                        pTraits<Type>::zero
                    );

                    label counter = 0;
                    forAll(procBndFaceData, procI)
                    {
                        if (procI != Pstream::myProcNo())
                        {
                            forAll(procBndFaceData[procI], faceI)
                            {
                                allBndFaceData[counter++] =
                                    procBndFaceData[procI][faceI];
                            }
                        }
                    }

                    fieldData.insert
                    (
                        curMeshPoint,
                        allBndFaceData
                    );
                }
            }
        }
#endif
    }
}

template<class Type>
void newLeastSquaresVolPointInterpolation::globalPointNgbProcCellFieldData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Map<Field<Type> >& fieldData
) const
{
    if (debug)
    {
        Info<< "newLeastSquaresVolPointInterpolation::"
            << "globalPointNgbProcCellFieldData("
            << "const GeometricField<Type, fvPatchField, volMesh>&) : "
            << "extracting cell data from ngb processors"
            << endl;
    }

    if (Pstream::parRun())
    {
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
        // notImplemented("globalPointNgbProcCellFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "globalPointNgbProcCellFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
#else
        const labelListList& ptCells = mesh().pointCells();

#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
        const Field<Type>& vfI = vf.primitiveField();
#else
        const Field<Type>& vfI = vf.internalField();
#endif

        // Global points
        if (mesh().globalData().nGlobalPoints())
        {
            const labelList& spLabels =
                mesh().globalData().sharedPointLabels();

            const labelList& spAddressing =
                mesh().globalData().sharedPointAddr();

            for (label k=0; k<mesh().globalData().nGlobalPoints(); k++)
            {
                List<List<Type> > procCellData(Pstream::nProcs());

                label curSpIndex = findIndex(spAddressing, k);

                if (curSpIndex != -1)
                {
                    label curMeshPoint = spLabels[curSpIndex];

                    const labelList& curCells = ptCells[curMeshPoint];

                    procCellData[Pstream::myProcNo()] =
                        List<Type>(curCells.size());

                    forAll (curCells, cellI)
                    {
                        procCellData[Pstream::myProcNo()][cellI] =
                            vfI[curCells[cellI]];
                    }
                }
                else
                {
                    procCellData[Pstream::myProcNo()] = List<Type>(0);
                }

                Pstream::gatherList(procCellData);
                Pstream::scatterList(procCellData);

                if (curSpIndex != -1)
                {
                    label curMeshPoint = spLabels[curSpIndex];

                    label nAllCells = 0;
                    forAll(procCellData, procI)
                    {
                        if (procI != Pstream::myProcNo())
                        {
                            nAllCells += procCellData[procI].size();
                        }
                    }

                    Field<Type> allCellData(nAllCells, pTraits<Type>::zero);

                    label counter = 0;
                    forAll(procCellData, procI)
                    {
                        if (procI != Pstream::myProcNo())
                        {
                            forAll(procCellData[procI], cellI)
                            {
                                allCellData[counter++] =
                                    procCellData[procI][cellI];
                            }
                        }
                    }

                    fieldData.insert
                    (
                        curMeshPoint,
                        allCellData
                    );
                }
            }
        }
#endif
    }
}


template<class Type>
Foam::tmp<Foam::FieldField<Foam::Field, Type> >
newLeastSquaresVolPointInterpolation::procCellsFieldData
(
    const Field<Type>& psi
) const
{
    tmp<FieldField<Field, Type> > tprocPsi
    (
        new FieldField<Field, Type>(Pstream::nProcs())
    );
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
    FieldField<Field, Type>& procPsi = tprocPsi.ref();
#else
    FieldField<Field, Type>& procPsi = tprocPsi();
#endif

    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                procCellCentres()[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

    if (Pstream::parRun())
    {
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
        //notImplemented("procCellsFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "procCellsFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
#else
        // Send
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not send empty lists
                if (!procCells()[procI].empty())
                {
                    Field<Type> curPsi(psi, procCells()[procI]);

                    // Parallel data exchange
                    OPstream::write
                    (
                        Pstream::blocking,
                        procI,
                        reinterpret_cast<const char*>(curPsi.begin()),
                        curPsi.byteSize()
                    );

//                     OPstream toProc
//                     (
//                         Pstream::blocking,
//                         procI
// //                         curPsi.size()*sizeof(Type)
//                     );

//                     toProc << curPsi;
                }
            }
        }

        // Receive
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not receive empty lists
                if (!procPsi[procI].empty())
                {
                    // Parallel data exchange
                    IPstream::read
                    (
                        Pstream::blocking,
                        procI,
                        reinterpret_cast<char*>(procPsi[procI].begin()),
                        procPsi[procI].byteSize()
                    );

//                     IPstream fromProc
//                     (
//                         Pstream::blocking,
//                         procI
// //                         procPsi[procI].size()*sizeof(Type)
//                     );

//                     fromProc >> procPsi[procI];
                }
            }
        }
#endif
    }

    return tprocPsi;
}


template<class Type>
Foam::tmp<Foam::FieldField<Foam::Field, Type> >
newLeastSquaresVolPointInterpolation::procBndFacesFieldData
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<FieldField<Field, Type> > tprocPsi
    (
        new FieldField<Field, Type>(Pstream::nProcs())
    );
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
    FieldField<Field, Type>& procPsi = tprocPsi.ref();
#else
    FieldField<Field, Type>& procPsi = tprocPsi();
#endif

    forAll (procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>
            (
                procBndFaceCentres()[procI].size(),
                pTraits<Type>::zero
            )
        );
    }

    if (Pstream::parRun())
    {
#if (defined(OPENFOAM) || defined(OPENFOAM_NOT_EXTEND))
        // notImplemented("procBndFacesFieldData");
        if (processorBoundariesExist_)
        {
            FatalErrorIn
            (
                "void newLeastSquaresVolPointInterpolation::"
                "procBndFacesFieldData"
            )   << "Currently only the 'none' decomposition is allowed in "
                << "parallel" << abort(FatalError);
        }
#else
        // Send
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not send empty lists
                if (!procBndFaces()[procI].empty())
                {
//                     Field<Type> curPsi(psi, procBndFaces()[procI]);
                    const labelList& curBndFaces = procBndFaces()[procI];
                    Field<Type> curPsi(curBndFaces.size());
                    forAll(curPsi, faceI)
                    {
                        label faceID = curBndFaces[faceI];
                        label patchID =
                            mesh().boundaryMesh().whichPatch(faceID);
                        label start = mesh().boundaryMesh()[patchID].start();
                        label localFaceID = faceID - start;

                        curPsi[faceI] =
                            vf.boundaryField()[patchID][localFaceID];
                    }

                    // Parallel data exchange
                    OPstream::write
                    (
                        Pstream::blocking,
                        procI,
                        reinterpret_cast<const char*>(curPsi.begin()),
                        curPsi.byteSize()
                    );

//                     OPstream toProc
//                     (
//                         Pstream::blocking,
//                         procI
// //                         curPsi.size()*sizeof(Type)
//                     );

//                     toProc << curPsi;
                }
            }
        }

        // Receive
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Do not receive empty lists
                if (!procPsi[procI].empty())
                {
                    // Parallel data exchange
                    IPstream::read
                    (
                        Pstream::blocking,
                        procI,
                        reinterpret_cast<char*>(procPsi[procI].begin()),
                        procPsi[procI].byteSize()
                    );

//                     IPstream fromProc
//                     (
//                         Pstream::blocking,
//                         procI
// //                         procPsi[procI].size()*sizeof(Type)
//                     );

//                     fromProc >> procPsi[procI];
                }
            }
        }
#endif
    }

    return tprocPsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
