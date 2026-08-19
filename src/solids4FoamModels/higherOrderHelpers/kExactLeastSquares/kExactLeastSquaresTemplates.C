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

#include "kExactLeastSquares.H"

namespace Foam
{

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class Type, class GlobalCells>
inline Type kExactLeastSquares::fieldValue
(
    const label globalCellID,
    const GlobalCells& globalCells,
    const Map<FixedList<label, 2>>& remoteLocation,
    const UList<Type>& localField,
    const List<Field<Type>>& remoteField
) const
{
    if (globalCells.isLocal(globalCellID))
    {
        return localField[globalCells.toLocal(globalCellID)];
    }

    // This is used in hot loop so we will do this check only for debugging
    if (debug)
    {
        if (!remoteLocation.found(globalCellID))
        {
            FatalErrorInFunction
                << "Remote location mapping not found for global cell ID "
                << globalCellID
                << abort(FatalError);
        }
    }

    const FixedList<label, 2>& loc = remoteLocation[globalCellID];
    return remoteField[loc[0]][loc[1]];
}


template<class Type>
void kExactLeastSquares::exchangeProcessorPatchField
(
    const processorFvPatch& procPatch,
    const Field<Type>& sendField,
    Field<Type>& receiveField
) const
{
    if (Pstream::myProcNo() < procPatch.neighbProcNo())
    {
        procPatch.send(Pstream::commsTypes::blocking, sendField);
        procPatch.receive(Pstream::commsTypes::blocking, receiveField);
    }
    else
    {
        procPatch.receive(Pstream::commsTypes::blocking, receiveField);
        procPatch.send(Pstream::commsTypes::blocking, sendField);
    }
}


template<class Type>
void kExactLeastSquares::fGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    CompactListList<typename outerProduct<vector, Type>::type>& result
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    // Preliminaries
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const Map<FixedList<label, 2>>& remoteLoc =  stencil().remoteCellLocation();
    const CompactListList<label>& stencils = faceGradStencil();
    const CompactListList<label>& cellStencils = stencil().cellsStencil();
    const List<CompactListList<vector>>& coeffs = faceGradCoeffs();
    const CompactListList<point>& faceQuadPoints =
        quadrature().faceQuadPoints();

    if (result.size() != faceQuadPoints.size())
    {
        FatalErrorInFunction
            << "Result has " << result.size() << " face rows, but "
            << faceQuadPoints.size() << " are required"
            << abort(FatalError);
    }

    // Get values from remote processors
    const UList<Type>& vfI = vf.internalField();
    const List<Field<Type>> remoteField = stencil().remoteFieldPerProc(vfI);

    // Initialise return field to zero
    forAll(result, faceI)
    {
        if (result[faceI].size() != faceQuadPoints[faceI].size())
        {
            FatalErrorInFunction
                << "Result face " << faceI << " has "
                << result[faceI].size() << " quadrature entries, but "
                << faceQuadPoints[faceI].size() << " are required"
                << abort(FatalError);
        }

        forAll(result[faceI], qpI)
        {
            result[faceI][qpI] = pTraits<GradType>::zero;
        }
    }

    // Loop for internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        const UList<label>& stencil = stencils[faceI];
        const CompactListList<vector>& faceCoeffs = coeffs[faceI];

        forAll(faceQuadPoints[faceI], qpI)
        {
            const UList<vector>& qpCoeffs = faceCoeffs[qpI];

            forAll(stencil, stencilI)
            {
                result[faceI][qpI] +=
                    qpCoeffs[stencilI]
                  * fieldValue
                    (
                        stencil[stencilI],
                        globalCells,
                        remoteLoc,
                        vfI,
                        remoteField
                    );
            }
        }
    }

    // Average owner-side reconstructions across processor patches. Exchange
    // quadrature-point coordinates as the point order can be reversed on the
    // neighbouring processor.
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (!isA<processorFvPatch>(patch))
        {
            continue;
        }

        const processorFvPatch& procPatch =
            refCast<const processorFvPatch>(patch);
        const labelUList& faceCells = patch.faceCells();
        const label patchStart = patch.start();

        labelField sendSizes(patch.size(), 0);
        label nSendValues = 0;
        forAll(patch, patchFaceI)
        {
            const label faceI = patchStart + patchFaceI;
            sendSizes[patchFaceI] = faceQuadPoints[faceI].size();
            nSendValues += sendSizes[patchFaceI];
        }

        pointField sendPoints(nSendValues);
        Field<GradType> sendGrad
        (
            nSendValues,
            pTraits<GradType>::zero
        );

        label sendI = 0;
        forAll(patch, patchFaceI)
        {
            const label faceI = patchStart + patchFaceI;
            const label ownCellI = faceCells[patchFaceI];
            const labelUList& cellStencil = cellStencils[ownCellI];
            vectorField ownCoeffs(cellStencil.size());

            forAll(faceQuadPoints[faceI], qpI)
            {
                const point& qp = faceQuadPoints[faceI][qpI];
                sendPoints[sendI] = qp;
                cellGradCoeffsAtPoint(ownCellI, qp, ownCoeffs);

                forAll(cellStencil, stencilI)
                {
                    sendGrad[sendI] +=
                        ownCoeffs[stencilI]
                      *
                        (
                            fieldValue
                            (
                                cellStencil[stencilI],
                                globalCells,
                                remoteLoc,
                                vfI,
                                remoteField
                            )
                          - vfI[ownCellI]
                        );
                }

                ++sendI;
            }
        }

        labelField receiveSizes(patch.size(), 0);
        exchangeProcessorPatchField(procPatch, sendSizes, receiveSizes);

        label nReceiveValues = 0;
        forAll(receiveSizes, patchFaceI)
        {
            nReceiveValues += receiveSizes[patchFaceI];
        }

        pointField receivePoints(nReceiveValues);
        Field<GradType> receiveGrad(nReceiveValues);
        exchangeProcessorPatchField(procPatch, sendPoints, receivePoints);
        exchangeProcessorPatchField(procPatch, sendGrad, receiveGrad);

        label sendStart = 0;
        label receiveStart = 0;
        forAll(patch, patchFaceI)
        {
            const label faceI = patchStart + patchFaceI;

            if (sendSizes[patchFaceI] != receiveSizes[patchFaceI])
            {
                FatalErrorInFunction
                    << "Processor patch " << patch.name() << " face "
                    << patchFaceI << " has " << sendSizes[patchFaceI]
                    << " local quadrature points but "
                    << receiveSizes[patchFaceI]
                    << " neighbour quadrature points"
                    << abort(FatalError);
            }

            boolList pointUsed(receiveSizes[patchFaceI], false);

            forAll(faceQuadPoints[faceI], qpI)
            {
                const point& qp = faceQuadPoints[faceI][qpI];
                scalar minDistance = GREAT;
                label matchI = -1;

                for
                (
                    label neighbourQpI = 0;
                    neighbourQpI < receiveSizes[patchFaceI];
                    ++neighbourQpI
                )
                {
                    if (!pointUsed[neighbourQpI])
                    {
                        const scalar distance = mag
                        (
                            qp
                          - receivePoints[receiveStart + neighbourQpI]
                        );

                        if (distance < minDistance)
                        {
                            minDistance = distance;
                            matchI = neighbourQpI;
                        }
                    }
                }

                const scalar matchTolerance =
                    1000.0*SMALL*max(scalar(1), mag(qp));

                if (matchI < 0 || minDistance > matchTolerance)
                {
                    FatalErrorInFunction
                        << "Cannot match quadrature point " << qpI
                        << " on processor patch " << patch.name()
                        << " face " << patchFaceI << nl
                        << "Minimum point distance = " << minDistance
                        << ", tolerance = " << matchTolerance
                        << abort(FatalError);
                }

                pointUsed[matchI] = true;
                result[faceI][qpI] = 0.5
                   *(
                        sendGrad[sendStart + qpI]
                      + receiveGrad[receiveStart + matchI]
                    );
            }

            sendStart += sendSizes[patchFaceI];
            receiveStart += receiveSizes[patchFaceI];
        }
    }
}


template<class Type>
tmp<GeometricField<
    typename outerProduct<vector, Type>::type,
    fvPatchField,
    volMesh
>> kExactLeastSquares::grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil().remoteCellLocation();

    const CompactListList<label>& stencils = stencil().cellsStencil();
    const CompactListList<vector>& gradCoeffs = this->cellGradCoeffs();

    const UList<Type>& vfI = vf.internalField();
    const List<Field<Type>> remoteField = stencil().remoteFieldPerProc(vfI);

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "grad(" + vf.name() + ")",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                dimless/dimLength * vf.dimensions(),
                pTraits<GradType>::zero
            ),
            "zeroGradient"
        )
    );

    GeometricField<GradType, fvPatchField, volMesh>& grad = tGrad.ref();

    forAll(stencils, cellI)
    {
        const UList<label>& stencil = stencils[cellI];
        const UList<vector>& coeffs = gradCoeffs[cellI];

        // Neighbour cell-average difference contribution
        forAll(stencil, cI)
        {
            grad[cellI] +=
                coeffs[cI]
              *
                (
                    fieldValue
                    (
                        stencil[cI],
                        globalCells,
                        remoteLoc,
                        vfI,
                        remoteField
                    )
                  - vfI[cellI]
                );
        }
    }

    grad.correctBoundaryConditions();

    return tGrad;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
