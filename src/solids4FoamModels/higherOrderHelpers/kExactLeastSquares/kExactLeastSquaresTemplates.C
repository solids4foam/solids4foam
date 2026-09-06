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
Type kExactLeastSquares::evaluateAtPoint
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const label cellID,
    const point& x
) const
{
    if (cellID >= mesh_.nCells())
    {
        FatalErrorInFunction
            << "Cell ID " << cellID << " is outside the local cell range 0 to "
            << mesh_.nCells() - 1 << abort(FatalError);
    }

    const globalIndex& globalCells = stencil().globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil().remoteCellLocation();
    auto& stencils =
        compactListListCRef(stencil().cellsStencil());

    // Construct all coefficient and moment tables collectively before
    // checking for the negative cell marker used by parallel point evaluation.
    (void)cellGradCoeffs();
    if (polynomialOrder() >= 2)
    {
        (void)cellSecondGradCoeffs();
        (void)secondOrderCellMoments();
    }
    if (polynomialOrder() >= 3)
    {
        (void)cellThirdGradCoeffs();
        (void)thirdOrderCellMoments();
    }

    const UList<Type>& vfI = vf.internalField();
    const List<Field<Type>> remoteField = stencil().remoteFieldPerProc(vfI);

    // A negative cell ID marks a point owned by another processor. All
    // processors still perform the coefficient setup and field exchange above
    // because these operations can require collective communication.
    if (cellID < 0)
    {
        return pTraits<Type>::zero;
    }

    const UList<label>& cellStencil = stencils[cellID];
    scalarField valueCoeffs(cellStencil.size() + 1);
    cellValueCoeffsAtPoint(cellID, x, valueCoeffs);
    Type value = valueCoeffs[cellStencil.size()]*vfI[cellID];

    forAll(cellStencil, cI)
    {
        value +=
            valueCoeffs[cI]
           *fieldValue
            (
                cellStencil[cI],
                globalCells,
                remoteLoc,
                vfI,
                remoteField
            );
    }

    return value;
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
    auto& stencils = compactListListCRef(faceGradStencil());
    auto& cellStencils =
        compactListListCRef(stencil().cellsStencil());
#ifdef FOAMEXTEND
    List<CompactListList<vector>>& coeffs =
        const_cast<List<CompactListList<vector>>&>(faceGradCoeffs());
#else
    const List<CompactListList<vector>>& coeffs = faceGradCoeffs();
#endif
    const List<List<FixedList<label, 2>>>& boundaryDataStencil =
        faceBoundaryDataStencil();
#ifdef FOAMEXTEND
    List<CompactListList<vector>>& boundaryDataCoeffs =
        const_cast<List<CompactListList<vector>>&>
        (
            faceBoundaryDataCoeffs()
        );
#else
    const List<CompactListList<vector>>& boundaryDataCoeffs =
        faceBoundaryDataCoeffs();
#endif
    auto& faceQuadPoints =
        compactListListCRef(quadrature().faceQuadPoints());

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
        auto& faceCoeffs = coeffs[faceI];

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

    // Loop over plain boundary faces (e.g. traction/Neumann patches) that
    // are neither coupled, symmetry, nor fixed-value: calcFaceGradCoeffs()
    // has already populated a one-sided reconstruction for these using the
    // owner cell's own polynomial, addressed exactly like an internal face's
    // owner side, so the same evaluation applies unmodified.
    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        const bool symmetryPatch =
            isA<symmetryPolyPatch>(pp)
#ifndef FOAMEXTEND
         || isA<symmetryPlanePolyPatch>(pp)
#endif
        ;

        if
        (
            pp.coupled()
         || symmetryPatch
         || includePatchInStencils_[patchI]
        )
        {
            continue;
        }

        const label patchStart = pp.start();

        forAll(pp, patchFaceI)
        {
            const label faceI = patchStart + patchFaceI;
            const UList<label>& stencil = stencils[faceI];
            auto& faceCoeffs = coeffs[faceI];

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
    }

    // Average owner-side reconstructions across processor patches. The face
    // quadrature object guarantees the same point order on both sides.
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
#ifdef FOAMEXTEND
        const label patchStart = patch.patch().start();
#else
        const label patchStart = patch.start();
#endif

        label nValues = 0;
        forAll(patch, patchFaceI)
        {
            const label faceI = patchStart + patchFaceI;
            nValues += faceQuadPoints[faceI].size();
        }

        Field<GradType> sendGrad
        (
            nValues,
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

        Field<GradType> receiveGrad(nValues);
        exchangeProcessorPatchField(procPatch, sendGrad, receiveGrad);

        label start = 0;
        forAll(patch, patchFaceI)
        {
            const label faceI = patchStart + patchFaceI;

            forAll(faceQuadPoints[faceI], qpI)
            {
                result[faceI][qpI] = 0.5
                   *(
                        sendGrad[start + qpI]
                      + receiveGrad[start + qpI]
                    );
            }

            start += faceQuadPoints[faceI].size();
        }
    }

    // Evaluate one-sided reconstructions on symmetry faces. Coefficients are
    // ordered as physical neighbours, mirrored owner, and mirrored
    // neighbours. Mirrored values are transformed with the same reflection
    // tensor used to construct the geometrical stencil.
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if
        (
           !isA<symmetryPolyPatch>(pp)
#ifndef FOAMEXTEND
         && !isA<symmetryPlanePolyPatch>(pp)
#endif
        )
        {
            continue;
        }

        const labelUList& faceCells = patch.faceCells();
        const vectorField patchNormals(patch.nf());

        forAll(patch, patchFaceI)
        {
#ifdef FOAMEXTEND
            const label faceI = patch.patch().start() + patchFaceI;
#else
            const label faceI = patch.start() + patchFaceI;
#endif
            const label ownCellI = faceCells[patchFaceI];
            const labelUList& cellStencil = cellStencils[ownCellI];
            const label Nn = cellStencil.size();
            const tensor R =
                I - 2.0*sqr(patchNormals[patchFaceI]);
            const Type& ownerValue = vfI[ownCellI];
            const Type mirroredOwnerValue = transform(R, ownerValue);
            auto& faceCoefficients =
                coeffs[faceI];

            forAll(faceQuadPoints[faceI], qpI)
            {
                const UList<vector>& qpCoefficients =
                    faceCoefficients[qpI];

                if (qpCoefficients.size() != 2*Nn + 1)
                {
                    FatalErrorInFunction
                        << "Symmetry face " << faceI << " has "
                        << qpCoefficients.size() << " coefficients, but "
                        << 2*Nn + 1 << " are required"
                        << abort(FatalError);
                }

                forAll(cellStencil, stencilI)
                {
                    const Type value = fieldValue
                    (
                        cellStencil[stencilI],
                        globalCells,
                        remoteLoc,
                        vfI,
                        remoteField
                    );

                    result[faceI][qpI] +=
                        qpCoefficients[stencilI]
                       *(value - ownerValue);

                    result[faceI][qpI] +=
                        qpCoefficients[Nn + 1 + stencilI]
                       *(transform(R, value) - ownerValue);
                }

                result[faceI][qpI] +=
                    qpCoefficients[Nn]
                   *(mirroredOwnerValue - ownerValue);
            }
        }
    }

    // Read all prescribed quadrature-point values once. Boundary-cell
    // reconstructions can depend on fixed-value points from more than one
    // boundary face or patch.
    List<Field<Type>> prescribedValues(mesh.nFaces());

    forAll(vf.boundaryField(), patchI)
    {
        if (!includePatchInStencils_[patchI])
        {
            continue;
        }

        const fvPatchField<Type>& pf = vf.boundaryField()[patchI];

        if (!pf.fixesValue())
        {
            FatalErrorInFunction
                << "Patch " << patchI << " was included as fixed-value "
                << "reconstruction data, but field " << vf.name()
                << " does not fix its value on that patch"
                << abort(FatalError);
        }

        autoPtr<CompactListList<Type>> patchValuesPtr =
            patchFaceQuadValues(vf, patchI);
#ifdef FOAMEXTEND
        CompactListList<Type>& patchValues = autoPtrRef(patchValuesPtr);
#else
        const CompactListList<Type>& patchValues = patchValuesPtr();
#endif
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if (patchValues.size() != pp.size())
        {
            FatalErrorInFunction
                << "Patch quadrature values have " << patchValues.size()
                << " faces but patch " << pp.name() << " has " << pp.size()
                << abort(FatalError);
        }

        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;

            if
            (
                patchValues[patchFaceI].size()
             != faceQuadPoints[faceI].size()
            )
            {
                FatalErrorInFunction
                    << "Patch quadrature-value count does not match face "
                    << faceI << abort(FatalError);
            }

            prescribedValues[faceI] = patchValues[patchFaceI];
        }
    }

    // Evaluate weighted one-sided reconstructions on fixed-value faces.
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (!includePatchInStencils_[patchI])
        {
            continue;
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;
            const labelUList& faceStencil = stencils[faceI];
            auto& faceCoefficients = coeffs[faceI];
            const List<FixedList<label, 2>>& dataStencil =
                boundaryDataStencil[faceI];
            auto& dataCoefficients =
                boundaryDataCoeffs[faceI];

            forAll(faceQuadPoints[faceI], qpI)
            {
                const UList<vector>& qpCoefficients =
                    faceCoefficients[qpI];
                const UList<vector>& qpDataCoefficients =
                    dataCoefficients[qpI];

                forAll(faceStencil, stencilI)
                {
                    result[faceI][qpI] +=
                        qpCoefficients[stencilI]
                       *fieldValue
                        (
                            faceStencil[stencilI],
                            globalCells,
                            remoteLoc,
                            vfI,
                            remoteField
                        );
                }

                forAll(dataStencil, dataI)
                {
                    const FixedList<label, 2>& address = dataStencil[dataI];
                    result[faceI][qpI] +=
                        qpDataCoefficients[dataI]
                       *prescribedValues[address[0]][address[1]];
                }
            }
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

    auto& stencils =
        compactListListCRef(stencil().cellsStencil());
    auto& gradCoeffs =
        compactListListCRef(this->cellGradCoeffs());

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

    GeometricField<GradType, fvPatchField, volMesh>& grad = tmpRef(tGrad);

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
