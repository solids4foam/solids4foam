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

#include "movingLeastSquares.H"

namespace Foam
{

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

template<class Type, class GlobalCells>
inline Type movingLeastSquares::fieldValue
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
Type movingLeastSquares::evaluateAtPoint
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

    // Construct all coefficient tables collectively before checking for the
    // negative cell marker used by parallel point evaluation.
    (void)cellGradCoeffs();
    if (polynomialOrder() >= 2)
    {
        (void)cellSecondGradCoeffs();
    }
    if (polynomialOrder() >= 3)
    {
        (void)cellThirdGradCoeffs();
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
autoPtr<CompactListList<Type>> movingLeastSquares::patchFaceQuadValues
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const label patchI
) const
{
    FatalErrorInFunction
        << "patchFaceQuadValues is not implemented for field type "
        << pTraits<Type>::typeName
        << " on patch " << patchI
        << abort(FatalError);

    return autoPtr<CompactListList<Type>>(nullptr);
}


template<class Type>
tmp<GeometricField<
    typename outerProduct<vector, Type>::type,
    fvPatchField,
    volMesh
>> movingLeastSquares::grad
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

        // Stencil contribution
        forAll(stencil, cI)
        {
            grad[cellI] +=
                coeffs[cI]
              * fieldValue
                (
                    stencil[cI],
                    globalCells,
                    remoteLoc,
                    vfI,
                    remoteField
                );
        }
        // Cell-centre contribution (last entry)
        grad[cellI] += coeffs[stencil.size()] * vfI[cellI];
    }

    grad.correctBoundaryConditions();

    return tGrad;
}


template<class Type>
autoPtr<CompactListList<typename outerProduct<vector, Type>::type>>
movingLeastSquares::fGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    autoPtr<CompactListList<GradType>> tResult
    (
        new CompactListList<GradType>(quadrature().faceQuadPoints().sizes())
    );

    this->fGrad(vf, autoPtrRef(tResult));

    return tResult;
}


template<class Type>
void movingLeastSquares::fGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    CompactListList<typename outerProduct<vector, Type>::type>& result
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    // Get references
    const fvMesh& mesh = mesh_;
    const globalIndex& globalCells = stencil().globalCells();
    const Map<FixedList<label, 2>>& remoteLoc = stencil().remoteCellLocation();
    auto& faceQuadPts =
        compactListListCRef(quadrature().faceQuadPoints());
    auto& stencils =
        compactListListCRef(faceGradStencil());
#ifdef FOAMEXTEND
    List<CompactListList<vector>>& fGradCoeffs =
        const_cast<List<CompactListList<vector>>&>(this->faceGradCoeffs());
#else
    const List<CompactListList<vector>>& fGradCoeffs = this->faceGradCoeffs();
#endif

    // Get values from remote processors
    const Field<Type>& vfI = vf.internalField();
    const List<Field<Type>> remoteField = stencil().remoteFieldPerProc(vfI);

    forAll(result, faceI)
    {
        forAll(result[faceI], qpI)
        {
            result[faceI][qpI] = pTraits<GradType>::zero;
        }
    }

    // Loop over internal faces
    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        const UList<label>& stencil = stencils[faceI];
        auto& coeffs = fGradCoeffs[faceI];

        forAll(faceQuadPts[faceI], qpI)
        {
            forAll(stencil, cI)
            {
                result[faceI][qpI] +=
                    coeffs[qpI][cI]
                  * fieldValue
                    (
                        stencil[cI],
                        globalCells,
                        remoteLoc,
                        vfI,
                        remoteField
                    );
            }
        }
    }

    // Loop over boundary
    forAll(vf.boundaryField(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        // Empty patch
        if (isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        const label patchStart = pp.start();

        // Symmetry patch
        if
        (
            isA<symmetryPolyPatch>(pp)
#ifndef FOAMEXTEND
         || isA<symmetryPlanePolyPatch>(pp)
#endif
        )
        {
            const vectorField patchNormals(mesh.boundary()[patchI].nf());

            forAll(pp, faceI)
            {
                const label globalFaceID = patchStart + faceI;
                const auto stencil = stencils[globalFaceID];
                auto& coeffs = fGradCoeffs[globalFaceID];

                const label stencilSize = stencil.size();
                const vector& faceNormal = patchNormals[faceI];
                const tensor R = I - 2.0*sqr(faceNormal);

                forAll(faceQuadPts[globalFaceID], qpI)
                {
                    forAll(stencil, cI)
                    {
                        // Physical stencil
                        const Type val = fieldValue
                            (
                                stencil[cI],
                                globalCells,
                                remoteLoc,
                                vfI,
                                remoteField
                            );

                        result[globalFaceID][qpI] +=
                            coeffs[qpI][cI] * val;

                        // Add mirrored cell value from mirrored stencil
                        const Type mirrorVal = transform(R, val);

                        result[globalFaceID][qpI] +=
                            coeffs[qpI][cI + stencilSize] * mirrorVal;
                    }
                }
            }
        }
        // All other patches
        else
        {
            forAll(pp, faceI)
            {
                const label globalFaceID = patchStart + faceI;
                const auto stencil = stencils[globalFaceID];
                auto& coeffs = fGradCoeffs[globalFaceID];

                forAll(faceQuadPts[globalFaceID], qpI)
                {
                    forAll(stencil, cI)
                    {
                        result[globalFaceID][qpI] +=
                            coeffs[qpI][cI]
                          * fieldValue
                            (
                                stencil[cI],
                                globalCells,
                                remoteLoc,
                                vfI,
                                remoteField
                            );
                    }
                }
            }
        }

        // For fixed-value patches we need to add ghost point contribution
        // Ghost point is not in stencil!
        const fvPatchField<Type>& pf = vf.boundaryField()[patchI];

        if (pf.fixesValue())
        {
            if (!includePatchInStencils_[patchI])
            {
                FatalErrorInFunction
                    << "Patch " << patchI << " fixes value but it is not "
                    << "considered during coefficient calculations!"
                    << abort(FatalError);
            }

            // Values at patch faces quadrature points
            autoPtr<CompactListList<Type>> patchFaceQuadValsPtr =
                patchFaceQuadValues(vf, patchI);
#ifdef FOAMEXTEND
            CompactListList<Type>& quadVal = autoPtrRef(patchFaceQuadValsPtr);
#else
            const CompactListList<Type>& quadVal = patchFaceQuadValsPtr();
#endif

            forAll(pf, faceI)
            {
                const label globalFaceID = patchStart + faceI;
                const auto stencil = stencils[globalFaceID];
                auto& coeffs = fGradCoeffs[globalFaceID];
                const label ghostPointID = stencil.size();

                forAll(faceQuadPts[globalFaceID], qpI)
                {
                    result[globalFaceID][qpI] +=
                        coeffs[qpI][ghostPointID] * quadVal[faceI][qpI];
                }
            }
        }
    }
}

}

// ************************************************************************* //
