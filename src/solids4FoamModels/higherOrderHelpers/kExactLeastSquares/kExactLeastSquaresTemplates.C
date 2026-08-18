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
