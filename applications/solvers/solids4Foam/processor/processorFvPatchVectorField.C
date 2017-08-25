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

#include "processorFvPatchVectorField.H"

#include "solidPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<>
void processorFvPatchField<vector>::initInterfaceMatrixUpdate
(
    const Field<vector>& psiInternal,
    Field<vector>&,
    const BlockLduMatrix<vector>&,
    const CoeffField<vector>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    procPatch_.compressedSend
    (
        commsType,
        this->patch().patchInternalField(psiInternal)()
    );
}


template<>
void processorFvPatchField<vector>::updateInterfaceMatrix
(
    const Field<vector>& psiInternal,
    Field<vector>& result,
    const BlockLduMatrix<vector>&,
    const CoeffField<vector>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    vectorField pnf
    (
        procPatch_.compressedReceive<vector>(commsType, this->size())()
    );

    // Transform according to the transformation tensor
    // PC: disabled: check if it affects standard procs!
    //transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& faceCells = this->patch().faceCells();
    const tensorField& tensorCoeffs = coeffs.asSquare();

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= tensorCoeffs[elemI] & pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += tensorCoeffs[elemI] & pnf[elemI];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
