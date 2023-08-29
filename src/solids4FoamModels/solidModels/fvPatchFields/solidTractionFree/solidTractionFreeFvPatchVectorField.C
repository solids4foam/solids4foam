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

#include "solidTractionFreeFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
    traction() = vector::zero;
    pressure() = 0.0;
}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF)
{
    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    traction() = vector::zero;
    pressure() = 0.0;
}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper)
{
    traction() = vector::zero;
    pressure() = 0.0;
}

#ifndef OPENFOAMFOUNDATION
solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf
)
:
    solidTractionFvPatchVectorField(stpvf)
{
    traction() = vector::zero;
    pressure() = 0.0;
}
#endif

solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF)
{
    traction() = vector::zero;
    pressure() = 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFreeFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
