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

#include "solidSymmetryFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(p, iF)
{}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    symmetryFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (!isType<symmetryFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "solidSymmetryFvPatchVectorField::"
            "solidSymmetryFvPatchVectorField\n"
            "(\n"
            "    const solidSymmetryFvPatchVectorField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    symmetryFvPatchField<vector>(p, iF, dict)
{
    Info << "Symmetry boundary condition with non-orthogonal correction"
        << endl;

    if (!isType<symmetryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "solidSymmetryFvPatchVectorField::"
            "solidSymmetryFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<vector>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf
)
:
    symmetryFvPatchField<vector>(ptf)
{}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(ptf, iF)
{}


// return gradient at boundary
tmp<Field<vector> > solidSymmetryFvPatchVectorField::snGrad() const
{
    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        );

    vectorField UP = this->patchInternalField();
    UP += (k&gradU.patchInternalField());

    return
    (
        transform(I - 2.0*sqr(nHat), UP)
      - UP
    )*(this->patch().deltaCoeffs()/2.0);
}


// Evaluate the field on the patch
void solidSymmetryFvPatchVectorField::
evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        );

    vectorField UP = this->patchInternalField();
    UP += (k&gradU.patchInternalField());

    Field<vector>::operator=
    (
        (
            UP
          + transform(I - 2.0*sqr(nHat), UP)
        )/2.0
    );
}


// Write
void solidSymmetryFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidSymmetryFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
