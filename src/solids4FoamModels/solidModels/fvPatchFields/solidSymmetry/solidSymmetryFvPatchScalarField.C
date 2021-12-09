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

#include "solidSymmetryFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    symmetryFvPatchField<scalar>(p, iF)
{}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const solidSymmetryFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    symmetryFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (!isType<symmetryFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "solidSymmetryFvPatchScalarField::"
            "solidSymmetryFvPatchScalarField\n"
            "(\n"
            "    const solidSymmetryFvPatchScalarField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
#ifdef OPENFOAMESIORFOUNDATION
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
#else
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
#endif
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    symmetryFvPatchField<scalar>(p, iF, dict)
{
    Info << "Symmetry boundary condition with non-orthogonal correction"
        << endl;

    if (!isType<symmetryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "solidSymmetryFvPatchScalarField::"
            "solidSymmetryFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<scalar>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
#ifdef OPENFOAMESIORFOUNDATION
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
#else
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
#endif
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const solidSymmetryFvPatchScalarField& ptf
)
:
    symmetryFvPatchField<scalar>(ptf)
{}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const solidSymmetryFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    symmetryFvPatchField<scalar>(ptf, iF)
{}


// return gradient at boundary
tmp<Field<scalar> > solidSymmetryFvPatchScalarField::snGrad() const
{
    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    // Lookup the gradient of displacement field
    const fvPatchField<vector>& gradU =
        patch().lookupPatchField<volVectorField, vector>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Calculate the corrected patch internal field
    scalarField UP(patchInternalField());
    UP += (k & gradU.patchInternalField());

    return
    (
        transform(I - 2.0*sqr(n), UP)
      - UP
    )*(patch().deltaCoeffs()/2.0);
}


// Evaluate the field on the patch
void solidSymmetryFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    // Lookup the gradient of displacement field
    const fvPatchField<vector>& gradU =
        patch().lookupPatchField<volVectorField, vector>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Calculate the corrected patch internal field
    scalarField UP(patchInternalField());
    UP += (k & gradU.patchInternalField());

    Field<scalar>::operator=
    (
        (
            UP
          + transform(I - 2.0*sqr(n), UP)
        )/2.0
    );
}


// Write
void solidSymmetryFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, solidSymmetryFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
