/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    symmetryFvPatchField<vector>(p, iF),
    secondOrder_(false)
{}


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    symmetryFvPatchField<vector>(ptf, p, iF, mapper),
    secondOrder_(ptf.secondOrder_)
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


solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    symmetryFvPatchField<vector>(p, iF, dict),
    secondOrder_(false)
{
    Info << "Symmetry boundary condition with non-orthogonal correction"
        << endl;

    if (dict.found("secondOrder"))
    {
        secondOrder_ = Switch(dict.lookup("secondOrder"));
        Info<< "Second order correction: " << secondOrder_ << endl;
    }

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

#ifndef OPENFOAMFOUNDATION
solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf
)
:
    symmetryFvPatchField<vector>(ptf),
    secondOrder_(ptf.secondOrder_)
{}
#endif

solidSymmetryFvPatchVectorField::solidSymmetryFvPatchVectorField
(
    const solidSymmetryFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(ptf, iF),
    secondOrder_(ptf.secondOrder_)
{}


// return gradient at boundary
tmp<Field<vector> > solidSymmetryFvPatchVectorField::snGrad() const
{
    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    // Lookup the gradient of displacement field
    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Calculate the corrected patch internal field
    const vectorField DP
    (
        patchInternalField()
      + (k & gradD.patchInternalField())
    );

    if (secondOrder_)
    {
        // Normal component of patch internal gradient
        const vectorField nGradDP(n & gradD.patchInternalField());

        return
          2*(
                transform(I - 2.0*sqr(n), DP) - DP
            )*(patch().deltaCoeffs()/2.0)
          - transform(sqr(n), nGradDP);
    }
    else
    {
        return
        (
            transform(I - 2.0*sqr(n), DP)
          - DP
        )*(patch().deltaCoeffs()/2.0);
    }
}


// Evaluate the field on the patch
void solidSymmetryFvPatchVectorField::
evaluate(const Pstream::commsTypes)
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
    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Calculate the corrected patch internal field
    const vectorField DP
    (
        patchInternalField()
      + (k & gradD.patchInternalField())
    );

    if (secondOrder_)
    {
        const vectorField nGradDP(n&gradD.patchInternalField());

        Field<vector>::operator=
        (
            transform
            (
                I - sqr(n),
                DP + 0.5*nGradDP/patch().deltaCoeffs()
            )
        );
    }
    else
    {
        Field<vector>::operator=
        (
            (
                DP
              + transform(I - 2.0*sqr(n), DP)
            )/2.0
        );
    }

    transformFvPatchField<vector>::evaluate();
}


// Write
void solidSymmetryFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidSymmetryFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
