/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "symmetryDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

symmetryDisplacementFvPatchVectorField::symmetryDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(p, iF),
    secondOrder_(false)
{}


symmetryDisplacementFvPatchVectorField::symmetryDisplacementFvPatchVectorField
(
    const symmetryDisplacementFvPatchVectorField& ptf,
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
            "symmetryDisplacementFvPatchVectorField::"
            "symmetryDisplacementFvPatchVectorField\n"
            "(\n"
            "    const symmetryDisplacementFvPatchVectorField& ptf,\n"
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


symmetryDisplacementFvPatchVectorField::symmetryDisplacementFvPatchVectorField
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
        Info << "Second order correction: " << secondOrder_ << endl;
    }

    if (!isType<symmetryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "symmetryDisplacementFvPatchVectorField::"
            "symmetryDisplacementFvPatchVectorField\n"
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


symmetryDisplacementFvPatchVectorField::symmetryDisplacementFvPatchVectorField
(
    const symmetryDisplacementFvPatchVectorField& ptf
)
:
    symmetryFvPatchField<vector>(ptf),
    secondOrder_(ptf.secondOrder_)
{}


symmetryDisplacementFvPatchVectorField::symmetryDisplacementFvPatchVectorField
(
    const symmetryDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    symmetryFvPatchField<vector>(ptf, iF),
    secondOrder_(ptf.secondOrder_)
{}


// return gradient at boundary
tmp<Field<vector> > symmetryDisplacementFvPatchVectorField::snGrad() const
{
    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    word DName = this->dimensionedInternalField().name();

//     const fvMesh& msh = patch().boundaryMesh().mesh();
//     if (msh.found("quadraticReconstruction"))
//     {
//         Info << "symm-gadratic-recon" << endl;

//         const quadraticReconstruction& recon =
//             msh.lookupObject<quadraticReconstruction>
//             (
//                 "quadraticReconstruction"
//             );

//         const unallocLabelList& faceCells = patch().faceCells();
//         const volVectorField& U = msh.lookupObject<volVectorField>(UName);

//         vectorField UP =
//             this->patchInternalField()
//           + recon.difference(U, faceCells, k);

//         vectorField nGradUP = recon.derivative(U, faceCells, k, nHat);

//         return
//           2*(
//                 transform(I - 2.0*sqr(nHat), UP) - UP
//             )*(this->patch().deltaCoeffs()/2.0)
//           - transform(sqr(nHat), nGradUP);
//     }

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DName + ")"
        );

    vectorField DP = this->patchInternalField();
    DP += (k&gradD.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradDP = (nHat&gradD.patchInternalField());

        return
          2*(
                transform(I - 2.0*sqr(nHat), DP) - DP
            )*(this->patch().deltaCoeffs()/2.0)
          - transform(sqr(nHat), nGradDP);
    }
    else
    {
        return
        (
            transform(I - 2.0*sqr(nHat), DP)
          - DP
        )*(this->patch().deltaCoeffs()/2.0);
    }
}


// Evaluate the field on the patch
void symmetryDisplacementFvPatchVectorField::
evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    vectorField delta = patch().delta();
    vectorField k = delta - nHat*(nHat&delta);

    word DName = this->dimensionedInternalField().name();

//     const fvMesh& msh = patch().boundaryMesh().mesh();
//     if (msh.found("quadraticReconstruction"))
//     {
//         Info << "symm-gadratic-recon" << endl;

//         const quadraticReconstruction& recon =
//             msh.lookupObject<quadraticReconstruction>
//             (
//                 "quadraticReconstruction"
//             );

//         const unallocLabelList& faceCells = patch().faceCells();
//         const volVectorField& U = msh.lookupObject<volVectorField>(UName);

//         vectorField UP =
//             this->patchInternalField()
//           + recon.difference(U, faceCells, k);

//         vectorField nGradUP = recon.derivative(U, faceCells, k, nHat);

//         Field<vector>::operator=
//         (
//             transform
//             (
//                 I - sqr(nHat),
//                 UP + 0.5*nGradUP/this->patch().deltaCoeffs()
//             )
//         );
//     }


    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DName + ")"
        );

    vectorField DP = this->patchInternalField();
    DP += (k&gradD.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradDP = (nHat&gradD.patchInternalField());

        Field<vector>::operator=
        (
            transform
            (
                I - sqr(nHat),
                DP + 0.5*nGradDP/this->patch().deltaCoeffs()
            )
        );
    }
    else
    {
        Field<vector>::operator=
        (
            (
                DP
              + transform(I - 2.0*sqr(nHat), DP)
            )/2.0
        );
    }

    transformFvPatchField<vector>::evaluate();
}


// Write
void symmetryDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, symmetryDisplacementFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
