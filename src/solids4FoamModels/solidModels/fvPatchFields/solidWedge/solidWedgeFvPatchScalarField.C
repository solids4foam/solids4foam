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

#include "solidWedgeFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wedgeFvPatchField<scalar>(p, iF)
{}


solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField
(
    const solidWedgeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wedgeFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField\n"
            "(\n"
            "    const solidWedgeFvPatchScalarField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
#ifdef OPENFOAM_NOT_EXTEND
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
#else
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
#endif
            << exit(FatalIOError);
    }
}


solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wedgeFvPatchField<scalar>(p, iF, dict)
{
    if (!isType<wedgeFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<scalar>& field,\n"
            "    dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
#ifdef OPENFOAM_NOT_EXTEND
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
#else
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
#endif
            << exit(FatalIOError);
    }

    //evaluate();

    // to initialise, we evaluate without non-orthogonal correction
    // this is because the explicit gradient is need for the correction
    // and obviously cannot be created until the field is created
    fvPatchField<scalar>::operator==
    (
        transform
        (
            refCast<const wedgeFvPatch>(this->patch()).faceT(),
            this->patchInternalField()
        )
    );
}

#ifndef OPENFOAM_ORG
solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField
(
    const solidWedgeFvPatchScalarField& ptf
)
:
    wedgeFvPatchField<scalar>(ptf)
{}
#endif

solidWedgeFvPatchScalarField::solidWedgeFvPatchScalarField
(
    const solidWedgeFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wedgeFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<scalar> > solidWedgeFvPatchScalarField::snGrad() const
{
    //Info<< "solidWedgeFvPatchScalarField::snGrad()" << endl;
    // Method:
    // We project the patch normal back to the wedge centrePlane and find the
    // intersection point. We then extrapolate U from the Cn to this
    // intersection point (projC) using gradU looked up from the solver.
    // We can then calculate snGrad by getting the difference between projU and
    // transformed projU (i.e. across the wedge), and then dividing by the
    // magnitude of the distance between them.

    const wedgePolyPatch& wedgePatch =
        refCast<const wedgePolyPatch>(patch().patch());

    const vectorField patchC(patch().patch().faceCentres());
    const vectorField n(patch().nf());
    const vector& centreN = wedgePatch.centreNormal();
    const scalarField d(((patch().Cn() - patchC) & centreN)/(n & centreN));
    const vectorField projC(d*n + patchC);

    // Calculate correction vector which connects actual cell centre to the
    // transformed cell centre
    const vectorField k(projC - patch().Cn());

    const Field<scalar> pif(this->patchInternalField());

    const fvPatchField<vector>& gradU =
        patch().lookupPatchField<volVectorField, vector>
        (
#ifdef OPENFOAM_NOT_EXTEND
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    const Field<scalar> projU
    (
        this->patchInternalField() + (k & gradU.patchInternalField())
    );

    // Calculate delta coeffs from proj position on centre plane to transformed
    // projected position
    const scalarField projDeltaCoeff
    (
        1.0/mag(transform(wedgePatch.cellT(), projC) - projC)
    );

    return
    (
        transform(wedgePatch.cellT(), projU) - projU
    )*projDeltaCoeff;

    // old way without correction
    // return
    // (
    //  transform(refCast<const wedgeFvPatch>(this->patch()).cellT(), pif) - pif
    // )*(0.5*this->patch().deltaCoeffs());
}


void solidWedgeFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Standard non-orthogonal correction vectors are parallel to the patch, but
    // for axisymmetric cases we must use a correction vector that lies along
    // the centre plane. Then we can transform/rotate this centre plane value to
    // the wedge patch.
    const wedgeFvPatch& wedgePatch =
        refCast<const wedgeFvPatch>(this->patch());

    // Rotate patchC field back to centre plane to find transformed cell centres
    const vectorField patchC(patch().patch().faceCentres());
    const vectorField transC(wedgePatch.faceT().T() & patchC);

    // Calculate correction vector which connects actual cell centre to the
    // transformed cell centre
    const vectorField k(transC - patch().Cn());

    const fvPatchField<vector>& gradU =
        patch().lookupPatchField<volVectorField, vector>
        (
#ifdef OPENFOAM_NOT_EXTEND
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    Field<scalar> pif(this->patchInternalField());
    pif += (k & gradU.patchInternalField());

    Field<scalar>::operator=
    (
        transform(wedgePatch.faceT(), pif)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, solidWedgeFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
