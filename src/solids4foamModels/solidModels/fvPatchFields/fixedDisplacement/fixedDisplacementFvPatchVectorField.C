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

#include "fixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_()
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    totalDisp_(ptf.totalDisp_, mapper),
    dispSeries_(ptf.dispSeries_)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    totalDisp_("value", dict, p.size()),
    dispSeries_()
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        fvPatchField<vector>::operator==
        (
            dispSeries_(this->db().time().timeOutputValue())
        );
    }
}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    totalDisp_(pivpvf.totalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    totalDisp_(pivpvf.totalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    totalDisp_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const fixedDisplacementFvPatchVectorField& dmptf =
        refCast<const fixedDisplacementFvPatchVectorField>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp = totalDisp_;

    if (dispSeries_.size())
    {
        disp = dispSeries_(this->db().time().timeOutputValue());
    }

    bool incremental = bool(dimensionedInternalField().name() == "DU");

    if (incremental)
    {
        if (patch().boundaryMesh().mesh().foundObject<volVectorField>("U_0"))
        {
            const fvPatchField<vector>& Uold =
              patch().lookupPatchField<volVectorField, vector>("U_0");

            disp -= Uold;
        }
        else if (patch().boundaryMesh().mesh().foundObject<volVectorField>("U"))
        {
            const fvPatchField<vector>& Uold =
              patch().lookupPatchField<volVectorField, vector>("U");

            disp -= Uold;
        }
    }

    fvPatchField<vector>::operator==
    (
        disp
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> >
fixedDisplacementFvPatchVectorField::snGrad() const
{
    //- fixedValue snGrad with no correction
    //  return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();

    //- correction vector
    vectorField k = delta - n*(n&delta);

    return
    (
        *this
        - (patchInternalField() + (k & gradField.patchInternalField()))
    )*this->patch().deltaCoeffs();
}

tmp<Field<vector> >
fixedDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();

    //- correction vector
    vectorField k = delta - n*(n&delta);

    return
    (
        this->patch().deltaCoeffs()
       *(*this - (k & gradField.patchInternalField()))
    );
}

void fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
