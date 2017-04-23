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

#include "normalDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

normalDisplacementFvPatchVectorField::normalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    normalDisp_(p.size(), 0.0),
    dispSeries_()
{}


normalDisplacementFvPatchVectorField::normalDisplacementFvPatchVectorField
(
    const normalDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    normalDisp_(ptf.normalDisp_, mapper),
    dispSeries_(ptf.dispSeries_)
{}


normalDisplacementFvPatchVectorField::normalDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    normalDisp_(p.size(), 0.0),
    dispSeries_()
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    normal displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<scalar>(dict.subDict("displacementSeries"));

        fvPatchField<vector>::operator==
        (
	    patch().nf()*dispSeries_(this->db().time().timeOutputValue())
        );
    }
    else
    {
        normalDisp_ = scalarField("normalDisplacement", dict, p.size());

        fvPatchField<vector>::operator==
        (
	    patch().nf()*normalDisp_ 
        );
    }
}


normalDisplacementFvPatchVectorField::normalDisplacementFvPatchVectorField
(
    const normalDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    normalDisp_(pivpvf.normalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


normalDisplacementFvPatchVectorField::normalDisplacementFvPatchVectorField
(
    const normalDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    normalDisp_(pivpvf.normalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void normalDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    normalDisp_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void normalDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const normalDisplacementFvPatchVectorField& dmptf =
        refCast<const normalDisplacementFvPatchVectorField>(ptf);

    normalDisp_.rmap(dmptf.normalDisp_, addr);
}


void normalDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField nDisp = normalDisp_;

    if (dispSeries_.size())
    {
        nDisp = dispSeries_(this->db().time().timeOutputValue());
    }

    vectorField disp = nDisp*patch().nf();

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
normalDisplacementFvPatchVectorField::snGrad() const
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
normalDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
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

void normalDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        normalDisp_.writeEntry("normalDisp", os);
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    normalDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
