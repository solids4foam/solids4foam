/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "momentBeamRotationNRFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupBeamModel.H"
#include "pseudoVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

momentBeamRotationNRFvPatchVectorField::
momentBeamRotationNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    moment_(p.size(), vector::zero),
    momentSeries_()
{
    fvPatchVectorField::operator==(patchInternalField());
}


momentBeamRotationNRFvPatchVectorField::
momentBeamRotationNRFvPatchVectorField
(
    const momentBeamRotationNRFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(tdpvf, p, iF, mapper),
    moment_(tdpvf.moment_, mapper),
    momentSeries_(tdpvf.momentSeries_)
{}


momentBeamRotationNRFvPatchVectorField::
momentBeamRotationNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    moment_(p.size(), vector::zero),
    momentSeries_()
{
    if (dict.found("value"))
    {
        fvPatchVectorField::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator==
        (
            patchInternalField()
        );
    }

    // Check if traction is time-varying
    if (dict.found("momentSeries"))
    {
        Info<< "moment is time-varying" << endl;
        momentSeries_ =
            interpolationTable<vector>(dict.subDict("momentSeries"));
    }
    else
    {
        moment_ = vectorField("moment", dict, p.size());
	//Info << "Moment value: " << moment_ << endl;
    }

}


momentBeamRotationNRFvPatchVectorField::
momentBeamRotationNRFvPatchVectorField
(
    const momentBeamRotationNRFvPatchVectorField& tdpvf
)
:
    fixedValueFvPatchVectorField(tdpvf),
    moment_(tdpvf.moment_),
    momentSeries_(tdpvf.momentSeries_)
{}


momentBeamRotationNRFvPatchVectorField::
momentBeamRotationNRFvPatchVectorField
(
    const momentBeamRotationNRFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tdpvf, iF),
    moment_(tdpvf.moment_),
    momentSeries_(tdpvf.momentSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void momentBeamRotationNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    moment_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void momentBeamRotationNRFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const momentBeamRotationNRFvPatchVectorField& dmptf =
        refCast<const momentBeamRotationNRFvPatchVectorField>(ptf);

    moment_.rmap(dmptf.moment_, addr);
}


// Update the coefficients associated with the patch field
void momentBeamRotationNRFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (momentSeries_.size())
    {
        moment_ = momentSeries_(this->db().time().timeOutputValue());
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void momentBeamRotationNRFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    fvPatchField<vector>& DTheta =
        const_cast<fvPatchField<vector>& >
        (
            patch().lookupPatchField<volVectorField, vector>("DTheta")
        );

    const fvPatchField<tensor>& pLambda =
        patch().lookupPatchField<volTensorField, tensor>("Lambda");

    const tensorField& CMTheta =
        patch().lookupPatchField<surfaceTensorField, tensor>("CMTheta");

    const tensorField& CMTheta2 =
        patch().lookupPatchField<surfaceTensorField, tensor>("CMTheta2");

    const vectorField& explicitM =
        patch().lookupPatchField<surfaceVectorField, vector>("explicitM");

    const scalarField delta (1.0/patch().deltaCoeffs());

    const tensorField invA (inv(CMTheta/delta + CMTheta2));

    const vectorField newDTheta
	(
        (invA & (moment_ - explicitM)) +
        ((invA & (CMTheta/delta)) & DTheta.patchInternalField())
	);

    fixedValueFvPatchField<vector>::operator==
    (
        (*this)
      + (pLambda.T() & newDTheta)
    );

    DTheta = newDTheta;

    fixedValueFvPatchVectorField::evaluate();
}


// Write
void momentBeamRotationNRFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    if (momentSeries_.size())
    {
        os.writeKeyword("momentSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        momentSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    moment_.writeEntry("moment", os);

    // writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, momentBeamRotationNRFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
