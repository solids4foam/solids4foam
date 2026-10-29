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

#include "forceBeamDisplacementNRFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupBeamModel.H"
#include "pseudoVector.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

forceBeamDisplacementNRFvPatchVectorField::
forceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    curForce_(p.size(), vector::zero),
    forceSeries_(),
    transitionSeries_()
{
    // fixedValueFvPatchVectorField::operator==(patchInternalField());
    // gradient() = vector::zero;
}


forceBeamDisplacementNRFvPatchVectorField::
forceBeamDisplacementNRFvPatchVectorField
(
    const forceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(tdpvf, p, iF, mapper),
    force_(tdpvf.force_, mapper),
    curForce_(tdpvf.curForce_, mapper),
    forceSeries_(tdpvf.forceSeries_),
    transitionSeries_(tdpvf.transitionSeries_)
{}


forceBeamDisplacementNRFvPatchVectorField::
forceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    curForce_(p.size(), vector::zero),
    forceSeries_(),
    transitionSeries_()
{
    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchVectorField::operator==
        (
            patchInternalField()
        );
    }

    // Check if force is time-varying
    if (dict.found("transitionSeries"))
    {
        Info<< "force is time-varying using transition function" << endl;
        transitionSeries_ =
            interpolationTable<scalar>(dict.subDict("transitionSeries"));

        force_ = vectorField("force", dict, p.size());

        curForce_ =
            transitionSeries_(this->db().time().timeOutputValue())
           *force_;
    }
    else if (dict.found("forceSeries"))
    {
        Info<< "force is time-varying" << endl;
        forceSeries_ =
            interpolationTable<vector>(dict.subDict("forceSeries"));

        curForce_ = forceSeries_(this->db().time().timeOutputValue());
        force_ = curForce_;
    }
    else
    {
        force_ = vectorField("force", dict, p.size());
        curForce_ = force_;
    }
}


forceBeamDisplacementNRFvPatchVectorField::
forceBeamDisplacementNRFvPatchVectorField
(
    const forceBeamDisplacementNRFvPatchVectorField& tdpvf
)
:
    fixedValueFvPatchVectorField(tdpvf),
    force_(tdpvf.force_),
    curForce_(tdpvf.curForce_),
    forceSeries_(tdpvf.forceSeries_),
    transitionSeries_(tdpvf.transitionSeries_)
{}


forceBeamDisplacementNRFvPatchVectorField::
forceBeamDisplacementNRFvPatchVectorField
(
    const forceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tdpvf, iF),
    force_(tdpvf.force_),
    curForce_(tdpvf.curForce_),
    forceSeries_(tdpvf.forceSeries_),
    transitionSeries_(tdpvf.transitionSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forceBeamDisplacementNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    force_.autoMap(m);
    curForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void forceBeamDisplacementNRFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const forceBeamDisplacementNRFvPatchVectorField& dmptf =
        refCast<const forceBeamDisplacementNRFvPatchVectorField>(ptf);

    force_.rmap(dmptf.force_, addr);
    curForce_.rmap(dmptf.force_, addr);
}


// Update the coefficients associated with the patch field
void forceBeamDisplacementNRFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (forceSeries_.size())
    {
        curForce_ = forceSeries_(this->db().time().timeOutputValue());
        force_ = curForce_;

        // Info << "force: " << force_ << endl;
    }
    else if(transitionSeries_.size())
    {
        curForce_ =
            transitionSeries_(this->db().time().timeOutputValue())
           *force_;
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void forceBeamDisplacementNRFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchField<vector>& DW =
        const_cast<fvPatchField<vector>& >
        (
            patch().lookupPatchField<volVectorField, vector>("DW")
        );

    const fvPatchField<vector>& DTheta =
        patch().lookupPatchField<volVectorField, vector>("DTheta");

    const tensorField& CQW =
        patch().lookupPatchField<surfaceTensorField, tensor>("CQW");

    const tensorField& CQTheta =
        patch().lookupPatchField<surfaceTensorField, tensor>("CQTheta");

    const vectorField explicitQ =
        patch().lookupPatchField<surfaceVectorField, vector>("explicitQ");

    const scalarField delta(1.0/patch().deltaCoeffs());

    const tensorField invA(inv(CQW/delta));

    const vectorField newDW
    (
        (invA & (force() - explicitQ)) -
        (invA & (CQTheta & DTheta))
      + DW.patchInternalField()
    );

    DW = newDW;

    // Setting W_ field
    fixedValueFvPatchField<vector>::operator==((*this) + newDW);

    fixedValueFvPatchVectorField::evaluate();
}


// Write
void forceBeamDisplacementNRFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    if (forceSeries_.size())
    {
        os.writeKeyword("forceSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        forceSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    curForce_.writeEntry("force", os);

    // writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, forceBeamDisplacementNRFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
