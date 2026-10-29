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

#include "fixedBeamDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
// #include "fvcMeshPhi.H"
// #include "pointMesh.H"
// #include "pointFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedBeamDisplacementFvPatchVectorField::fixedBeamDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_()
{}


fixedBeamDisplacementFvPatchVectorField::fixedBeamDisplacementFvPatchVectorField
(
    const fixedBeamDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    totalDisp_(ptf.totalDisp_, mapper),
    dispSeries_(ptf.dispSeries_)
{}


fixedBeamDisplacementFvPatchVectorField::fixedBeamDisplacementFvPatchVectorField
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


fixedBeamDisplacementFvPatchVectorField::fixedBeamDisplacementFvPatchVectorField
(
    const fixedBeamDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    totalDisp_(pivpvf.totalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


fixedBeamDisplacementFvPatchVectorField::fixedBeamDisplacementFvPatchVectorField
(
    const fixedBeamDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    totalDisp_(pivpvf.totalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedBeamDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    totalDisp_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedBeamDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const fixedBeamDisplacementFvPatchVectorField& dmptf =
        refCast<const fixedBeamDisplacementFvPatchVectorField>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedBeamDisplacementFvPatchVectorField::updateCoeffs()
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

    // Info << patch().name() << ", " << disp << endl;

    fvPatchField<vector>::operator==(disp);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void fixedBeamDisplacementFvPatchVectorField::write(Ostream& os) const
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
    fixedBeamDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
