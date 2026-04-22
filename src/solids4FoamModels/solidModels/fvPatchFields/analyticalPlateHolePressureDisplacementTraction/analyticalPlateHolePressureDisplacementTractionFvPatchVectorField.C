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

#include "analyticalPlateHolePressureDisplacementTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "plateHoleAnalyticalFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::
analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionPressureDisplacementFvPatchVectorField(p, iF),
    T_(0.0),
    holeR_(0.0)
{}


analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::
analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionPressureDisplacementFvPatchVectorField(p, iF, dict),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius")))
{}


analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::
analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
(
    const analyticalPlateHolePressureDisplacementTractionFvPatchVectorField& ptpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionPressureDisplacementFvPatchVectorField(ptpvf, p, iF, mapper),
    T_(ptpvf.T_),
    holeR_(ptpvf.holeR_)
{}


analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::
analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
(
    const analyticalPlateHolePressureDisplacementTractionFvPatchVectorField& ptpvf
)
:
    tractionPressureDisplacementFvPatchVectorField(ptpvf),
    T_(ptpvf.T_),
    holeR_(ptpvf.holeR_)
{}


analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::
analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
(
    const analyticalPlateHolePressureDisplacementTractionFvPatchVectorField& ptpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionPressureDisplacementFvPatchVectorField(ptpvf, iF),
    T_(ptpvf.T_),
    holeR_(ptpvf.holeR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField& trac = traction();

    const vectorField n(patch().nf());
    const vectorField& Cf = patch().Cf();

    forAll(trac, faceI)
    {
        const vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);

        trac[faceI] =
            plateHoleAnalyticalFields::traction(curC, n[faceI], T_, holeR_);
    }

    tractionPressureDisplacementFvPatchVectorField::updateCoeffs();
}


void analyticalPlateHolePressureDisplacementTractionFvPatchVectorField::write
(
    Ostream& os
) const
{
    tractionPressureDisplacementFvPatchVectorField::write(os);

    os.writeKeyword("farFieldTractionX")
        << T_ << token::END_STATEMENT << nl;

    os.writeKeyword("holeRadius")
        << holeR_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    analyticalPlateHolePressureDisplacementTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
