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

#include "neoHookeanMMSFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

neoHookeanMMSFvPatchVectorField::neoHookeanMMSFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    mmsPtr_()
{}


neoHookeanMMSFvPatchVectorField::neoHookeanMMSFvPatchVectorField
(
    const neoHookeanMMSFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    mmsPtr_()
{}


neoHookeanMMSFvPatchVectorField::neoHookeanMMSFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    mmsPtr_()
{}


#ifndef OPENFOAM_ORG
neoHookeanMMSFvPatchVectorField::neoHookeanMMSFvPatchVectorField
(
    const neoHookeanMMSFvPatchVectorField& pivpvf
)
:
    fixedDisplacementFvPatchVectorField(pivpvf),
    mmsPtr_()
{}
#endif


neoHookeanMMSFvPatchVectorField::neoHookeanMMSFvPatchVectorField
(
    const neoHookeanMMSFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(pivpvf, iF),
    mmsPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void neoHookeanMMSFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (!mmsPtr_.valid())
    {
        mmsPtr_.reset(new neoHookeanMMS(patch().boundaryMesh().mesh()));
    }

    // Set the displacement at each patch face centre at the current time
    const scalar t = this->db().time().value();
    vectorField& disp = totalDisp();
    const vectorField& Cf = patch().Cf();
    forAll(disp, faceI)
    {
        disp[faceI] = mmsPtr_->displacement(Cf[faceI], t);
    }

    fixedDisplacementFvPatchVectorField::updateCoeffs();
}


void neoHookeanMMSFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    neoHookeanMMSFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
