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

#include "neoHookeanManufacturedSolutionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

neoHookeanManufacturedSolutionFvPatchVectorField::
neoHookeanManufacturedSolutionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    mmsPtr_()
{}


neoHookeanManufacturedSolutionFvPatchVectorField::
neoHookeanManufacturedSolutionFvPatchVectorField
(
    const neoHookeanManufacturedSolutionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    mmsPtr_()
{}


neoHookeanManufacturedSolutionFvPatchVectorField::
neoHookeanManufacturedSolutionFvPatchVectorField
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
neoHookeanManufacturedSolutionFvPatchVectorField::
neoHookeanManufacturedSolutionFvPatchVectorField
(
    const neoHookeanManufacturedSolutionFvPatchVectorField& pivpvf
)
:
    fixedDisplacementFvPatchVectorField(pivpvf),
    mmsPtr_()
{}
#endif


neoHookeanManufacturedSolutionFvPatchVectorField::
neoHookeanManufacturedSolutionFvPatchVectorField
(
    const neoHookeanManufacturedSolutionFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(pivpvf, iF),
    mmsPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void neoHookeanManufacturedSolutionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (!mmsPtr_.valid())
    {
        mmsPtr_.reset
        (
            new neoHookeanManufacturedSolution(patch().boundaryMesh().mesh())
        );
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


autoPtr<CompactListList<vector>>
neoHookeanManufacturedSolutionFvPatchVectorField::evaluateQuadrature() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const solidModel& solMod = lookupSolidModel(mesh);

    // Quadrature points are indexed by global face labels
    auto& faceQuadPoints = compactListListCRef
    (
        solMod.displacementLeastSquares().quadrature().faceQuadPoints()
    );

    labelList nQpPerFace(this->size(), 0);
    const label start = this->patch().patch().start();
    forAll(nQpPerFace, faceI)
    {
        nQpPerFace[faceI] = faceQuadPoints[faceI + start].size();
    }

    autoPtr<CompactListList<vector>> tQuadPointsValue
    (
        new CompactListList<vector>(nQpPerFace)
    );
    CompactListList<vector>& quadPointsValue = tQuadPointsValue();

    if (!mmsPtr_.valid())
    {
        mmsPtr_.reset(new neoHookeanManufacturedSolution(mesh));
    }

    const scalar t = this->db().time().value();

    forAll(*this, faceI)
    {
        const label globalFaceID = faceI + start;

        forAll(faceQuadPoints[globalFaceID], pointI)
        {
            quadPointsValue[faceI][pointI] =
                mmsPtr_->displacement(faceQuadPoints[globalFaceID][pointI], t);
        }
    }

    return tQuadPointsValue;
}


void neoHookeanManufacturedSolutionFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    neoHookeanManufacturedSolutionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
