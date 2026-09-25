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

#include "manufacturedSolutionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

manufacturedSolutionFvPatchVectorField::manufacturedSolutionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    mmsPtr_(),
    dict_()
{}


manufacturedSolutionFvPatchVectorField::manufacturedSolutionFvPatchVectorField
(
    const manufacturedSolutionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    mmsPtr_(),
    dict_(ptf.dict_)
{}


manufacturedSolutionFvPatchVectorField::manufacturedSolutionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    mmsPtr_(),
    dict_(dict)
{}

#ifndef OPENFOAM_ORG
manufacturedSolutionFvPatchVectorField::manufacturedSolutionFvPatchVectorField
(
    const manufacturedSolutionFvPatchVectorField& pivpvf
)
:
    fixedDisplacementFvPatchVectorField(pivpvf),
    mmsPtr_(),
    dict_(pivpvf.dict_)
{}
#endif

manufacturedSolutionFvPatchVectorField::manufacturedSolutionFvPatchVectorField
(
    const manufacturedSolutionFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(pivpvf, iF),
    mmsPtr_(),
    dict_(pivpvf.dict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void manufacturedSolutionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Set the MMS object, if required
    if (!mmsPtr_.valid())
    {
        mmsPtr_.reset
        (
            new manufacturedSolution(patch().boundaryMesh().mesh(), dict_)
        );
    }

    // Set the displacement at each patch face
    vectorField& disp = totalDisp();
    const vectorField& Cf = patch().Cf();
    forAll(disp, faceI)
    {
        disp[faceI] = mmsPtr_->calculateDisplacement(Cf[faceI]);
    }

    fixedDisplacementFvPatchVectorField::updateCoeffs();
}


autoPtr<CompactListList<vector>>
manufacturedSolutionFvPatchVectorField::evaluateQuadrature
() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const solidModel& solMod = lookupSolidModel(mesh);

    // Quadrature points are indexed by global face labels
    const CompactListList<point>& faceQuadPoints =
        solMod.displacementLeastSquares().quadrature().faceQuadPoints();

    labelList nQpPerFace(this->size(), 0);
    const label start = this->patch().start();
    forAll(nQpPerFace, faceI)
    {
        const label globalFaceID = faceI + start;
        nQpPerFace[faceI] = faceQuadPoints[globalFaceID].size();
    }

    autoPtr<CompactListList<vector>> tQuadPointsValue
    (
        new CompactListList<vector>(nQpPerFace)
    );

    // Get a reference to the actual data for easier access
    CompactListList<vector>& quadPointsValue = tQuadPointsValue();

    // Set the MMS object, if required
    if (!mmsPtr_.valid())
    {
        mmsPtr_.reset
        (
            new manufacturedSolution(patch().boundaryMesh().mesh(), dict_)
        );
    }

    // Loop over faces
    forAll(*this, faceI)
    {
        const label globalFaceID = faceI + start;

        // Get the number of quadrature points for this face
        const label nPoints = faceQuadPoints[globalFaceID].size();

        // Assign the values to face quadrature points
        for (label pointI = 0; pointI < nPoints; ++pointI)
        {
            const point quadPoint = faceQuadPoints[globalFaceID][pointI];
            quadPointsValue[faceI][pointI] =
                mmsPtr_->calculateDisplacement(quadPoint);
        }
    }

    return tQuadPointsValue;
}


void manufacturedSolutionFvPatchVectorField::write(Ostream& os) const
{
    // os.writeKeyword("rate")
    //     << rate_ << token::END_STATEMENT << nl;
    // initialC_.writeEntry("initialC", os);

    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    manufacturedSolutionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
