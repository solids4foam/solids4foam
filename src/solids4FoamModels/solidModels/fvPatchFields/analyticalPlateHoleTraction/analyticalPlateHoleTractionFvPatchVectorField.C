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

#include "analyticalPlateHoleTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mechanicalModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "lookupSolidModel.H"
#include "plateHoleAnalyticalFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

symmTensor analyticalPlateHoleTractionFvPatchVectorField::plateHoleSolution
(
    const vector& C
) const
{
    return plateHoleAnalyticalFields::stress(C, T_, holeR_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    T_(0.0),
    holeR_(0.0)
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius")))
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}

#ifndef OPENFOAM_ORG
analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf
)
:
    solidTractionFvPatchVectorField(stpvf),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}
#endif

analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF),
    T_(stpvf.T_),
    holeR_(stpvf.holeR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalPlateHoleTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void analyticalPlateHoleTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void analyticalPlateHoleTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Patch unit normals
    vectorField n(patch().nf());

    // Patch face centres
    const vectorField& Cf = patch().Cf();

    // Set the patch traction

    vectorField& trac = traction();

    forAll(traction(), faceI)
    {
        vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
        vector curN = n[faceI];

        if (patch().name() == "hole")
        {
            curC /= mag(curC);
            curC *= holeR_;

            curN = -curC/mag(curC);
        }

        trac[faceI] =
            plateHoleAnalyticalFields::traction(curC, curN, T_, holeR_);
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}

#ifndef FOAMEXTEND
autoPtr<CompactListList<vector>>
analyticalPlateHoleTractionFvPatchVectorField::evaluateQuadrature() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const solidModel& solMod = lookupSolidModel(mesh);

    // faceQuadPoints is list for the  whole mesh
    const CompactListList<point>& faceQuadPoints =
        solMod.displacementMLS().quadrature().faceQuadPoints();

    labelList nQpPerFace(this->size(), 0);
    const label start = this->patch().start();
    forAll(nQpPerFace, faceI)
    {
        const label globalFaceID = faceI + start;
        nQpPerFace[faceI]=faceQuadPoints[globalFaceID].size();
    }

    autoPtr<CompactListList<vector>> tQuadPointsValue
    (
        new CompactListList<vector>(nQpPerFace)
    );

    // Get a reference to the actual data for easier access
    CompactListList<vector>& quadPointsValue = tQuadPointsValue();

    // Patch unit normals
    vectorField n(patch().nf());

    forAll(*this, faceI)
    {
        const label globalFaceID = faceI + start;

        // Get the number of quadrature points for this face
        const label nPoints = faceQuadPoints[globalFaceID].size();

        // Assign the same value to all quadrature points on this face
        // We assume constant distribution of traction!
        for (label pointI = 0; pointI < nPoints; ++pointI)
        {
            const point quadPoint = faceQuadPoints[globalFaceID][pointI];
            quadPointsValue[faceI][pointI] =
                plateHoleAnalyticalFields::traction
                (
                    quadPoint,
                    n[faceI],
                    T_,
                    holeR_
                );
        }
    }

    return tQuadPointsValue;
}
#endif


// Write
void analyticalPlateHoleTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("farFieldTractionX")
        << T_ << token::END_STATEMENT << nl;

    os.writeKeyword("holeRadius")
        << holeR_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    analyticalPlateHoleTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
