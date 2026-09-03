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

InClass
    cantileverTractionFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "cantileverTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "GeometricFields.H"
#include "lookupSolidModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cantileverTractionFvPatchVectorField::cantileverTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    analyticalSol_(),
    curTimeIndex_(-1)
{}


Foam::cantileverTractionFvPatchVectorField::cantileverTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    analyticalSol_(dict),
    curTimeIndex_(-1)
{
    Info<< "Creating " << cantileverTractionFvPatchVectorField::typeName
        << " patch" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }
}


Foam::cantileverTractionFvPatchVectorField::cantileverTractionFvPatchVectorField
(
    const cantileverTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}

#ifndef OPENFOAM_ORG
Foam::cantileverTractionFvPatchVectorField::cantileverTractionFvPatchVectorField
(
    const cantileverTractionFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif


Foam::cantileverTractionFvPatchVectorField::cantileverTractionFvPatchVectorField
(
    const cantileverTractionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::cantileverTractionFvPatchVectorField::
~cantileverTractionFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cantileverTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


void Foam::cantileverTractionFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


void Foam::cantileverTractionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = db().time().timeIndex();

        // Patch face coordinates
        const vectorField& f = patch().Cf();

        // Set the face traction field
        traction() = vectorField(analyticalSol_.traction(f));
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


#ifndef FOAMEXTEND
Foam::autoPtr<Foam::CompactListList<Foam::vector>>
Foam::cantileverTractionFvPatchVectorField::evaluateQuadrature
() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const solidModel& solMod = lookupSolidModel(mesh);

    // faceQuadPoints is list for the  whole mesh
    const CompactListList<point>& faceQuadPoints =
        solMod.displacementLeastSquares().quadrature().faceQuadPoints();

    // faceQuadPoints is list for whole mesh.
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
                analyticalSol_.traction(quadPoint);
        }
    }

    return tQuadPointsValue;
}
#endif


void Foam::cantileverTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    analyticalSol_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cantileverTractionFvPatchVectorField
    );
}

// ************************************************************************* //
