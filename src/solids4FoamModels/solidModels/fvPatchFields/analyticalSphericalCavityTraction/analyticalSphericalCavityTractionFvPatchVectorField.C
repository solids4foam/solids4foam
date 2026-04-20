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

#include "analyticalSphericalCavityTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "sphericalCavityStressDisplacement.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalSphericalCavityTractionFvPatchVectorField::
analyticalSphericalCavityTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    T0_(0.0),
    cavityR_(0.0),
    nu_(0.0)
{}


analyticalSphericalCavityTractionFvPatchVectorField::
analyticalSphericalCavityTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    T0_(readScalar(dict.lookup("farFieldTractionZ"))),
    cavityR_(readScalar(dict.lookup("cavityRadius"))),
    nu_(readScalar(dict.lookup("nu")))
{}


analyticalSphericalCavityTractionFvPatchVectorField::
analyticalSphericalCavityTractionFvPatchVectorField
(
    const analyticalSphericalCavityTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
    T0_(stpvf.T0_),
    cavityR_(stpvf.cavityR_),
    nu_(stpvf.nu_)
{}

#ifndef OPENFOAM_ORG
analyticalSphericalCavityTractionFvPatchVectorField::
analyticalSphericalCavityTractionFvPatchVectorField
(
    const analyticalSphericalCavityTractionFvPatchVectorField& stpvf
)
:
    solidTractionFvPatchVectorField(stpvf),
    T0_(stpvf.T0_),
    cavityR_(stpvf.cavityR_),
    nu_(stpvf.nu_)
{}
#endif

analyticalSphericalCavityTractionFvPatchVectorField::
analyticalSphericalCavityTractionFvPatchVectorField
(
    const analyticalSphericalCavityTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF),
    T0_(stpvf.T0_),
    cavityR_(stpvf.cavityR_),
    nu_(stpvf.nu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalSphericalCavityTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void analyticalSphericalCavityTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void analyticalSphericalCavityTractionFvPatchVectorField::updateCoeffs()
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
        trac[faceI] =
            (
                n[faceI] & sphericalCavityStress(T0_, nu_, cavityR_, Cf[faceI])
            );
       }

    solidTractionFvPatchVectorField::updateCoeffs();
}


#ifndef FOAMEXTEND
autoPtr<CompactListList<vector>>
analyticalSphericalCavityTractionFvPatchVectorField::evaluateQuadrature() const
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
                n[faceI] & sphericalCavityStress(T0_, nu_, cavityR_, quadPoint);
        }
    }

    return tQuadPointsValue;
}
#endif


// Write
void analyticalSphericalCavityTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("farFieldTractionZ")
        << T0_ << token::END_STATEMENT << nl;

    os.writeKeyword("cavityRadius")
        << cavityR_ << token::END_STATEMENT << nl;

    os.writeKeyword("nu")
        << nu_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    analyticalSphericalCavityTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
