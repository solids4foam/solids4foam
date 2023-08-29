/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "rigidCylinderContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rigidCylinderContactFvPatchVectorField::
rigidCylinderContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    radius_(0.0),
    cylinderCentre_(),
    penaltyStiffness_(0),
    relaxFactor_(0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


rigidCylinderContactFvPatchVectorField::
rigidCylinderContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    radius_(readScalar(dict.lookup("cylinderRadius"))),
    cylinderCentre_(dict.subDict("cylinderCentre")),
    penaltyStiffness_(readScalar(dict.lookup("penaltyStiffness"))),
    relaxFactor_(readScalar(dict.lookup("relaxFactor")))
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


rigidCylinderContactFvPatchVectorField::
rigidCylinderContactFvPatchVectorField
(
    const rigidCylinderContactFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(tdpvf, p, iF, mapper),
    radius_(tdpvf.radius_),
    cylinderCentre_(tdpvf.cylinderCentre_),
    penaltyStiffness_(tdpvf.penaltyStiffness_),
    relaxFactor_(tdpvf.relaxFactor_)
{}

#ifndef OPENFOAM_ORG
rigidCylinderContactFvPatchVectorField::
rigidCylinderContactFvPatchVectorField
(
    const rigidCylinderContactFvPatchVectorField& tdpvf
)
:
    solidTractionFvPatchVectorField(tdpvf),
    radius_(tdpvf.radius_),
    cylinderCentre_(tdpvf.cylinderCentre_),
    penaltyStiffness_(tdpvf.penaltyStiffness_),
    relaxFactor_(tdpvf.relaxFactor_)
{}
#endif

rigidCylinderContactFvPatchVectorField::
rigidCylinderContactFvPatchVectorField
(
    const rigidCylinderContactFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(tdpvf, iF),
    radius_(tdpvf.radius_),
    cylinderCentre_(tdpvf.cylinderCentre_),
    penaltyStiffness_(tdpvf.penaltyStiffness_),
    relaxFactor_(tdpvf.relaxFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rigidCylinderContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void rigidCylinderContactFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);

//     const rigidCylinderContactFvPatchVectorField& dmptf =
//         refCast<const rigidCylinderContactFvPatchVectorField>(ptf);
}


// Update the coefficients associated with the patch field
void rigidCylinderContactFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Let us first calculate the deformed position of the patch face centres
    // Note: patch().Cf() are the initial undeformed face centre positions and
    // *this are the current displacements
    const vectorField deformedCf(patch().Cf() + *this);

    // Calculate the cylinder centre position fo rthe current time
    const vector curCylinderCentre =
        cylinderCentre_(this->db().time().timeOutputValue());

    // Store radius square for efficiency
    const scalar radiusSqr = sqr(radius_);

    // Take a reference to the traction field
    vectorField& traction_ = traction();

    // Take a copy of previous traction field so we can apply under-relaxation
    // later
    const vectorField prevTraction(traction_);

    // Loop through all the faces in the patch
    forAll(deformedCf, faceI)
    {
        const vector& curFaceCf = deformedCf[faceI];

        // Check if the face is inside the cylinder
        // Note: we will assume that the cylinder axis points in the z direction
        if
        (
            sqr(curFaceCf.x() - curCylinderCentre.x())
          + sqr(curFaceCf.y() - curCylinderCentre.y())
          < radiusSqr
        )
        {
            // This face is in contact

            // Calculate the radial vector from the cylinder centre to face
            // position
            vector radialVec = curFaceCf - curCylinderCentre;

            // Calculate the overlap distance (the difference between the radius
            // and the length or the radialVec)
            scalar overlap = radius_ - mag(radialVec);

            // Calculate radial direction (direction the force will be applied)
            vector radialDir = radialVec/mag(radialVec);

            // Calculate the traction using the penalty method (equation of a
            // spring)
            traction_[faceI] = penaltyStiffness_*overlap*radialDir;
        }
        else
        {
            // The face is not in contact so we will set the traction to zero
            traction_[faceI] = vector::zero;
        }
    }

    // Apply under-relaxation
    traction_ = relaxFactor_*traction_ + (1.0 - relaxFactor_)*prevTraction;

    // Apply the traction
    solidTractionFvPatchVectorField::updateCoeffs();
}


// Write
void rigidCylinderContactFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
    os.writeKeyword("cylinderRadius")
        << radius_ << token::END_STATEMENT << nl;
    os.writeKeyword("penaltyStiffness")
        << penaltyStiffness_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxFactor")
        << relaxFactor_ << token::END_STATEMENT << nl;
    os.writeKeyword("cylinderCentre") << nl;
    os << token::BEGIN_BLOCK << nl;
    cylinderCentre_.write(os);
    os << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, rigidCylinderContactFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
