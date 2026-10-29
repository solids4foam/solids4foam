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

#include "linearSpringForceBeamDisplacementNRFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupBeamModel.H"
#include "pseudoVector.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    springStiffness_(0.0),
    springDirection_(vector::zero),
    Ks_(p.size(), tensor::zero),
    offsetForce_(p.size(), vector::zero),
    offsetForceSeries_()
{
    // fixedValueFvPatchVectorField::operator==(patchInternalField());
    // gradient() = vector::zero;
}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const linearSpringForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, p, iF, mapper),
    springStiffness_(tdpvf.springStiffness_),
    springDirection_(tdpvf.springDirection_),
    Ks_(tdpvf.Ks_, mapper),
    offsetForce_(tdpvf.offsetForce_, mapper),
    offsetForceSeries_(tdpvf.offsetForceSeries_)
{}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    springStiffness_(0.0),
    springDirection_(vector::zero),
    Ks_(p.size(), tensor::zero),
    offsetForce_(p.size(), vector::zero),
    offsetForceSeries_()
{
    // // NOT SURE ABOUT THIS
    // if (dict.found("value"))
    // {
    //     fixedValueFvPatchVectorField::operator==
    //     (
    //         vectorField("value", dict, p.size())
    //     );
    // }
    // else
    // {
        // fixedValueFvPatchVectorField::operator==
        // (
        //     patchInternalField()
        // );
    // }

    // Constructing the spring dictionary
    const dictionary& springDict = dict.subDict("linearSpringCoeffs");

    // Read the spring stiffness
    springStiffness_ = readScalar(springDict.lookup("springStiffness"));

    // Read the spring direction
    springDirection_ = vector(springDict.lookup("springDirection"));

    // Normalise the spring direction to unit vector
    if (mag(springDirection_) > SMALL)
    {
        springDirection_ /= mag(springDirection_);
    }

    // Spring stiffness tensor
    Ks_ = springStiffness_*(springDirection_*springDirection_);

    if (dict.found("offsetForceSeries"))
    {
        Info<< "Offset force found at spring BC of beam solver"
            << " and the force is time-varying" << endl;
        offsetForceSeries_ =
            interpolationTable<vector>(dict.subDict("offsetForceSeries"));

        offsetForce_ = offsetForceSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        Info<< "Offset force not found! skipping" << endl;
    }
}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const linearSpringForceBeamDisplacementNRFvPatchVectorField& tdpvf
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf),
    springStiffness_(tdpvf.springStiffness_),
    springDirection_(tdpvf.springDirection_),
    Ks_(tdpvf.Ks_),
    offsetForce_(tdpvf.offsetForce_),
    offsetForceSeries_(tdpvf.offsetForceSeries_)
{}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const linearSpringForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, iF),
    springStiffness_(tdpvf.springStiffness_),
    springDirection_(tdpvf.springDirection_),
    Ks_(tdpvf.Ks_),
    offsetForce_(tdpvf.offsetForce_),
    offsetForceSeries_(tdpvf.offsetForceSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void linearSpringForceBeamDisplacementNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    forceBeamDisplacementNRFvPatchVectorField::autoMap(m);
    Ks_.autoMap(m);
    offsetForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void linearSpringForceBeamDisplacementNRFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    forceBeamDisplacementNRFvPatchVectorField::rmap(ptf, addr);

    const linearSpringForceBeamDisplacementNRFvPatchVectorField& dmptf =
        refCast<const linearSpringForceBeamDisplacementNRFvPatchVectorField>(ptf);

    Ks_.rmap(dmptf.Ks_, addr);
    offsetForce_.rmap(dmptf.offsetForce_, addr);
}


// Update the coefficients associated with the patch field
void linearSpringForceBeamDisplacementNRFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (offsetForceSeries_.size())
    {
        offsetForce_ = offsetForceSeries_(this->db().time().timeOutputValue());
    }

    forceBeamDisplacementNRFvPatchVectorField::updateCoeffs();
}

void linearSpringForceBeamDisplacementNRFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    // Value of displacement at the boundary
    const vectorField& pW = *this;

    // Current Spring Force calculated from previous Newton iteration
    forAll(pW, pFace)
    {
        // Spring force magnitude calculated from displacement
        // component in the spring direction
        scalar curSpringForceMag =
            springStiffness_*(pW[pFace] & springDirection_);

        // Spring force vector directly updated into base
        // forceBeamDisplacementNR class.
        // This force() is also the used in the source vector
        // of global AX = B system.

        // This force is compressive and acts opposite of spring direction
        force()[pFace] =
            -curSpringForceMag*springDirection_
          + offsetForce_[pFace];
    }

    // Other necessary fields to enforce the boundary condition
    fvPatchField<vector>& DW =
        const_cast<fvPatchField<vector>&>
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

    // This is the old coefficient in force boundary condition
    // whose inverse was being multiplied
    const tensorField A(CQW/delta);

    // This coefficent changes due to the spring
    // The sign is +ve infront of Ks because the spring force
    // acting on the beam is negative (opp) to that acting on
    // the spring itself
    const tensorField invAmodified
    (
        inv
        (
            CQW/delta + Ks_
        )
    );

    // Modified displacement at boundary due to spring
    const vectorField newDW
    (
     (invAmodified & (force() - explicitQ))
      - (invAmodified & (CQTheta & DTheta))
      + (invAmodified & (A & DW.patchInternalField()))
    );

    DW = newDW;

    // Setting W_ field
    fixedValueFvPatchField<vector>::operator==((*this) + newDW);

    fixedValueFvPatchVectorField::evaluate();
}


// Write
void linearSpringForceBeamDisplacementNRFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    os.writeKeyword("linearSpringCoeffs") << nl;
    os << token::BEGIN_BLOCK << nl;
    os.writeKeyword("springStiffness") << springStiffness_ << token::END_STATEMENT << nl;
    os.writeKeyword("springDirection") << springDirection_ << token::END_STATEMENT << nl;
    os.writeKeyword("springForce") << force() << token::END_STATEMENT << nl;
    os << token::END_BLOCK << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    linearSpringForceBeamDisplacementNRFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
