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

#include "paraboloidInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

paraboloidInletVelocityFvPatchVectorField::
paraboloidInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(0.0),
    yMax_(0.0),
    zMax_(0.0),
    timeVaryingEndTime_(0.0),
    timeAtMaxVelocity_(0.0)
{}


paraboloidInletVelocityFvPatchVectorField::
paraboloidInletVelocityFvPatchVectorField
(
    const paraboloidInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umax_(ptf.Umax_),
    yMax_(ptf.yMax_),
    zMax_(ptf.zMax_),
    timeVaryingEndTime_(ptf.timeVaryingEndTime_),
    timeAtMaxVelocity_(ptf.timeAtMaxVelocity_)
{}


paraboloidInletVelocityFvPatchVectorField::
paraboloidInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    Umax_(readScalar(dict.lookup("maxVelocity"))),
    yMax_(readScalar(dict.lookup("yWidth"))),
    zMax_(readScalar(dict.lookup("zWidth"))),
    timeVaryingEndTime_(readScalar(dict.lookup("timeVaryingEndTime"))),
    timeAtMaxVelocity_
    (
        bool(timeVaryingEndTime_ > 0.0)
      ? readScalar(dict.lookup("timeAtMaxVelocity"))
      : 0.0
    )
{}


#ifndef OPENFOAM_ORG
paraboloidInletVelocityFvPatchVectorField::
paraboloidInletVelocityFvPatchVectorField
(
    const paraboloidInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    Umax_(pivpvf.Umax_),
    yMax_(pivpvf.yMax_),
    zMax_(pivpvf.zMax_),
    timeVaryingEndTime_(pivpvf.timeVaryingEndTime_),
    timeAtMaxVelocity_(pivpvf.timeAtMaxVelocity_)
{}
#endif


paraboloidInletVelocityFvPatchVectorField::
paraboloidInletVelocityFvPatchVectorField
(
    const paraboloidInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    Umax_(pivpvf.Umax_),
    yMax_(pivpvf.yMax_),
    zMax_(pivpvf.zMax_),
    timeVaryingEndTime_(pivpvf.timeVaryingEndTime_),
    timeAtMaxVelocity_(pivpvf.timeAtMaxVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void paraboloidInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField::subField Cf = patch().patch().faceCentres();

    // Calculate the max velocity for the current time
    scalar curMaxU = Umax_;
    if (db().time().value() < timeVaryingEndTime_)
    {
        curMaxU =
            Umax_*(1.0 - cos(M_PI*db().time().value()/timeAtMaxVelocity_))/2.0;
    }

    // Calculate the patch velocity
    vectorField patchU(patch().size(), vector::zero);
    forAll(patchU, faceI)
    {
        const scalar y = Cf[faceI].y();
        const scalar z = Cf[faceI].z();

        patchU[faceI] =
            //curMaxU*y*(zMax_ - y)*(sqr(zMax_) - sqr(z))
            curMaxU*y*(2*yMax_ - y)*(sqr(zMax_) - sqr(z)) // fix: A. Shay Jun-19
           *vector(1, 0, 0)/(sqr(yMax_)*sqr(zMax_));
    }

    fvPatchField<vector>::operator==(patchU);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> >
paraboloidInletVelocityFvPatchVectorField::snGrad() const
{
    bool secondOrder_ = false;

#ifdef OPENFOAM_NOT_EXTEND
    const word& UName = internalField().name();
#else
    const word& UName = dimensionedInternalField().name();
#endif

    if (db().foundObject<volTensorField>("grad(" + UName + ")"))
    {
        tmp<Field<vector> > tnGradU
        (
            new vectorField(this->patch().size(), vector::zero)
        );

        const fvPatchField<tensor>& gradU =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + UName + ")"
            );

        const vectorField n(patch().nf());
        const vectorField delta(patch().delta());
        const vectorField k((I - sqr(n)) & delta);

        if (secondOrder_)
        {
            const vectorField dUP(k & gradU.patchInternalField());
            const vectorField nGradUP(n & gradU.patchInternalField());

            #ifdef OPENFOAM_NOT_EXTEND
            tnGradU.ref() =
            #else
            tnGradU() =
            #endif
                2
               *(
                    *this
                  - (patchInternalField() + dUP)
                )*this->patch().deltaCoeffs()
              - nGradUP;

            return tnGradU;
        }

        // First order
        const vectorField dUP(k & gradU.patchInternalField());

        #ifdef OPENFOAM_NOT_EXTEND
        tnGradU.ref() =
        #else
        tnGradU() =
        #endif
        (
            *this
          - (patchInternalField() + dUP)
        )*patch().deltaCoeffs();

        return tnGradU;
    }

     return
     (
         *this - patchInternalField()
     )*patch().deltaCoeffs();
}

tmp<Field<vector> > paraboloidInletVelocityFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    bool secondOrder_ = false;

#ifdef OPENFOAM_NOT_EXTEND
    const word& UName = internalField().name();
#else
    const word& UName = dimensionedInternalField().name();
#endif

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    const vectorField n(patch().nf());
    const vectorField delta(patch().delta());
    const vectorField k((I - sqr(n)) & delta);

    if (secondOrder_)
    {
        const vectorField dUP(k & gradU.patchInternalField());
        const vectorField nGradUP(n & gradU.patchInternalField());

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dUP)
              - patchInternalField()
            )
          - nGradUP;
    }

    // First order
    const vectorField dUP(k & gradU.patchInternalField());

    return
        this->patch().deltaCoeffs()
       *(
           *this - dUP
        );
}


void paraboloidInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("maxVelocity")
        << Umax_ << token::END_STATEMENT << nl;

    os.writeKeyword("yWidth")
        << yMax_ << token::END_STATEMENT << nl;

    os.writeKeyword("zWidth")
        << zMax_ << token::END_STATEMENT << nl;

    os.writeKeyword("timeVaryingEndTime")
        << timeVaryingEndTime_ << token::END_STATEMENT << nl;

    os.writeKeyword("timeAtMaxVelocity")
        << timeAtMaxVelocity_ << token::END_STATEMENT << nl;

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    paraboloidInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
