/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "movingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero)
{}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    myTimeIndex_(ptf.myTimeIndex_),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero)
{}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    const fvMesh& mesh = dimensionedInternalField().mesh();
    const pointField& points = mesh.allPoints();

    forAll(Fc_, i)
    {
        Fc_[i] = patch().patch()[i].centre(points);
    }

    oldFc_ = Fc_;
    oldoldFc_ = Fc_;
}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.Fc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_)
{}


movingWallVelocityFvPatchVectorField::movingWallVelocityFvPatchVectorField
(
    const movingWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.oldFc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingWallVelocityFvPatchVectorField::updateCoeffs()
{
//     Info << "movingWallVelocityFvPatchVectorField::updateCoeffs" << endl;

    if (updated())
    {
        return;
    }

    const fvMesh& mesh = dimensionedInternalField().mesh();
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();
    const pointField& points = mesh.allPoints();

    vectorField Up(p.size(), vector::zero);

    //const pointField& oldPoints = mesh.oldPoints();
    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    word ddtScheme
    (
        mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
    );

    if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
    {
        if (myTimeIndex_ < mesh.time().timeIndex())
        {
            oldoldFc_ = oldFc_;
            oldFc_ = Fc_;
            // Fc_ = pp.faceCentres();

            myTimeIndex_ = mesh.time().timeIndex();
        }

        forAll(Fc_, i)
        {
            Fc_[i] = pp[i].centre(points);
        }

        scalar deltaT = mesh.time().deltaT().value();
        scalar deltaT0 = mesh.time().deltaT0().value();
        if
        (
            U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex()
         || U.oldTime().oldTime().timeIndex() < 0
        )
        {
            deltaT0 = GREAT;
        }

        //Set coefficients based on deltaT and deltaT0
        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
//         scalar coefft0  = coefft + coefft00;

//         Up = (coefft*Fc_ - coefft0*oldFc_ + coefft00*oldoldFc_)
//            /mesh.time().deltaT().value();

        Up = coefft*(Fc_ - oldFc_)/deltaT
          - coefft00*(oldFc_ - oldoldFc_)/deltaT;

//         Info << max(mag(Up)) << endl;
    }
    else // Euler
    {
//         Info << "void movingWallVelocityFvPatchVectorField::updateCoeffs() - "
//             << "Euler"
//             << endl;

        if (myTimeIndex_ < mesh.time().timeIndex())
        {
            oldoldFc_ = oldFc_;
            oldFc_ = Fc_;

//             Fc_ = pp.faceCentres();
            myTimeIndex_ = mesh.time().timeIndex();
        }

        forAll(Fc_, i)
        {
            Fc_[i] = pp[i].centre(points);
        }

        Up = (Fc_ - oldFc_)/mesh.time().deltaT().value();
    }

    scalarField phip =
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

    vectorField n = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    Up += n*(Un - (n & Up));

    vectorField::operator=(Up);

//     Info << "mwvuc " << max(mag(Up)) << ", " << average(mag(Up)) << endl;

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> > movingWallVelocityFvPatchVectorField::
snGrad() const
{
    bool secondOrder_ = false;

    word UName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField dUP = (k&gradU.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradUP = (n&gradU.patchInternalField());

        tmp<Field<vector> > tnGradU
        (
            new vectorField(this->patch().size(), vector::zero)
        );

        tnGradU() =
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs()
          - nGradUP;

        tnGradU() -= n*(n&tnGradU());

        return tnGradU;

//         return
//             2
//            *(
//                 *this
//               - (patchInternalField() + dUP)
//             )*this->patch().deltaCoeffs()
//           - nGradUP;
    }


    // First order
//     vectorField dUP = (k&gradU.patchInternalField());

    tmp<Field<vector> > tnGradU
    (
        new vectorField(this->patch().size(), vector::zero)
    );

    tnGradU() =
        (
            *this
          - (patchInternalField() + dUP)
        )*this->patch().deltaCoeffs();

    tnGradU() -= n*(n&tnGradU());

    return tnGradU;

//     return
//     (
//         *this
//       - (patchInternalField() + (k&gradU.patchInternalField()))
//     )*this->patch().deltaCoeffs();
}


tmp<Field<vector> > movingWallVelocityFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    bool secondOrder_ = false;

    word UName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField dUP = (k&gradU.patchInternalField());

    if (secondOrder_)
    {
        vectorField nGradUP = (n&gradU.patchInternalField());

        vectorField nGradU =
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs()
          - nGradUP;

        vectorField nGradUn = n*(n&nGradU);

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dUP)
              - patchInternalField()
            )
          - nGradUP
          - nGradUn;

//         return
//             this->patch().deltaCoeffs()
//            *(
//                 2*(*this - dUP)
//               - patchInternalField()
//             )
//           - nGradUP;
    }


    // First order
//     vectorField dUP = (k&gradU.patchInternalField());

    vectorField nGradU =
        (
            *this
          - (patchInternalField() + dUP)
        )*this->patch().deltaCoeffs();

    vectorField nGradUn = n*(n&nGradU);

    return
        this->patch().deltaCoeffs()
       *(
           *this - dUP
        )
      - nGradUn;

//     return this->patch().deltaCoeffs()
//        *(*this - (k&gradU.patchInternalField()));
}



void movingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingWallVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
