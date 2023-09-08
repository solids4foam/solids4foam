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

#include "newMovingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

newMovingWallVelocityFvPatchVectorField::newMovingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
#ifdef OPENFOAM_NOT_EXTEND
    myTimeIndex_(internalField().mesh().time().timeIndex()),
#else
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldoldFc_(p.patch().size(), vector::zero),
    acceleration_(p.patch().size(), vector::zero)
{}


newMovingWallVelocityFvPatchVectorField::newMovingWallVelocityFvPatchVectorField
(
    const newMovingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    myTimeIndex_(ptf.myTimeIndex_),
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldoldFc_(p.patch().size(), vector::zero),
    acceleration_(p.patch().size(), vector::zero)
{}


newMovingWallVelocityFvPatchVectorField::newMovingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
#ifdef OPENFOAM_NOT_EXTEND
    myTimeIndex_(internalField().mesh().time().timeIndex()),
#else
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    Fc_(p.patch().size(), vector::zero),
    oldFc_(p.patch().size(), vector::zero),
    oldoldFc_(p.patch().size(), vector::zero),
    acceleration_(p.patch().size(), vector::zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
    const pointField& points = mesh.points();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const pointField& points = mesh.allPoints();
#endif

    forAll(Fc_, i)
    {
        Fc_[i] = patch().patch()[i].centre(points);
    }

    oldFc_ = Fc_;
    oldoldFc_ = Fc_;
}


#ifndef OPENFOAM_ORG
newMovingWallVelocityFvPatchVectorField::newMovingWallVelocityFvPatchVectorField
(
    const newMovingWallVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.Fc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_),
    acceleration_(pivpvf.acceleration_)
{}
#endif


newMovingWallVelocityFvPatchVectorField::newMovingWallVelocityFvPatchVectorField
(
    const newMovingWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.oldFc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_),
    acceleration_(pivpvf.acceleration_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void newMovingWallVelocityFvPatchVectorField::updateCoeffs()
{
//     Info << "newMovingWallVelocityFvPatchVectorField::updateCoeffs" << endl;

    if (updated())
    {
        return;
    }

#ifdef OPENFOAM_NOT_EXTEND
    const fvMesh& mesh = internalField().mesh();
    const pointField& points = mesh.points();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const pointField& points = mesh.allPoints();
#endif
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();

    vectorField Up(p.size(), vector::zero);

    //const pointField& oldPoints = mesh.oldPoints();
    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
#ifdef OPENFOAM_NOT_EXTEND
            internalField().name()
#else
            dimensionedInternalField().name()
#endif
        );

    word ddtScheme
    (
#ifdef OPENFOAM_NOT_EXTEND
        mesh.ddtScheme("ddt(" + U.name() +')')
#else
        mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
#endif
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
//         Info << "void newMovingWallVelocityFvPatchVectorField::updateCoeffs() - "
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

    const vectorField n(p.nf());
    const scalarField& magSf = p.magSf();
    const scalarField Un(phip/(magSf + VSMALL));

    Up += n*(Un - (n & Up));

    vectorField::operator=(Up);

//     Info << "mwvuc " << max(mag(Up)) << ", " << average(mag(Up)) << endl;

    // Update acceleration
    if
    (
        (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
     && (myTimeIndex_ > 1)
    )
    {
        scalar deltaT = mesh.time().deltaT().value();
        scalar deltaT0 = mesh.time().deltaT0().value();

        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        scalar coefft0  = coefft + coefft00;

        acceleration_ =
            (
                coefft*Up
              - coefft0
               *U.oldTime().boundaryField()[this->patch().index()]
              + coefft00
               *U.oldTime().oldTime().boundaryField()[this->patch().index()]
            )
           /mesh.time().deltaT().value();
    }
    else
    {
        acceleration_ =
            (Up - U.oldTime().boundaryField()[this->patch().index()])/
            mesh.time().deltaT().value();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> > newMovingWallVelocityFvPatchVectorField::
snGrad() const
{
    bool secondOrder_ = false;

#ifdef OPENFOAM_NOT_EXTEND
    word UName = internalField().name();
#else
    word UName = dimensionedInternalField().name();
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

        const vectorField dUP(k & gradU.patchInternalField());

        if (secondOrder_)
        {
            const vectorField nGradUP(n & gradU.patchInternalField());

            tmp<Field<vector> > tnGradU
            (
                new vectorField(this->patch().size(), vector::zero)
            );

    #ifdef OPENFOAM_NOT_EXTEND
            tnGradU.ref() =
                2
               *(
                    *this
                  - (patchInternalField() + dUP)
                )*this->patch().deltaCoeffs()
              - nGradUP;

            tnGradU.ref() -= n*(n&tnGradU());
    #else
            tnGradU() =
                2
               *(
                    *this
                  - (patchInternalField() + dUP)
                )*this->patch().deltaCoeffs()
              - nGradUP;

            tnGradU() -= n*(n&tnGradU());
    #endif

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
    // vectorField dUP = (k&gradU.patchInternalField());

    #ifdef OPENFOAM_NOT_EXTEND
        tnGradU.ref() =
            (
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs();

        tnGradU.ref() -= n*(n&tnGradU());
    #else
        tnGradU() =
            (
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs();

        tnGradU() -= n*(n&tnGradU());
    #endif

        return tnGradU;
    }

    return
    (
        (I - sqr(patch().nf()))
      & (
         *this - patchInternalField()
         )*this->patch().deltaCoeffs()
    );
}


tmp<Field<vector> > newMovingWallVelocityFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    bool secondOrder_ = false;

#ifdef OPENFOAM_NOT_EXTEND
    word UName = internalField().name();
#else
    word UName = dimensionedInternalField().name();
#endif

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    const vectorField n(patch().nf());
    const vectorField delta(patch().delta());
    const vectorField k((I - sqr(n)) & delta);

    const vectorField dUP(k & gradU.patchInternalField());

    if (secondOrder_)
    {
        const vectorField nGradUP(n & gradU.patchInternalField());

        const vectorField nGradU
        (
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*this->patch().deltaCoeffs()
          - nGradUP
        );

        const vectorField nGradUn(sqr(n) & nGradU);

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

    const vectorField nGradU
    (
        (
            *this
          - (patchInternalField() + dUP)
        )*this->patch().deltaCoeffs()
    );

    const vectorField nGradUn(sqr(n) & nGradU);

    return
        this->patch().deltaCoeffs()
       *(
           *this - dUP
        )
      - nGradUn;

//     return this->patch().deltaCoeffs()
//        *(*this - (k&gradU.patchInternalField()));
}



void newMovingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
#ifdef OPENFOAM_ORG
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    newMovingWallVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
