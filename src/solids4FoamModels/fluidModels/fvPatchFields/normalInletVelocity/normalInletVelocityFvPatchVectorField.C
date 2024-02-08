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

#include "normalInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

normalInletVelocityFvPatchVectorField::normalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


normalInletVelocityFvPatchVectorField::normalInletVelocityFvPatchVectorField
(
    const normalInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


normalInletVelocityFvPatchVectorField::normalInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{
//     fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

// #ifdef OPENFOAM_NOT_EXTEND
//     const fvMesh& mesh = internalField().mesh();
//     const pointField& points = mesh.points();
// #else
//     const fvMesh& mesh = dimensionedInternalField().mesh();
//     const pointField& points = mesh.allPoints();
// #endif

//     forAll(Fc_, i)
//     {
//         Fc_[i] = patch().patch()[i].centre(points);
//     }

//     oldFc_ = Fc_;
//     oldoldFc_ = Fc_;
}


#ifndef OPENFOAM_ORG
normalInletVelocityFvPatchVectorField::normalInletVelocityFvPatchVectorField
(
    const normalInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf)
{}
#endif


normalInletVelocityFvPatchVectorField::normalInletVelocityFvPatchVectorField
(
    const normalInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void normalInletVelocityFvPatchVectorField::updateCoeffs()
// {
//     if (updated())
//     {
//         return;
//     }

// #ifdef OPENFOAM_NOT_EXTEND
//     const fvMesh& mesh = internalField().mesh();
//     const pointField& points = mesh.points();
// #else
//     const fvMesh& mesh = dimensionedInternalField().mesh();
//     const pointField& points = mesh.allPoints();
// #endif
//     const fvPatch& p = patch();
//     const polyPatch& pp = p.patch();

//     vectorField Up(p.size(), vector::zero);

//     //const pointField& oldPoints = mesh.oldPoints();
//     const volVectorField& U =
//         mesh.lookupObject<volVectorField>
//         (
//             dimensionedInternalField().name()
//         );

//     word ddtScheme
//     (
//         mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
//     );

//     if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
//     {
//         Info << "void normalInletVelocityFvPatchVectorField::updateCoeffs()"
//             << endl;

//         if(myTimeIndex_ < mesh.time().timeIndex())
//         {
//             oldoldFc_ = oldFc_;
//             oldFc_ = Fc_;
// //             Fc_ = pp.faceCentres();

//             forAll(Fc_, i)
//             {
//                 Fc_[i] = pp[i].centre(points);
//             }

//             myTimeIndex_ = mesh.time().timeIndex();
//         }

//         scalar deltaT = mesh.time().deltaT().value();
//         scalar deltaT0 = mesh.time().deltaT0().value();
//         if
//         (
//             U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex()
//          || U.oldTime().oldTime().timeIndex() < 0
//         )
//         {
//             deltaT0 = GREAT;
//         }

//         //Set coefficients based on deltaT and deltaT0
//         scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
//         scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
// //         scalar coefft0  = coefft + coefft00;

// //         Up = (coefft*Fc_ - coefft0*oldFc_ + coefft00*oldoldFc_)
// //            /mesh.time().deltaT().value();


//         Up = coefft*(Fc_ - oldFc_)/deltaT
//           + coefft00*(oldFc_ - oldoldFc_)/deltaT;

// //         Info << max(mag(Up)) << endl;
//     }
//     else // Euler
//     {
//         Info << "void normalInletVelocityFvPatchVectorField::updateCoeffs()"
//             << endl;

//         if(myTimeIndex_ < mesh.time().timeIndex())
//         {
//             oldFc_ = Fc_;

//             forAll(Fc_, i)
//             {
//                 Fc_[i] = pp[i].centre(points);
//             }

// //             Fc_ = pp.faceCentres();
//             myTimeIndex_ = mesh.time().timeIndex();
//         }

//         Up = (Fc_ - oldFc_)/mesh.time().deltaT().value();
//     }

//     scalarField phip =
//         p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

//     vectorField n = p.nf();
//     const scalarField& magSf = p.magSf();
//     scalarField Un = phip/(magSf + VSMALL);

// //     Info << max(mag(Un)) << endl;

//     Up += n*(Un - (n & Up));

// //     Info << max(mag(Up)) << endl;

//     vectorField::operator=(Up);

//     fixedValueFvPatchVectorField::updateCoeffs();
// }


Foam::tmp<Foam::Field<vector> > normalInletVelocityFvPatchVectorField::
snGrad() const
{
//     Info << "snGrad - in" << endl;

    bool secondOrder_ = false;

#ifdef OPENFOAM_NOT_EXTEND
    word UName = this->internalField().name();
#else
    word UName = this->dimensionedInternalField().name();
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

        tmp<Field<vector> > tnGradU
        (
            new vectorField(this->patch().size(), vector::zero)
        );

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

//         tnGradU() -= n*(n&tnGradU());

        return tnGradU;

//         return
//             2
//            *(
//                 *this
//               - (patchInternalField() + dUP)
//             )*this->patch().deltaCoeffs()
//           - nGradUP;
    }

//     Info << "snGrad - out" << endl;

//     return fixedValueFvPatchVectorField::snGrad();



    // First order
    const vectorField dUP(k & gradU.patchInternalField());

    tmp<Field<vector> > tnGradU
    (
        new vectorField(this->patch().size(), vector::zero)
    );

#ifdef OPENFOAM_NOT_EXTEND
    tnGradU.ref() =
#else
    tnGradU() =
#endif
        (
            *this
          - (patchInternalField() + dUP)
        )*this->patch().deltaCoeffs();

//     tnGradU() -= n*(n&tnGradU());

    return tnGradU;



//     return
//     (
//         *this
//       - (patchInternalField() + (k&gradU.patchInternalField()))
//     )*this->patch().deltaCoeffs();
}


tmp<Field<vector> > normalInletVelocityFvPatchVectorField::
gradientBoundaryCoeffs() const
{
//     Info << "gradientBoundaryCoeffs - in" << endl;

    bool secondOrder_ = false;

#ifdef OPENFOAM_NOT_EXTEND
    word UName = this->internalField().name();
#else
    word UName = this->dimensionedInternalField().name();
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

//         vectorField nGradU =
//             2
//            *(
//                 *this
//               - (patchInternalField() + dUP)
//             )*this->patch().deltaCoeffs()
//           - nGradUP;

//         vectorField nGradUn = n*(n&nGradU);

        return
            this->patch().deltaCoeffs()
           *(
                2*(*this - dUP)
              - patchInternalField()
            )
          - nGradUP;
//           - nGradUn;

//         return
//             this->patch().deltaCoeffs()
//            *(
//                 2*(*this - dUP)
//               - patchInternalField()
//             )
//           - nGradUP;
    }

//     tmp<Field<vector> > tResult
//     (
//         new vectorField(this->patch().size(), vector::zero)
//     );

//     tResult() = this->patch().deltaCoeffs()
//        *(*this - (k&gradU.patchInternalField()));

//     Info << "gradientBoundaryCoeffs - out" << endl;

//     return tResult;

//     return fixedValueFvPatchVectorField::gradientBoundaryCoeffs();



    // First order
    const vectorField dUP(k & gradU.patchInternalField());

//     vectorField nGradU =
//         (
//             *this
//           - (patchInternalField() + dUP)
//         )*this->patch().deltaCoeffs();

//     vectorField nGradUn = n*(n&nGradU);

    return
        this->patch().deltaCoeffs()
       *(
           *this - dUP
        );
//       - nGradUn;


//     return this->patch().deltaCoeffs()
//        *(*this - (k&gradU.patchInternalField()));
}



void normalInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
//     writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    normalInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
