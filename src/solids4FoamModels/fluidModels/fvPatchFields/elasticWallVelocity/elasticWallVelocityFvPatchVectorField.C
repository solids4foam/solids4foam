/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "elasticWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

elasticWallVelocityFvPatchVectorField::elasticWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
#ifdef OPENFOAMESIORFOUNDATION
    timeIndex_(internalField().mesh().time().timeIndex()),
#else
    timeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
#endif
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldOldFc_(p.patch().size(),vector::zero)
{}


elasticWallVelocityFvPatchVectorField::elasticWallVelocityFvPatchVectorField
(
    const elasticWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    timeIndex_(ptf.timeIndex_),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldOldFc_(p.patch().size(),vector::zero)
{}


elasticWallVelocityFvPatchVectorField::elasticWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    timeIndex_(-1),
    //timeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldOldFc_(p.patch().size(),vector::zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

#ifdef OPENFOAMESIORFOUNDATION
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
    oldOldFc_ = Fc_;
}


#ifndef OPENFOAMFOUNDATION
elasticWallVelocityFvPatchVectorField::elasticWallVelocityFvPatchVectorField
(
    const elasticWallVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    timeIndex_(pivpvf.timeIndex_),
    Fc_(pivpvf.Fc_),
    oldFc_(pivpvf.oldFc_),
    oldOldFc_(pivpvf.oldOldFc_)
{}
#endif


elasticWallVelocityFvPatchVectorField::elasticWallVelocityFvPatchVectorField
(
    const elasticWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    timeIndex_(pivpvf.timeIndex_),
    Fc_(pivpvf.oldFc_),
    oldFc_(pivpvf.oldFc_),
    oldOldFc_(pivpvf.oldOldFc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void elasticWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

#ifdef OPENFOAMESIORFOUNDATION
    const fvMesh& mesh = internalField().mesh();
#else
    const fvMesh& mesh = dimensionedInternalField().mesh();
#endif

    const pointField& points = mesh.points();
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();

    vectorField Up(p.size(), vector::zero);

    //const pointField& oldPoints = mesh.oldPoints();
    const volVectorField& U =
        mesh.lookupObject<volVectorField>
        (
#ifdef OPENFOAMESIORFOUNDATION
            internalField().name()
#else
            dimensionedInternalField().name()
#endif
        );

    word ddtScheme
    (
#ifdef OPENFOAMESIORFOUNDATION
        mesh.ddtScheme("ddt(" + U.name() +')')
#else
        mesh.schemesDict().ddtScheme("ddt(" + U.name() +')')
#endif
    );

    // Compute velocity vector
    if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
    {
        if (timeIndex_ < mesh.time().timeIndex())
        {
            oldOldFc_ = oldFc_;
            oldFc_ = Fc_;
            //Fc_ = pp.faceCentres();

            timeIndex_ = mesh.time().timeIndex();
        }

        forAll(Fc_, i)
        {
            Fc_[i] = pp[i].centre(points);
        }

        const scalar deltaT = mesh.time().deltaT().value();
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
        const scalar coefft = 1 + deltaT/(deltaT + deltaT0);
        const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));

        Up = coefft*(Fc_ - oldFc_)/deltaT
             - coefft00*(oldFc_ - oldOldFc_)/deltaT;
    }
    else // Assume Euler
    {
        if (timeIndex_ < mesh.time().timeIndex())
        {
            oldOldFc_ = oldFc_;
            oldFc_ = Fc_;

            timeIndex_ = mesh.time().timeIndex();
        }

        forAll(Fc_, i)
        {
            Fc_[i] = pp[i].centre(points);
        }

        Up = (Fc_ - oldFc_)/mesh.time().deltaT().value();
    }

    // Compute normal velocity component (from mesh flux
    const scalarField phip
    (
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
    );

    const vectorField n(p.nf());
    const scalarField& magSf = p.magSf();
    scalarField Un(phip/(magSf + VSMALL));

    if (mesh.foundObject<surfaceScalarField>("phi"))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>("phi");

        const scalarField phipNew
        (
            phi.boundaryField()[this->patch().index()]
        );

        Un = phipNew/(magSf + VSMALL);
    }
    else
    {
        const volScalarField& pressure =
            mesh.lookupObject<volScalarField>("p");

        const scalarField nGradP
        (
            pressure.boundaryField()[patch().index()].snGrad()
        );

        const IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        const dimensionedScalar rho
        (
            transportProperties.lookup("rho")
        );

        const vectorField UOld(U.oldTime().boundaryField()[patch().index()]);

        Un = (n & UOld) - nGradP*mesh.time().deltaT().value();
        //Un = (n & UOld) - nGradP*mesh.time().deltaT().value()/rho.value();
    }

    // Compute tangential + normal components (from flux)
    Up += n*(Un - (n & Up));

    vectorField::operator=(Up);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> > elasticWallVelocityFvPatchVectorField::
snGrad() const
{
    bool secondOrder_ = false;

#ifdef OPENFOAMESIORFOUNDATION
    const word UName = internalField().name();
#else
    const word UName = dimensionedInternalField().name();
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

        tmp<Field<vector> > tnGradU
        (
            new vectorField(this->patch().size(), vector::zero)
        );

#ifdef OPENFOAMESIORFOUNDATION
        tnGradU.ref() =
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*patch().deltaCoeffs()
          - nGradUP;

        tnGradU.ref() -= (sqr(n) & tnGradU());
#else
        tnGradU() =
            2
           *(
                *this
              - (patchInternalField() + dUP)
            )*patch().deltaCoeffs()
          - nGradUP;

        tnGradU() -= (sqr(n) & tnGradU());
#endif

        return tnGradU;

        //return
        //    2
        //   *(
        //        *this
        //      - (patchInternalField() + dUP)
        //    )*this->patch().deltaCoeffs()
        //  - nGradUP;
    }


    // First order
    //vectorField dUP = (k&gradU.patchInternalField());

    tmp<Field<vector> > tnGradU
    (
        new vectorField(this->patch().size(), vector::zero)
    );

#ifdef OPENFOAMESIORFOUNDATION
    tnGradU.ref() =
        (
            *this
          - (patchInternalField() + dUP)
        )*patch().deltaCoeffs();

    tnGradU.ref() -= (sqr(n) & tnGradU());
#else
    tnGradU() =
        (
            *this
          - (patchInternalField() + dUP)
        )*patch().deltaCoeffs();

    tnGradU() -= (sqr(n) & tnGradU());
#endif

    return tnGradU;

    //return
    //(
    //    *this
    //  - (patchInternalField() + (k&gradU.patchInternalField()))
    //)*this->patch().deltaCoeffs();
}


tmp<Field<vector> > elasticWallVelocityFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    bool secondOrder_ = false;

#ifdef OPENFOAMESIORFOUNDATION
    const word UName = internalField().name();
#else
    const word UName = dimensionedInternalField().name();
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
    }


    // First order
    //const vectorField dUP(k & gradU.patchInternalField());

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
}



void elasticWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    elasticWallVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
