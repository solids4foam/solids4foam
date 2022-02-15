/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void fixedDisplacementFvPatchVectorField::makeInterp() const
{
    if (interpPtr_.valid())
    {
        FatalErrorIn
        (
            "void fixedDisplacementFvPatchVectorField::makeInterp() const"
        ) << "pointer already set" << abort(FatalError);
    }

    interpPtr_.set(new primitivePatchInterpolation(patch().patch()));
}


primitivePatchInterpolation& fixedDisplacementFvPatchVectorField::interp()
{
    if (interpPtr_.empty())
    {
        makeInterp();
    }

    return interpPtr_();
}


// * * * * * * * * *  Protected Member Functions  * * * * * * * * * * * * * * //

void fixedDisplacementFvPatchVectorField::setPointDisplacement
(
    const vectorField& faceDisp
)
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if
    (
        mesh.foundObject<pointVectorField>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "point" + internalField().name()
#else
            "point" + dimensionedInternalField().name()
#endif
        )
    )
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                "point" + internalField().name()
#else
                "point" + dimensionedInternalField().name()
#endif
            );

        // Check if the boundary is fixedValue
        if
        (
            pointD.boundaryField()[patch().index()].type()
         == fixedValuePointPatchVectorField::typeName
        )
        {
            // Use const_cast to set boundary condition
            fixedValuePointPatchVectorField& patchPointD =
                refCast<fixedValuePointPatchVectorField>
                (
                    const_cast<pointVectorField&>
                    (
                        pointD
#ifdef OPENFOAMESIORFOUNDATION
                    ).boundaryFieldRef()[patch().index()]
#else
                    ).boundaryField()[patch().index()]
#endif
                );

            // Interpolate face values to the points
            patchPointD == interp().faceToPointInterpolate(faceDisp);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_(),
    interpPtr_()
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    totalDisp_(mapper(ptf.totalDisp_)),
#else
    totalDisp_(ptf.totalDisp_, mapper),
#endif
    dispSeries_(ptf.dispSeries_),
    interpPtr_()
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    totalDisp_("value", dict, p.size()),
    dispSeries_(),
    interpPtr_()
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        fvPatchField<vector>::operator==
        (
            dispSeries_(this->db().time().timeOutputValue())
        );
    }
}


// fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
// (
//     const fixedDisplacementFvPatchVectorField& pivpvf
// )
// :
//     fixedValueFvPatchVectorField(pivpvf),
//     totalDisp_(pivpvf.totalDisp_),
//     dispSeries_(pivpvf.dispSeries_),
//     interpPtr_()
// {}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    totalDisp_(pivpvf.totalDisp_),
    dispSeries_(pivpvf.dispSeries_),
    interpPtr_()
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //


fixedDisplacementFvPatchVectorField::~fixedDisplacementFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(totalDisp_, totalDisp_);;
#else
    totalDisp_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const fixedDisplacementFvPatchVectorField& dmptf =
        refCast<const fixedDisplacementFvPatchVectorField>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp = totalDisp_;

    if (dispSeries_.size())
    {
        disp = dispSeries_(this->db().time().timeOutputValue());
    }

#ifdef OPENFOAMESIORFOUNDATION
    if (internalField().name() == "DD")
#else
    if (dimensionedInternalField().name() == "DD")
#endif
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    fvPatchField<vector>::operator==(disp);

    fixedValueFvPatchVectorField::updateCoeffs();

    // If the corresponding point displacement field has a fixedValue type
    // boundary condition, then we wil update it
    setPointDisplacement(disp);
}


Foam::tmp<Foam::Field<vector> >
fixedDisplacementFvPatchVectorField::snGrad() const
{
    //- fixedValue snGrad with no correction
    //  return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    //- Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    return
    (
        *this
        - (patchInternalField() + (k & gradField.patchInternalField()))
    )*patch().deltaCoeffs();
}

tmp<Field<vector> >
fixedDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Unit normals
    vectorField n(patch().nf());

    // Delta vectors
    vectorField delta(patch().delta());

    //- Non-orthogonal correction vectors
    vectorField k((I - sqr(n)) & delta);

    return
    (
        patch().deltaCoeffs()
       *(*this - (k & gradField.patchInternalField()))
    );
}

void fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
