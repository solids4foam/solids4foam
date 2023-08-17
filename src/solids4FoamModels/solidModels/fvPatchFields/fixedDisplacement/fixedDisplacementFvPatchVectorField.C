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

#include "fixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"
#include "patchCorrectionVectors.H"

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
    nonOrthogonalCorrections_(true),
    totalDisp_(p.size(), vector::zero),
    dispSeries_(),
    interpPtr_()
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
#ifdef OPENFOAMFOUNDATION
    totalDisp_(mapper(pvf.totalDisp_)),
#else
    totalDisp_(pvf.totalDisp_, mapper),
#endif
    dispSeries_(pvf.dispSeries_),
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
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    ),
    totalDisp_(*this),
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

#ifndef OPENFOAMFOUNDATION
fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    totalDisp_(pvf.totalDisp_),
    dispSeries_(pvf.dispSeries_),
    interpPtr_()
{}
#endif

fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_),
    totalDisp_(pvf.totalDisp_),
    dispSeries_(pvf.dispSeries_),
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
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);

    const fixedDisplacementFvPatchVectorField& rpvf =
        refCast<const fixedDisplacementFvPatchVectorField>(pvf);

    totalDisp_.rmap(rpvf.totalDisp_, addr);
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
    if (nonOrthogonalCorrections_)
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

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            *this
          - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        // fixedValue snGrad with no correction
        return (*this - patchInternalField())*patch().deltaCoeffs();
    }
}

tmp<Field<vector> >
fixedDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
    if (nonOrthogonalCorrections_)
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

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            patch().deltaCoeffs()
           *(*this - (k & gradField.patchInternalField()))
        );
    }
    else
    {
        return (patch().deltaCoeffs()*(*this));        
    }
}

void fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("nonOrthogonalCorrections")
        << nonOrthogonalCorrections_ << token::END_STATEMENT << nl;

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
