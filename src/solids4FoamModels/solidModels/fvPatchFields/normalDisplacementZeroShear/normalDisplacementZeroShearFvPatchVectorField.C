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

#include "normalDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

normalDisplacementZeroShearFvPatchVectorField::
normalDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    normalDisp_(p.size(), 0.0),
    dispSeries_(),
    forceZeroShearGrad_(false)
{}


normalDisplacementZeroShearFvPatchVectorField::
normalDisplacementZeroShearFvPatchVectorField
(
    const normalDisplacementZeroShearFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(pvf, p, iF, mapper),
#ifdef OPENFOAM_ORG
    normalDisp_(mapper(pvf.normalDisp_)),
#else
    normalDisp_(pvf.normalDisp_, mapper),
#endif
    dispSeries_(pvf.dispSeries_),
    forceZeroShearGrad_(pvf.forceZeroShearGrad_)
{}


normalDisplacementZeroShearFvPatchVectorField::
normalDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    normalDisp_(p.size(), 0.0),
    dispSeries_(),
    forceZeroShearGrad_
    (
        dict.lookupOrDefault<Switch>("forceZeroShearGrad", false)
    )
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<scalar>(dict.subDict("displacementSeries"));

        refValue() =
            patch().nf()*dispSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        normalDisp_ = scalarField("normalDisplacement", dict, p.size());

        refValue() = patch().nf()*normalDisp_;
    }

    this->refGrad() = vector::zero;

    this->valueFraction() = sqr(patch().nf());

    Field<vector> normalValue
    (
        transform(valueFraction(), refValue())
    );

    Field<vector> gradValue
    (
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs()
    );

    Field<vector> transformGradValue
    (
        transform(I - valueFraction(), gradValue)
    );

    Field<vector>::operator=(normalValue + transformGradValue);
}


normalDisplacementZeroShearFvPatchVectorField::
normalDisplacementZeroShearFvPatchVectorField
(
    const normalDisplacementZeroShearFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(pvf, iF),
    normalDisp_(pvf.normalDisp_),
    dispSeries_(pvf.dispSeries_),
    forceZeroShearGrad_(pvf.forceZeroShearGrad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void normalDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);

#ifdef OPENFOAM_ORG
    m(normalDisp_, normalDisp_);;
#else
    normalDisp_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void normalDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(pvf, addr);

    const normalDisplacementZeroShearFvPatchVectorField& rpvf =
        refCast<const normalDisplacementZeroShearFvPatchVectorField>(pvf);

    normalDisp_.rmap(rpvf.normalDisp_, addr);
}


void normalDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField nDisp = normalDisp_;

    if (dispSeries_.size())
    {
        nDisp = dispSeries_(this->db().time().timeOutputValue());
    }

    vectorField disp(nDisp*patch().nf());

#ifdef OPENFOAM_NOT_EXTEND
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

    // Set displacement
    refValue() = disp;

    // Set gradient to zero to force zero shear traction
    if (forceZeroShearGrad_)
    {
        refGrad() = vector::zero;
    }
    else
    {
        // Calculate the shear gradient such that the shear traction is zero

        // Lookup the solidModel object
        const solidModel& solMod =
            lookupSolidModel(patch().boundaryMesh().mesh());

        // Set gradient to force zero shear traction
        refGrad() =
            solMod.tractionBoundarySnGrad
            (
                vectorField(patch().size(), vector::zero),
                scalarField(patch().size(), 0.0),
                patch()
            );
    }

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void normalDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAM_ORG
        writeEntry(os, "normalDisp", normalDisp_);
#else
        normalDisp_.writeEntry("normalDisp", os);
#endif
    }

    os.writeKeyword("forceZeroShearGrad")
        << forceZeroShearGrad_ << token::END_STATEMENT << nl;

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    normalDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
