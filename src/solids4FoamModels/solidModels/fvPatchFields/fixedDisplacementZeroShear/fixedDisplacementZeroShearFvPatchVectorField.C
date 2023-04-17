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

#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_(),
    forceZeroShearGrad_(false)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(pvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    totalDisp_(mapper(pvf.totalDisp_)),
#else
    totalDisp_(pvf.totalDisp_, mapper),
#endif
    dispSeries_(pvf.dispSeries_),
    forceZeroShearGrad_(pvf.forceZeroShearGrad_)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    totalDisp_("value", dict, p.size()),
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
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        refValue() = dispSeries_(this->db().time().timeOutputValue());
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


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(pvf, iF),
    totalDisp_(pvf.totalDisp_),
    dispSeries_(pvf.dispSeries_),
    forceZeroShearGrad_(pvf.forceZeroShearGrad_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(totalDisp_, totalDisp_);;
#else
    totalDisp_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(pvf, addr);

    const fixedDisplacementZeroShearFvPatchVectorField& rpvf =
        refCast<const fixedDisplacementZeroShearFvPatchVectorField>(pvf);

    totalDisp_.rmap(rpvf.totalDisp_, addr);
}


void fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
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
void fixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    os.writeKeyword("forceZeroShearGrad")
        << forceZeroShearGrad_ << token::END_STATEMENT << nl;

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
