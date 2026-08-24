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

#include "tractionPressureDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fixedGradientFvPatchFields.H"
#include "lookupSolidModel.H"
#include "coupledPressureDisplacementSolid.H"
// #include "mixedNonLinGeomTotLagSolid.H"
// #include "newMixedNonLinGeomTotLagSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<tensorField> tractionPressureDisplacementFvPatchVectorField::FM
(
    const bool refConfig
) const
{
    tmp<tensorField> tFm
    (
        new tensorField(this->size(), tensor::I)
    );

    // const solidModel& solMod =
    //     lookupSolidModel(patch().boundaryMesh().mesh());

    if (!refConfig)
    {
        tensorField& Fm = tFm();

        const fvPatchField<tensor>& GradDD =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(DD)"
            );

        const fvPatchField<tensor>& GradD =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(D)"
            );

        const tensorField GradDm = GradD - 0.5*GradDD;

        Fm = I + GradDm.T();
    }

    return tFm;
}


tmp<tensorField> tractionPressureDisplacementFvPatchVectorField::intFM
(
    const bool refConfig
) const
{
    tmp<tensorField> tIntFm
    (
        new tensorField(this->size(), tensor::I)
    );

    // const solidModel& solMod =
    //     lookupSolidModel(patch().boundaryMesh().mesh());

    if (!refConfig)
    {
        tensorField& intFm = tIntFm();

        const fvPatchField<tensor>& GradDD =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(DD)"
            );

        const fvPatchField<tensor>& GradD =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(D)"
            );

        const tensorField GradDm =
            GradD.patchInternalField()
          - 0.5*GradDD.patchInternalField();

        intFm = I + GradDm.T();
    }

    return tIntFm;
}


tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
dispCorr(const bool refConfig) const
{
    tmp<vectorField> tDispCorr
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& DispCorr = tDispCorr();

    // Normal vector
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Non-orthogonal correction vectors
    vectorField k = ((I - sqr(n)) & delta);

    word DDName =
        this->dimensionedInternalField().name();

    // Disp incr gradient calc on ref configuration
    const fvPatchField<tensor>& GradDD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DDName + ")"
        );

    const tensorField invFm = inv(intFM(refConfig));

    // Disp incr gradient calculated at mean configuration
    tensorField gradDDP =
    (
        invFm.T() & GradDD.patchInternalField()
    );

    DispCorr = (k & gradDDP);

    return tDispCorr;
}


tmp<vectorField> tractionPressureDisplacementFvPatchVectorField::
patchInternalSnGrad(const bool refConfig) const
{
    tmp<vectorField> tPatchInternalSnGrad
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& pIntSnGrad = tPatchInternalSnGrad();

    // Normal vector
    const vectorField n = patch().nf();

    word DDName =
        this->dimensionedInternalField().name();

    const fvPatchField<tensor>& pGradDD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DDName + ")"
        );

    const tensorField invFm = inv(intFM(refConfig));

    tensorField gradDDP =
    (
        invFm.T() & pGradDD.patchInternalField()
    );

    pIntSnGrad = (n & gradDDP);

    return tPatchInternalSnGrad;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondaryGradient_(p.size(), vector::zero),
    rhieChowDisplCorr_(p.size(), 0),
    blockCoupled_(false),
    iPointsIndices_(),
    iPoints_(),
    iMatrices_(),
    PP_(p.size(), -1),
    frozenTraction_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const tractionPressureDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper),
    tractionSeries_(tdpvf.tractionSeries_),
    pressureSeries_(tdpvf.pressureSeries_),
    secondaryGradient_(tdpvf.secondaryGradient_, mapper),
    rhieChowDisplCorr_(tdpvf.rhieChowDisplCorr_, mapper),
    blockCoupled_(false),
    iPointsIndices_(),
    iPoints_(),
    iMatrices_(),
    PP_(p.size(), -1),
    frozenTraction_(false)
{}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondaryGradient_(p.size(), vector::zero),
    rhieChowDisplCorr_(p.size(), 0),
    blockCoupled_(dict.lookupOrDefault<bool>("blockCoupled", false)),
    iPointsIndices_(),
    iPoints_(),
    iMatrices_(),
    PP_(p.size(), -1),
    frozenTraction_(dict.lookupOrDefault<bool>("frozenTraction", false))
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    // Check if traction is time-varying
    if (dict.found("tractionSeries"))
    {
        if (debug)
        {
            Info<< "    traction is time-varying" << endl;
        }

        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));

        traction_ =
            tractionSeries_(db().time().timeOutputValue());
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        if (debug)
        {
            Info<< "    pressure is time-varying" << endl;
        }

        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));

        pressure_ =
            pressureSeries_(db().time().timeOutputValue());
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const tractionPressureDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    tractionSeries_(tdpvf.tractionSeries_),
    pressureSeries_(tdpvf.pressureSeries_),
    secondaryGradient_(tdpvf.secondaryGradient_),
    rhieChowDisplCorr_(tdpvf.rhieChowDisplCorr_),
    blockCoupled_(tdpvf.blockCoupled_),
    iPointsIndices_(),
    iPoints_(),
    iMatrices_(),
    PP_(tdpvf.PP_),
    frozenTraction_(false)
{}


tractionPressureDisplacementFvPatchVectorField::
tractionPressureDisplacementFvPatchVectorField
(
    const tractionPressureDisplacementFvPatchVectorField& tdpvf,

    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    tractionSeries_(tdpvf.tractionSeries_),
    pressureSeries_(tdpvf.pressureSeries_),
    secondaryGradient_(tdpvf.secondaryGradient_),
    rhieChowDisplCorr_(tdpvf.rhieChowDisplCorr_),
    blockCoupled_(tdpvf.blockCoupled_),
    iPointsIndices_(),
    iPoints_(),
    iMatrices_(),
    PP_(tdpvf.PP_),
    frozenTraction_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionPressureDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
    secondaryGradient_.autoMap(m);
    rhieChowDisplCorr_.autoMap(m);
    PP_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionPressureDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionPressureDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionPressureDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);

    tractionSeries_ = dmptf.tractionSeries_;
    pressureSeries_ = dmptf.pressureSeries_;

    secondaryGradient_.rmap(dmptf.secondaryGradient_, addr);
    rhieChowDisplCorr_.rmap(dmptf.rhieChowDisplCorr_, addr);

    PP_.rmap(dmptf.PP_, addr);

    blockCoupled_ = dmptf.blockCoupled_;
}


// Update the coefficients associated with the patch field
void tractionPressureDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    const solidModel& solMod = lookupSolidModel(mesh);

    secondaryGradient_ =
        solMod.tractionBoundarySnGrad
        (
            traction_,
            pressure_,
            this->patch()
        );

    // Function updateCoeffs() will be called when
    // mesh is in initial (reference) configuration
    const bool refConfig = true;

    const vectorField n = patch().nf();

    const vectorField& refSf = patch().Sf();

    vectorField curSf = refSf;
    if
    (
        solMod.nonLinGeom() !=
        nonLinearGeometry::LINEAR_GEOMETRY
    )
    {
        const tensorField Fm = FM(false);
        const tensorField invFm = inv(Fm);
        const scalarField Jm = det(Fm);
        curSf = Jm*(invFm.T() & refSf);
    }
    const scalarField magCurSf = mag(curSf);
    const vectorField curN = curSf/magCurSf;

    // Take conservative displacement increment
    vectorField pDD =
        patch().lookupPatchField<surfaceVectorField, vector>("DDf");

    // Calculate disp increment gradient
    this->gradient() =
    (
        pDD - (patchInternalField() + dispCorr(refConfig))
    )*this->patch().deltaCoeffs();

    // Correct tangential component of gradient
    this->gradient() -= ((I - sqr(curN)) & this->gradient());
    this->gradient() += ((I - sqr(curN)) & secondaryGradient_);

    // Add 2nd order correction
    bool secondOrder = false;
    if (secondOrder)
    {
        this->gradient() *= 2;
        this->gradient() -= patchInternalSnGrad(refConfig);
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


tmp<vectorField>
tractionPressureDisplacementFvPatchVectorField::
transformedGradient() const
{
    // Transform derivative to current configuration

    tmp<vectorField> tTransGradient
    (
        new vectorField(this->size(), vector::zero)
    );
    vectorField& transGradient = tTransGradient();

    // const fvsPatchVectorField& refS =
    //     patch().lookupPatchField<surfaceVectorField, vector>("refSf");
    // const fvsPatchVectorField& refDelta =
    //     patch().lookupPatchField<surfaceVectorField, vector>("refDelta");
    // const vectorField N = refS/mag(refS);

    const vectorField n = patch().nf();
    scalarField deltaN = 1.0/patch().deltaCoeffs();
    // scalarField refDeltaN = (refDelta & N);

    // For linear model
    scalarField refDeltaN = deltaN;

    // Check if non-linear model is used and mesh is moved
    const solidModel& solMod =
        lookupSolidModel(patch().boundaryMesh().mesh());

    if
    (
        solMod.nonLinGeom() !=
        nonLinearGeometry::LINEAR_GEOMETRY
    )
    {
        const bool refConfig = false;
        const tensorField invFm = inv(FM(refConfig));
        refDeltaN = mag(invFm & (n*deltaN));
    }

    transGradient = gradient()*refDeltaN/deltaN;

    return tTransGradient;
}


//- Evaluate the patch field
void tractionPressureDisplacementFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    bool secondOrder = false;

    const vectorField transGradient = transformedGradient();

    // The evaluate() function is called on the current (mean) configuration
    const bool refConfig = false;

    if (secondOrder)
    {
        Field<vector>::operator=
        (
            patchInternalField() + dispCorr(refConfig)
          + 0.5*(gradient() + patchInternalSnGrad(refConfig))/
            patch().deltaCoeffs()
        );
    }
    else
    {
        Field<vector>::operator=
        (
            patchInternalField()
          + dispCorr(refConfig)
          + transGradient/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


tmp<Field<vector> > tractionPressureDisplacementFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    // return gradient();
    return secondaryGradient_;
}


tmp<Field<vector> > tractionPressureDisplacementFvPatchVectorField::
valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<vector> >
    (
        new Field<vector>(this->size(), pTraits<vector>::one)
    );
}


tmp<Field<vector> > tractionPressureDisplacementFvPatchVectorField::
valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    bool secondOrder = false;

    // Transform gradient() to current configuration
    const vectorField transGradient = transformedGradient();

    // The evaluate() function is called on the current (mean) configuration
    const bool refConfig = false;

    if (secondOrder)
    {
        return dispCorr(refConfig) +
            0.5*(transGradient + patchInternalSnGrad(refConfig))/
            this->patch().deltaCoeffs();
    }
    else
    {
        return dispCorr(refConfig) +
            transGradient/this->patch().deltaCoeffs();
    }
}


// Write
void tractionPressureDisplacementFvPatchVectorField::
write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    secondaryGradient_.writeEntry("secondaryGradient", os);

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionPressureDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
