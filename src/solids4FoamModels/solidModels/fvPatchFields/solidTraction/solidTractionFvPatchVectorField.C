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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
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
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(false),
    setEffectiveTraction_(false),
    relaxFac_(1.0),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
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
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    setEffectiveTraction_
    (
        dict.lookupOrDefault<Switch>("setEffectiveTraction", false)
    ),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    curTimeIndex_(-1)
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check how traction is defined
    if
    (
        (dict.found("traction") && dict.found("tractionSeries"))
     || (dict.found("traction") && dict.found("tractionField"))
     || (dict.found("tractionSeries") && dict.found("tractionField"))
    )
    {
        FatalErrorIn
        (
            "solidTractionFvPatchVectorField::solidTractionFvPatchVectorField"
        )   << "Only one of traction, tractionSeries or tractionField can be "
            << "specified!"
            << abort(FatalError);
    }
    else if (dict.found("tractionSeries"))
    {
        Info<< "    traction is time-varying" << endl;
        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else if (dict.found("tractionField"))
    {
        Info<< "    traction is specified as a field" << endl;
        tractionFieldPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    word(dict.lookup("tractionField")),
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check how pressure is defined
    if
    (
        (dict.found("pressure") && dict.found("pressureSeries"))
     || (dict.found("pressure") && dict.found("pressureField"))
     || (dict.found("pressureSeries") && dict.found("pressureField"))
    )
    {
        FatalErrorIn
        (
            "solidTractionFvPatchVectorField::solidTractionFvPatchVectorField"
        )   << "Only one of pressure, pressureSeries or pressureField can be "
            << "specified!"
            << abort(FatalError);
    }
    else if (dict.found("pressureSeries"))
    {
        Info<< "    pressure is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else if (dict.found("pressureField"))
    {
        Info<< "    pressure is specified as a field" << endl;
        pressureFieldPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    word(dict.lookup("pressureField")),
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

    if (secondOrder_)
    {
        Info<< "    second order correction" << endl;
    }

    if (setEffectiveTraction_)
    {
        Info<< "    set effective traction" << endl;
    }

    if (relaxFac_ < 1.0)
    {
        Info<< "    relaxation factor: " << relaxFac_ << endl;
    }
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    traction_(mapper(stpvf.traction_)),
    pressure_(mapper(stpvf.pressure_)),
#else
    traction_(stpvf.traction_, mapper),
    pressure_(stpvf.pressure_, mapper),
#endif
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(stpvf.secondOrder_),
    setEffectiveTraction_(stpvf.setEffectiveTraction_),
    relaxFac_(stpvf.relaxFac_),
    curTimeIndex_(stpvf.curTimeIndex_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(stpvf.secondOrder_),
    setEffectiveTraction_(stpvf.setEffectiveTraction_),
    relaxFac_(stpvf.relaxFac_),
    curTimeIndex_(stpvf.curTimeIndex_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    tractionFieldPtr_(),
    pressureFieldPtr_(),
    secondOrder_(stpvf.secondOrder_),
    setEffectiveTraction_(stpvf.setEffectiveTraction_),
    relaxFac_(stpvf.relaxFac_),
    curTimeIndex_(stpvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(traction_, traction_);
    m(pressure_, pressure_);
#else
    traction_.autoMap(m);
    pressure_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const solidTractionFvPatchVectorField& dmptf =
        refCast<const solidTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Called once per time-step

        if (pressureFieldPtr_.valid())
        {
            // Force the pressure field boundary conditions to update
            const_cast<volScalarField&>
            (
                pressureFieldPtr_()
            ).correctBoundaryConditions();
        }

        if (tractionFieldPtr_.valid())
        {
            // Force the traction field boundary conditions to update
            const_cast<volVectorField&>
            (
                tractionFieldPtr_()
            ).correctBoundaryConditions();
        }
    }

    if (tractionFieldPtr_.valid())
    {
        traction_ = tractionFieldPtr_().boundaryField()[patch().index()];
    }
    else if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    if (pressureFieldPtr_.valid())
    {
        pressure_ = pressureFieldPtr_().boundaryField()[patch().index()];
    }
    else if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    scalarField press(pressure_);
    if (setEffectiveTraction_)
    {
        const fvPatchField<scalar>& p =
            patch().lookupPatchField<volScalarField, scalar>("p");

        // Remove the dynamic pressure component: this will force the effective
        // traction to be enforced rather than the total traction
        press -= p;
    }

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    gradient() =
        relaxFac_*solMod.tractionBoundarySnGrad
        (
            traction_, press, patch()
        )
      + (1.0 - relaxFac_)*gradient();

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Face unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Non-orthogonal correction vectors
    const vectorField k((I - sqr(n)) & delta);

    if (secondOrder_)
    {
        const vectorField dUP(k & gradField.patchInternalField());
        const vectorField nGradUP(n & gradField.patchInternalField());

        Field<vector>::operator=
        (
            patchInternalField()
          + dUP
          + 0.5*(gradient() + nGradUP)/patch().deltaCoeffs()
        );
    }
    else
    {

        Field<vector>::operator=
        (
            patchInternalField()
          + (k & gradField.patchInternalField())
          + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


void solidTractionFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //fixedGradientFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

    if (tractionFieldPtr_.valid())
    {
        os.writeKeyword("tractionField")
            << tractionFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "traction", traction_);
#else
        traction_.writeEntry("traction", os);
#endif
    }

    if (pressureFieldPtr_.valid())
    {
        os.writeKeyword("pressureField")
            << pressureFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "pressure", pressure_);
#else
        pressure_.writeEntry("pressure", os);
#endif
    }

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
    os.writeKeyword("setEffectiveTraction")
        << setEffectiveTraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << nl;

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
    writeEntry(os, "gradient", gradient());
#else
    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
