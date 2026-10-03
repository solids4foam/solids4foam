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

#include "solidSpringDashpotFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "compatibilityFunctions.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * //

namespace
{

// Dictionary for the base class, which requires refValue, refGradient and
// valueFraction; they are set in updateCoeffs here
dictionary baseDict(const dictionary& dict)
{
    dictionary d
    (
        IStringStream
        (
            "refValue uniform (0 0 0);"
            "refGradient uniform (0 0 0);"
            "valueFraction uniform (0 0 0 0 0 0);"
        )()
    );
    d.merge(dict);

    return d;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidSpringDashpotFvPatchVectorField::solidSpringDashpotFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    kNormal_(p.size(), 0.0),
    kTangential_(p.size(), 0.0),
    cNormal_(p.size(), 0.0),
    cTangential_(p.size(), 0.0),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{}


solidSpringDashpotFvPatchVectorField::solidSpringDashpotFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF, baseDict(dict)),
    kNormal_("kNormal", dict, p.size()),
    kTangential_(p.size(), 0.0),
    cNormal_(p.size(), 0.0),
    cTangential_(p.size(), 0.0),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // valueFraction below assumes the first-order boundary value
    if (dict.lookupOrDefault<Switch>("secondOrder", false))
    {
        FatalIOErrorIn
        (
            "solidSpringDashpotFvPatchVectorField::"
            "solidSpringDashpotFvPatchVectorField(...)",
            dict
        )   << "secondOrder is not supported on patch " << patch().name()
            << exit(FatalIOError);
    }

    if (dict.found("kTangential"))
    {
        kTangential_ = scalarField("kTangential", dict, p.size());
    }

    if (dict.found("cNormal"))
    {
        cNormal_ = scalarField("cNormal", dict, p.size());
    }

    if (dict.found("cTangential"))
    {
        cTangential_ = scalarField("cTangential", dict, p.size());
    }

    if (dict.found("traction"))
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    if (dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

    if (dict.found("weightField"))
    {
        const fvMesh& mesh = patch().boundaryMesh().mesh();

        const volScalarField weight
        (
            IOobject
            (
                word(dict.lookup("weightField")),
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        Info<< "    coefficients weighted by " << weight.name() << endl;

        const scalarField& w = weight.boundaryField()[patch().index()];

        kNormal_ *= w;
        kTangential_ *= w;
        cNormal_ *= w;
        cTangential_ *= w;
    }

    if
    (
        min(kNormal_) < 0 || min(kTangential_) < 0
     || min(cNormal_) < 0 || min(cTangential_) < 0
    )
    {
        FatalErrorIn
        (
            "solidSpringDashpotFvPatchVectorField::"
            "solidSpringDashpotFvPatchVectorField(...)"
        )   << "Spring and dashpot coefficients must be non-negative on patch "
            << patch().name() << abort(FatalError);
    }

    // The dashpot velocity is (D - D.oldTime())/deltaT, so deltaT must be a
    // physical time step
    if
    (
        (gMax(cNormal_) > 0 || gMax(cTangential_) > 0)
     && word(d2dt2SchemeCompat(patch().boundaryMesh().mesh(), "d2dt2(D)"))
     == "steadyState"
    )
    {
        WarningIn
        (
            "solidSpringDashpotFvPatchVectorField::"
            "solidSpringDashpotFvPatchVectorField(...)"
        )   << "Patch " << patch().name() << ": the dashpot uses "
            << "(D - D.oldTime())/deltaT, which is a physical velocity only "
            << "if deltaT is a physical time step. Set cNormal and "
            << "cTangential to 0 for pseudo-time load stepping." << endl;
    }

    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=(patchInternalField());
    }
}


solidSpringDashpotFvPatchVectorField::solidSpringDashpotFvPatchVectorField
(
    const solidSpringDashpotFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(pvf, p, iF, mapper),
#ifdef OPENFOAM_ORG
    kNormal_(mapper(pvf.kNormal_)),
    kTangential_(mapper(pvf.kTangential_)),
    cNormal_(mapper(pvf.cNormal_)),
    cTangential_(mapper(pvf.cTangential_)),
    traction_(mapper(pvf.traction_)),
    pressure_(mapper(pvf.pressure_))
#else
    kNormal_(pvf.kNormal_, mapper),
    kTangential_(pvf.kTangential_, mapper),
    cNormal_(pvf.cNormal_, mapper),
    cTangential_(pvf.cTangential_, mapper),
    traction_(pvf.traction_, mapper),
    pressure_(pvf.pressure_, mapper)
#endif
{}


solidSpringDashpotFvPatchVectorField::solidSpringDashpotFvPatchVectorField
(
    const solidSpringDashpotFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(pvf, iF),
    kNormal_(pvf.kNormal_),
    kTangential_(pvf.kTangential_),
    cNormal_(pvf.cNormal_),
    cTangential_(pvf.cTangential_),
    traction_(pvf.traction_),
    pressure_(pvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> solidSpringDashpotFvPatchVectorField::springDashpotTraction
(
    const vectorField& nPressure,
    const scalarField& areaRatio
) const
{
    const vectorField n(patch().nf());
    const symmTensorField nn(sqr(n));
    const scalar rDeltaT = 1.0/db().time().deltaTValue();

    const vectorField& DP = *this;
    const vectorField& DoldP =
        db().lookupObject<volVectorField>("D").oldTime().boundaryField()
        [
            patch().index()
        ];

    return tmp<vectorField>
    (
        new vectorField
        (
            traction_ - nPressure*pressure_
          - areaRatio
           *(
                ((kNormal_*nn + kTangential_*(I - nn)) & DP)
              + rDeltaT*((cNormal_*nn + cTangential_*(I - nn)) & (DP - DoldP))
            )
        )
    );
}


void solidSpringDashpotFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);

#ifdef OPENFOAM_ORG
    m(kNormal_, kNormal_);
    m(kTangential_, kTangential_);
    m(cNormal_, cNormal_);
    m(cTangential_, cTangential_);
    m(traction_, traction_);
    m(pressure_, pressure_);
#else
    kNormal_.autoMap(m);
    kTangential_.autoMap(m);
    cNormal_.autoMap(m);
    cTangential_.autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
#endif
}


void solidSpringDashpotFvPatchVectorField::rmap
(
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(pvf, addr);

    const solidSpringDashpotFvPatchVectorField& rpvf =
        refCast<const solidSpringDashpotFvPatchVectorField>(pvf);

    kNormal_.rmap(rpvf.kNormal_, addr);
    kTangential_.rmap(rpvf.kTangential_, addr);
    cNormal_.rmap(rpvf.cNormal_, addr);
    cTangential_.rmap(rpvf.cTangential_, addr);
    traction_.rmap(rpvf.traction_, addr);
    pressure_.rmap(rpvf.pressure_, addr);
}


void solidSpringDashpotFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    if
    (
        solMod.incremental()
     || solMod.nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN
    )
    {
        FatalErrorIn("solidSpringDashpotFvPatchVectorField::updateCoeffs()")
            << "solidSpringDashpot acts on the total displacement D with "
            << "reference-configuration normals, and does not support "
            << "incremental or updated Lagrangian solid models (patch "
            << patch().name() << ")" << abort(FatalError);
    }

    // Reference unit normal (mesh is not moved)
    const vectorField n(patch().nf());
    const symmTensorField nn(sqr(n));

    // Reference-to-current area ratio, as K and C act per reference area
    scalarField areaRatio(patch().size(), 1.0);

    if (solMod.nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        const tensorField& Finv =
            patch().lookupPatchField<volTensorField, tensor>("Finv");
        const scalarField& J =
            patch().lookupPatchField<volScalarField, scalar>("J");

        areaRatio = patch().magSf()/mag(J*Finv.T() & patch().Sf());
    }

    const scalar rDeltaT = 1.0/db().time().deltaTValue();

    // Old-time patch displacement
    const volVectorField& Dold =
        db().lookupObject<volVectorField>("D").oldTime();
    const vectorField& DoldP = Dold.boundaryField()[patch().index()];

    // Explicit part of the dashpot traction, with dD/dt = (D - Dold)/dt
    const vectorField tractionExp
    (
        traction_
      + rDeltaT*areaRatio*((cNormal_*nn + cTangential_*(I - nn)) & DoldP)
    );

    const scalarField kEffNormal
    (
        areaRatio*(kNormal_ + rDeltaT*cNormal_)
    );
    const scalarField kEffTangential
    (
        areaRatio*(kTangential_ + rDeltaT*cTangential_)
    );

    // Gradient for the explicit traction
    refGrad() =
        solMod.tractionBoundarySnGrad(tractionExp, pressure_, patch());

    // 1/impK, as tractionBoundarySnGrad is affine in the traction
    const scalarField rImpK
    (
        n
      & (
            solMod.tractionBoundarySnGrad(tractionExp + n, pressure_, patch())
          - refGrad()
        )
    );

    // Implicit spring-dashpot: valueFraction = Keff/(Keff + impK*deltaCoeffs)
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    valueFraction() =
        kEffNormal*rImpK/(kEffNormal*rImpK + deltaCoeffs)*nn
      + kEffTangential*rImpK/(kEffTangential*rImpK + deltaCoeffs)*(I - nn);

    refValue() = vector::zero;

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


tmp<Field<vector> >
solidSpringDashpotFvPatchVectorField::snGradTransformDiag() const
{
    const symmTensorField& f = valueFraction();

    vectorField diag(f.size());
    diag.replace(vector::X, f.component(symmTensor::XX));
    diag.replace(vector::Y, f.component(symmTensor::YY));
    diag.replace(vector::Z, f.component(symmTensor::ZZ));

    return tmp<Field<vector> >(new vectorField(diag));
}


void solidSpringDashpotFvPatchVectorField::write(Ostream& os) const
{
#ifdef OPENFOAM_ORG
    writeEntry(os, "kNormal", kNormal_);
    writeEntry(os, "kTangential", kTangential_);
    writeEntry(os, "cNormal", cNormal_);
    writeEntry(os, "cTangential", cTangential_);
    writeEntry(os, "traction", traction_);
    writeEntry(os, "pressure", pressure_);
#else
    kNormal_.writeEntry("kNormal", os);
    kTangential_.writeEntry("kTangential", os);
    cNormal_.writeEntry("cNormal", os);
    cTangential_.writeEntry("cTangential", os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
#endif

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidSpringDashpotFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
