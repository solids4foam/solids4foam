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
#include "volFields.H"
#include "IStringStream.H"
#include "patchCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dictionary solidSpringDashpotFvPatchVectorField::tractionDict
(
    const dictionary& dict
)
{
    if (dict.found("tractionSeries") || dict.found("tractionField"))
    {
        FatalIOErrorInFunction(dict)
            << "tractionSeries and tractionField are not allowed for a "
            << typeName << " condition: the traction is what it computes."
            << exit(FatalIOError);
    }

    dictionary d(dict);

    // A written field carries the last computed traction; it is recomputed
    // anyway, so it is discarded rather than refused, which keeps restarts
    // and decomposed cases working
    d.remove("traction");
    d.remove("useUndeformedArea");

    string extra("traction uniform (0 0 0); useUndeformedArea yes;");

    if
    (
        !d.found("pressure")
     && !d.found("pressureSeries")
     && !d.found("pressureField")
    )
    {
        extra += " pressure uniform 0;";
    }

    IStringStream is(extra);
    d.merge(dictionary(is));

    return d;
}


tmp<vectorField>
solidSpringDashpotFvPatchVectorField::extrapolatedDisplacement() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

#ifdef OPENFOAM_NOT_EXTEND
    const word gradName("grad(" + internalField().name() + ")");
#else
    const word gradName("grad(" + dimensionedInternalField().name() + ")");
#endif

    tmp<vectorField> tDb(new vectorField(patchInternalField()));

    // Before the solid model has created its gradient, e.g. while the fields
    // are being read, the cell value is the best available
    if (mesh.foundObject<volTensorField>(gradName))
    {
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>(gradName);

        const tensorField gradP(gradField.patchInternalField());
        const vectorField n(patch().nf());
        const vectorField k(patchCorrectionVectors(patch()));

        // Same extrapolation as solidTraction with extrapolateValue, which
        // the PETSc SNES solid models switch on for traction boundaries
#ifdef OPENFOAM_NOT_EXTEND
        tDb.ref() +=
#else
        tDb() +=
#endif
            (k & gradP) + (n & gradP)/patch().deltaCoeffs();
    }

    return tDb;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidSpringDashpotFvPatchVectorField::
solidSpringDashpotFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    stiffness_(0),
    damping_(0),
    normalOnly_(false)
{}


solidSpringDashpotFvPatchVectorField::
solidSpringDashpotFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF, tractionDict(dict)),
    stiffness_(readScalar(dict.lookup("stiffness"))),
    damping_(dict.lookupOrDefault<scalar>("damping", 0.0)),
    normalOnly_(dict.lookupOrDefault<Switch>("normalOnly", false))
{
    if (stiffness_ < 0 || damping_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "stiffness and damping must not be negative."
            << exit(FatalIOError);
    }

    Info<< "    spring-dashpot: stiffness = " << stiffness_
        << " Pa/m, damping = " << damping_ << " Pa s/m"
        << (normalOnly_ ? ", normal direction only" : "") << endl;
}


solidSpringDashpotFvPatchVectorField::
solidSpringDashpotFvPatchVectorField
(
    const solidSpringDashpotFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(pvf, p, iF, mapper),
    stiffness_(pvf.stiffness_),
    damping_(pvf.damping_),
    normalOnly_(pvf.normalOnly_)
{}


#ifndef OPENFOAM_ORG
solidSpringDashpotFvPatchVectorField::
solidSpringDashpotFvPatchVectorField
(
    const solidSpringDashpotFvPatchVectorField& pvf
)
:
    solidTractionFvPatchVectorField(pvf),
    stiffness_(pvf.stiffness_),
    damping_(pvf.damping_),
    normalOnly_(pvf.normalOnly_)
{}
#endif


solidSpringDashpotFvPatchVectorField::
solidSpringDashpotFvPatchVectorField
(
    const solidSpringDashpotFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(pvf, iF),
    stiffness_(pvf.stiffness_),
    damping_(pvf.damping_),
    normalOnly_(pvf.normalOnly_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidSpringDashpotFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Current and old-time boundary displacement
    const vectorField Db(extrapolatedDisplacement());

    const volVectorField& D =
        db().lookupObject<volVectorField>
        (
#ifdef OPENFOAM_NOT_EXTEND
            internalField().name()
#else
            dimensionedInternalField().name()
#endif
        );

    const vectorField& Db0 = D.oldTime().boundaryField()[patch().index()];

    const scalar dt = db().time().deltaTValue();

    // Velocity by a backward difference; zero before the first time step
    vectorField Ub(Db.size(), vector::zero);
    if (dt > VSMALL && damping_ > 0)
    {
        Ub = (Db - Db0)/dt;
    }

    // Nominal traction of the spring and dashpot
    if (normalOnly_)
    {
        const vectorField N(patch().nf());

        traction() = -N*(stiffness_*(N & Db) + damping_*(N & Ub));
    }
    else
    {
        traction() = -(stiffness_*Db + damping_*Ub);
    }

    // solidTraction applies it, together with any pressure
    solidTractionFvPatchVectorField::updateCoeffs();
}


void solidSpringDashpotFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("stiffness")
        << stiffness_ << token::END_STATEMENT << nl;
    os.writeKeyword("damping")
        << damping_ << token::END_STATEMENT << nl;
    os.writeKeyword("normalOnly")
        << normalOnly_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidSpringDashpotFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
