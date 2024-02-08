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

#include "solidForceFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidForceFvPatchVectorField::
solidForceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceFieldPtr_(),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidForceFvPatchVectorField::
solidForceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceFieldPtr_(),
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

    // Check how force is defined
    if (dict.found("force") && dict.found("forceField"))
    {
        FatalErrorIn
        (
            "solidForceFvPatchVectorField::solidForceFvPatchVectorField"
        )   << "Only force or forceField can be "
            << "specified, not both!"
            << abort(FatalError);
    }
    else if (dict.found("forceField"))
    {
        Info<< "    force is specified as a field" << endl;
        forceFieldPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    word(dict.lookup("forceField")),
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
        force_ = vectorField("force", dict, p.size());
    }
}


solidForceFvPatchVectorField::
solidForceFvPatchVectorField
(
    const solidForceFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAM_ORG
    force_(mapper(stpvf.force_)),
#else
    force_(stpvf.force_, mapper),
#endif
    forceFieldPtr_(),
    curTimeIndex_(stpvf.curTimeIndex_)
{}

#ifndef OPENFOAM_ORG
solidForceFvPatchVectorField::solidForceFvPatchVectorField
(
    const solidForceFvPatchVectorField& stpvf
)
:
    solidTractionFvPatchVectorField(stpvf),
    force_(stpvf.force_),
    forceFieldPtr_(),
    curTimeIndex_(stpvf.curTimeIndex_)
{}
#endif

solidForceFvPatchVectorField::solidForceFvPatchVectorField
(
    const solidForceFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(stpvf, iF),
    force_(stpvf.force_),
    forceFieldPtr_(),
    curTimeIndex_(stpvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidForceFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);

#ifdef OPENFOAM_ORG
    m(force_, force_);
#else
    force_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidForceFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);

    const solidForceFvPatchVectorField& dmptf =
        refCast<const solidForceFvPatchVectorField>(ptf);

    force_.rmap(dmptf.force_, addr);
}


// Update the coefficients associated with the patch field
void solidForceFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Called once per time-step
        if (forceFieldPtr_.valid())
        {
            // Force the traction field boundary conditions to update
            const_cast<volVectorField&>
            (
                forceFieldPtr_()
            ).correctBoundaryConditions();
        }
    }

    if (forceFieldPtr_.valid())
    {
        force_ = forceFieldPtr_().boundaryField()[patch().index()];
    }

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // Convert the force field to a traction field

    // Check if it is a linear or nonlinear geometry case
    if (solMod.nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
    {
        traction() = force_/patch().magSf();
    }
    else if (solMod.nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        // Patch area vectors
        const vectorField& patchSf = patch().Sf();

        // Lookup the inverse of the deformation gradient
        const tensorField& Finv =
            patch().lookupPatchField<volTensorField, tensor>("Finv");

        // Lookup the Jacobian
        const scalarField& J =
            patch().lookupPatchField<volScalarField, scalar>("J");

        // Calculate area vectors in the deformed configuration
        const scalarField patchDeformMagSf(mag(J*Finv.T() & patchSf));

        traction() = force_/patchDeformMagSf;
    }
    else if (solMod.nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Patch area vectors
        const vectorField& patchSf = patch().Sf();

        // Lookup the inverse of the relative deformation gradient
        const tensorField& relFinv =
            patch().lookupPatchField<volTensorField, tensor>("relFinv");

        // Lookup the relative Jacobian
        const scalarField& relJ =
            patch().lookupPatchField<volScalarField, scalar>("relJ");

        // Calculate area vectors in the deformed configuration
        const scalarField patchDeformMagSf
        (
            mag(relJ*relFinv.T() & patchSf)
        );

        traction() = force_/patchDeformMagSf;
    }
    else
    {
        FatalErrorIn("solidForceFvPatchVectorField::updateCoeffs()")
            << "Unknown solidModel nonLinGeom type = "
            << solMod.nonLinGeom() << abort(FatalError);
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}

void solidForceFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //solidTractionFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

    if (forceFieldPtr_.valid())
    {
        os.writeKeyword("forceField")
            << forceFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else
    {
#ifdef OPENFOAM_ORG
        writeEntry(os, "force", force_);
#else
        force_.writeEntry("force", os);
#endif
    }

#ifdef OPENFOAM_ORG
    writeEntry(os, "value", *this);
    writeEntry(os, "gradient", gradient());
#else
    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidForceFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
