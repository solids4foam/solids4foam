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

#include "solidVelocityFvPatchVectorField.H"
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    velocity_(p.size(), vector::zero),
    velocitySeries_()
{}


solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    velocity_(ptf.velocity_),
    velocitySeries_(ptf.velocitySeries_)
{}


solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    velocity_(p.size(), vector::zero),
    velocitySeries_()
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // Read velocity
    if (dict.found("velocity"))
    {
        velocity_ = vectorField("velocity", dict, p.size());
    }
    else if (dict.found("velocitySeries"))
    {
        Info<< "    velocity is time-varying" << endl;
        velocitySeries_ =
            interpolationTable<vector>(dict.subDict("velocitySeries"));

        fvPatchField<vector>::operator==
        (
            velocitySeries_(this->db().time().timeOutputValue())
        );
    }
    else
    {
        FatalErrorIn(type() + "::solidVelocityFvPatchVectorField(...)")
            << "Either 'velocity' or 'velocitySeries' should be specified!"
            << abort(FatalError);
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}

#ifndef OPENFOAMFOUNDATION
solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& pivpvf
)
:
    fixedDisplacementFvPatchVectorField(pivpvf),
    velocity_(pivpvf.velocity_),
    velocitySeries_(pivpvf.velocitySeries_)
{}
#endif

solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(pivpvf, iF),
    velocity_(pivpvf.velocity_),
    velocitySeries_(pivpvf.velocitySeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedDisplacementFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(velocity_, velocity_);
#else
    velocity_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedDisplacementFvPatchVectorField::rmap(ptf, addr);

    const solidVelocityFvPatchVectorField& dmptf =
       refCast<const solidVelocityFvPatchVectorField>(ptf);

    velocity_.rmap(dmptf.velocity_, addr);
}


void solidVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Check if the velocity is time-varying
    if (velocitySeries_.size())
    {
        velocity_ = velocitySeries_(this->db().time().timeOutputValue());
    }

    vectorField disp = vectorField(patch().size(), vector::zero);

#ifdef OPENFOAMESIORFOUNDATION
    if (internalField().name() == "DD")
#else
    if (dimensionedInternalField().name() == "DD")
#endif
    {
        // Incremental approach, so we wil set the increment of displacement for
        // this time-step
        disp = velocity_*db().time().deltaTValue();
    }
    else
    {
        // Lookup the old time total displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        // The new total displacement is equal to Dold plus the increment of
        // displacement based on the current velocity and time-step
        disp =
            Dold.boundaryField()[patch().index()]
          + velocity_*db().time().deltaTValue();
    }

    // Set the displacement (or displacement increment) on the patch
    fvPatchField<vector>::operator==(disp);
    fixedValueFvPatchVectorField::updateCoeffs();

    // If the corresponding point displacement field has a fixedValue type
    // boundary condition, then we wil update it
    fixedDisplacementFvPatchVectorField::setPointDisplacement(disp);
}

void solidVelocityFvPatchVectorField::write(Ostream& os) const
{
    if (velocitySeries_.size())
    {
        os.writeKeyword("velocitySeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        velocitySeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "velocity", velocity_);
#else
        velocity_.writeEntry("velocity", os);
#endif
    }

    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    solidVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
