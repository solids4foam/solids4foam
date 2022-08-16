/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "fsiInterfaceSolidFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fsiInterfaceSolidFvPatchVectorField::updateForce()
{
    // Check if coupling switch needs to be updated
    if (!coupled_)
    {
        updateCoupled();
    }

    Info<< "Setting traction on solid interfaces" << endl;

    // for (label interfaceI = 0; interfaceI < nGlobalPatches_; interfaceI++)
    {
        // Take references to zones
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        // Calculate total traction of fluid zone
        vectorField fluidZoneTotalTraction
        (
            fluid().faceZoneViscousForce(interfaceI)
          - fluid().faceZonePressureForce(interfaceI)*fluidZone.faceNormals()
        );

        // Initialise the solid zone traction field that is to be interpolated
        // from the fluid zone
        vectorField solidZoneTotalTraction(solidZone.size(), vector::zero);

        // Transfer the field frm the fluid interface to the solid interface
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone,                 // from zone
            solidZone,                 // to zone
            fluidZoneTotalTraction,    // from field
            solidZoneTotalTraction     // to field
        );

        // Flip traction sign after transferring from fluid to solid
        solidZoneTotalTraction = -solidZoneTotalTraction;

        // Set traction on solid
        if (coupled())
        {
            solid().setTraction
            (
                interfaceI,
                solidPatchIndices()[interfaceI],
                solidZoneTotalTraction
            );
        }

        // Print total force on solid and fluid interfaces
        Info<< "Total force on fluid interface " << interfaceI << ": "
            << totalForceOnInterface(fluidZone, fluidZoneTotalTraction) << nl
            << "Total force on solid interface " << interfaceI << ": "
            << totalForceOnInterface(solidZone, solidZoneTotalTraction) << nl
            << endl;

        // Set interface pressure for elasticWallPressure boundary condition
        const label fluidPatchID = fluidPatchIndices()[interfaceI];
        if
        (
            isA<elasticWallPressureFvPatchScalarField>
            (
                fluid().p().boundaryField()[fluidPatchID]
            )
        )
        {
            scalarField& prevPressure =
                const_cast<elasticWallPressureFvPatchScalarField&>
                (
                    refCast<const elasticWallPressureFvPatchScalarField>
                    (
                        fluid().p().boundaryField()[fluidPatchID]
                    )
                ).prevPressure();

            if (coupled())
            {
                prevPressure = fluid().patchPressureForce(fluidPatchID);
            }
            else
            {
                prevPressure = 0;
            }
        }
    }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fsiInterfaceSolidFvPatchVectorField::fsiInterfaceSolidFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF)
{}


Foam::fsiInterfaceSolidFvPatchVectorField::fsiInterfaceSolidFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF)
{}


Foam::fsiInterfaceSolidFvPatchVectorField::fsiInterfaceSolidFvPatchVectorField
(
    const fsiInterfaceSolidFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::fsiInterfaceSolidFvPatchVectorField::fsiInterfaceSolidFvPatchVectorField
(
    const fsiInterfaceSolidFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf)
{}


Foam::fsiInterfaceSolidFvPatchVectorField::fsiInterfaceSolidFvPatchVectorField
(
    const fsiInterfaceSolidFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::fsiInterfaceSolidFvPatchVectorField::
~fsiInterfaceSolidFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fsiInterfaceSolidFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


void Foam::fsiInterfaceSolidFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);

    // const fsiInterfaceSolidFvPatchVectorField& dmptf =
    //     refCast<const fsiInterfaceSolidFvPatchVectorField>(ptf);
}

void Foam::fsiInterfaceSolidFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // if (curTimeIndex_ != db().time().timeIndex())
    // {
    //     // Update old quantities at the start of a new time-step
    //     curTimeIndex_ = db().time().timeIndex();
    // }

    solidTractionFvPatchVectorField::updateCoeffs();
}


void Foam::fsiInterfaceSolidFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fsiInterfaceSolidFvPatchVectorField
    );
}


// ************************************************************************* //
