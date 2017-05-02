/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

InClass
    solidContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::solidContactFvPatchVectorField::movingMesh() const
{
    // If the deformation gradient "F" and the displacement increment DD" are
    // found then we can assume it is a moving mesh (updated Lagrangian) case
    if (db().foundObject<volVectorField>("DD") && nonLinearGeometry())
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::solidContactFvPatchVectorField::nonLinearGeometry() const
{
    // If the deformation gradient "F" is found then we will assume the case to
    // be nonlinear geometry (i.e. finite strain)
    if (db().foundObject<volTensorField>("F"))
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::solidContactFvPatchVectorField::moveZonesToDeformedConfiguration()
{
    // Only the master moves the zones
    if (!master_)
    {
        return;
    }

    // Reference to mesh for tidiness
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Method
    // We will interpolate the patch face displacements to the patch vertices
    // and then add these vertex/point displacements to the initial patch
    // points
    // We need to take care in parallel, and also realise that the solidModel
    // might have a moving or stationary mesh

    // Assemble the zone face displacement field to move the zones
    vectorField zoneD(zone().size(), vector::zero);
    vectorField shadowZoneD(shadowZone().size(), vector::zero);

    // For a non-moving mesh, we will move the zones by the total
    // displacement, whereas for a moving mesh (updated Lagrangian), we will
    // move the zones by the displacement increment

    if (movingMesh())
    {
        // Updated Lagrangian, so we will move the zones by the displacement
        // increment

        // Lookup the current total displacement field
        const volVectorField& DD = db().lookupObject<volVectorField>("DD");

        // Take a reference to the patch face displacement increment field
        const vectorField& patchDD =
            DD.boundaryField()[patch().index()];
        const vectorField& shadowPatchDD =
            DD.boundaryField()[shadowPatchIndex()];

        zoneD =
            zoneField(zoneIndex(), patch().index(), patchDD);
        shadowZoneD =
            zoneField(shadowZoneIndex(), shadowPatchIndex(), shadowPatchDD);
    }
    else
    {
        // Non-moving mesh: we will move the zones by the total displacement

        // Lookup the current total displacement field
        const volVectorField& D = db().lookupObject<volVectorField>("D");

        // Take a reference to the patch face total displacement field
        const vectorField& patchD =
            D.boundaryField()[patch().index()];
        const vectorField& shadowPatchD =
            D.boundaryField()[shadowPatchIndex()];

        zoneD =
            zoneField(zoneIndex(), patch().index(), patchD);
        shadowZoneD =
            zoneField(shadowZoneIndex(), shadowPatchIndex(), shadowPatchD);
    }

    // Interpolate the zone face field to the zone points
    const pointField zonePointD =
        zoneFaceToPointInterpolate(zoneIndex(), zoneD);
    const pointField shadowZonePointD =
        zoneFaceToPointInterpolate(shadowZoneIndex(), shadowZoneD);

    // The zone deformed points are the initial position plus the
    // displacement
    const pointField zoneNewPoints =
        mesh.faceZones()[zoneIndex()]().localPoints()
      + zonePointD;
    const pointField shadowZoneNewPoints =
        mesh.faceZones()[shadowZoneIndex()]().localPoints()
      + shadowZonePointD;

    // Move the zones

    // Remove zones weights
    zone().movePoints(zoneNewPoints);
    shadowZone().movePoints(shadowZoneNewPoints);

    // We need to use const_cast to move the standAlonePatch points as the
    // movePoints function only clears weights
    // Also, be careful to move the points are opposed to the localPoints
    const_cast<pointField&>(zone().points()) = zoneNewPoints;
    const_cast<pointField&>(shadowZone().points()) = shadowZoneNewPoints;
}


void Foam::solidContactFvPatchVectorField::calcZoneIndex() const
{
    if (zoneIndexPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::calcZoneIndex() const"
        )   << "zoneIndexPtr_ already set" << abort(FatalError);
    }

    // It is assumed that the faceZone correspondng to the patch has the same
    // name as the patch with the "FaceZone" sufix
    const word zoneName = patch().name() + "FaceZone";
    const faceZoneID zone(zoneName, patch().boundaryMesh().mesh().faceZones());

    if (!zone.active())
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::calcZoneIndex() const"
        )   << "Face zone name " << zoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    zoneIndexPtr_ = new label(zone.index());
}


Foam::label Foam::solidContactFvPatchVectorField::zoneIndex() const
{
    if (!zoneIndexPtr_)
    {
        calcZoneIndex();
    }

    return *zoneIndexPtr_;
}


void Foam::solidContactFvPatchVectorField::calcShadowZoneIndex() const
{
    if (shadowZoneIndexPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcShadowZoneIndex() const"
        )   << "shadowZoneIndexPtr_ already set" << abort(FatalError);
    }

    // It is assumed that the faceZone correspondng to the patch has the same
    // name as the patch with the "FaceZone" sufix
    const word shadowZoneName = shadowPatchName_ + "FaceZone";
    const faceZoneID shadowZone
    (
        shadowZoneName, patch().boundaryMesh().mesh().faceZones()
    );

    if (!shadowZone.active())
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcShadowZoneIndex() const"
        )   << "Face zone name " << shadowZoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    shadowZoneIndexPtr_ = new label(shadowZone.index());
}


Foam::label Foam::solidContactFvPatchVectorField::shadowZoneIndex() const
{
    if (!shadowZoneIndexPtr_)
    {
        calcShadowZoneIndex();
    }

    return *shadowZoneIndexPtr_;
}


void Foam::solidContactFvPatchVectorField::calcNormalModel() const
{
    if (normalModelPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcNormalModel() const"
        )   << "pointer already set" << abort(FatalError);
    }

    normalModelPtr_ =
        normalContactModel::New
        (
            word(dict().lookup("normalContactModel")),
            patch(),
            dict(),
            patch().index(),       // master
            shadowPatchIndex(),    // slave
            zoneIndex(),           // master face zone ID
            shadowZoneIndex(),     // slave face zone ID
            zone(),
            shadowZone()
        ).ptr();
}


Foam::normalContactModel& Foam::solidContactFvPatchVectorField::normalModel()
{
    if (master())
    {
        if (!normalModelPtr_)
        {
            calcNormalModel();
        }

        return *normalModelPtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.normalModel();
}


const Foam::normalContactModel&
Foam::solidContactFvPatchVectorField::normalModel() const
{
    if (master())
    {
        if (!normalModelPtr_)
        {
            calcNormalModel();
        }

        return *normalModelPtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.normalModel();
}


void Foam::solidContactFvPatchVectorField::calcFrictionModel() const
{
    if (frictionModelPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::"
            "calcFrictionModel() const"
        )   << "pointer already set" << abort(FatalError);
    }

    frictionModelPtr_ =
        frictionContactModel::New
        (
            word(dict().lookup("frictionContactModel")),
            patch(),
            dict(),
            patch().index(),       // master
            shadowPatchIndex(),    // slave
            zoneIndex(),           // master face zone ID
            shadowZoneIndex()      // slave face zone ID
        ).ptr();
}


Foam::frictionContactModel&
Foam::solidContactFvPatchVectorField::frictionModel()
{
    if (master())
    {
        if (!frictionModelPtr_)
        {
            calcFrictionModel();
        }

        return *frictionModelPtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.frictionModel();
}


const Foam::frictionContactModel&
Foam::solidContactFvPatchVectorField::frictionModel() const
{
    if (master())
    {
        if (!frictionModelPtr_)
        {
            calcFrictionModel();
        }

        return *frictionModelPtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.frictionModel();
}


void Foam::solidContactFvPatchVectorField::calcZone() const
{
    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::calcZone() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ =
        new standAlonePatch
        (
            mesh.faceZones()[zoneIndex()]().localFaces(),
            mesh.faceZones()[zoneIndex()]().localPoints()
        );
}


const Foam::standAlonePatch&
Foam::solidContactFvPatchVectorField::zone() const
{
    if (master())
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.zone();
}


Foam::standAlonePatch& Foam::solidContactFvPatchVectorField::zone()
{
    if (master())
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.zone();
}


void Foam::solidContactFvPatchVectorField::calcShadowZone() const
{
    if (shadowZonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactFvPatchVectorField::calcShadowZone() const"
        )   << "shadowZonePtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    shadowZonePtr_ =
        new standAlonePatch
        (
            mesh.faceZones()[shadowZoneIndex()]().localFaces(),
            mesh.faceZones()[shadowZoneIndex()]().localPoints()
        );
}


const Foam::standAlonePatch&
Foam::solidContactFvPatchVectorField::shadowZone() const
{
    if (master())
    {
        if (!shadowZonePtr_)
        {
            calcShadowZone();
        }

        return *shadowZonePtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.shadowZone();
}


Foam::standAlonePatch& Foam::solidContactFvPatchVectorField::shadowZone()
{
    if (master())
    {
        if (!shadowZonePtr_)
        {
            calcShadowZone();
        }

        return *shadowZonePtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.shadowZone();
}


void Foam::solidContactFvPatchVectorField::calcZoneToZone() const
{
    // Create zone-to-zone interpolation
    if (zoneToZonePtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::calcZoneToZone() const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    // Check master and slave patch
    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    if (master())
    {
        if (shadowPatchField.master() == true)
        {
            FatalErrorIn("solidContactFvPatchScalarField")
                << "There are two master patches"
                << abort(FatalError);
        }
    }
    else
    {
        if (shadowPatchField.master() == false)
        {
            FatalErrorIn("solidContactFvPatchScalarField")
                << "There is no master patch"
                << abort(FatalError);
        }
    }

    if (master())
    {
        // Create interpolation for patches
        zoneToZonePtr_ =
            new extendedGgiStandAlonePatchInterpolation
            (
                zone(),
                shadowZone(),
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation. Bug fix
                0,              // Master non-overlapping face tolerances
                0,              // Slave non-overlapping face tolerances
                true,           // Rescale weighting factors
                quickReject_ //,
                //regionOfInterest_
            );

        //Info<< "Region of interest: " << regionOfInterest_ << endl;
    }
    else
    {
        FatalErrorIn
        (
            "void solidContactFvPatchVectorField::calcZoneToZone() const"
        )   << "Attempting to create GGIInterpolation on a slave"
            << abort(FatalError);
    }
}


const Foam::extendedGgiStandAlonePatchInterpolation&
Foam::solidContactFvPatchVectorField::zoneToZone() const
{
    if (master())
    {
        if (!zoneToZonePtr_)
        {
            if (debug)
            {
                word zoneName =
                    patch().boundaryMesh().mesh().faceZones()
                    [
                        zoneIndex()
                    ].name();

                word shadowZoneName =
                    patch().boundaryMesh().mesh()
                    .faceZones()[shadowZoneIndex()].name();

                Info<< "Initializing the GGI interpolator between "
                    << "master/shadow zones: "
                    << zoneName << "/" << shadowZoneName
                    << endl;
            }

            calcZoneToZone();
        }

        return *zoneToZonePtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.zoneToZone();
}


Foam::extendedGgiStandAlonePatchInterpolation&
Foam::solidContactFvPatchVectorField::zoneToZone()
{
    if (master())
    {
        if (!zoneToZonePtr_)
        {
            word zoneName =
                patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

            word shadowZoneName =
                patch().boundaryMesh().mesh()
               .faceZones()[shadowZoneIndex()].name();

            if (debug)
            {
                Info<< "Initializing the GGI interpolator between "
                    << "master/shadow zones: "
                    << zoneName << "/" << shadowZoneName
                    << endl;
            }

            calcZoneToZone();
        }

        return *zoneToZonePtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    solidContactFvPatchVectorField& shadowPatchField =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndex_]
            )
        );

    return shadowPatchField.zoneToZone();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    master_(false),
    writeZoneVTK_(false),
    shadowPatchName_("undefined"),
    shadowPatchIndex_(-1),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(false),
    dict_(NULL),
    normalModelPtr_(NULL),
    frictionModelPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_(Foam::extendedGgiInterpolation::AABB),
    //regionOfInterest_(vector::min, vector::max),
    curTimeIndex_(-1)
{}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModelPtr_(NULL),
    frictionModelPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_(ptf.quickReject_),
    //regionOfInterest_(ptf.regionOfInterest_),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Copy pointer objects

    if (ptf.zoneIndexPtr_)
    {
        zoneIndexPtr_ = new label(*ptf.zoneIndexPtr_);
    }

    if (ptf.shadowZoneIndexPtr_)
    {
        shadowZoneIndexPtr_ = new label(*ptf.shadowZoneIndexPtr_);
    }

    if (ptf.normalModelPtr_)
    {
        normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    }

    if (ptf.frictionModelPtr_)
    {
        frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new standAlonePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    master_(dict.lookup("master")),
    writeZoneVTK_(dict.lookupOrDefault<Switch>("writeZoneVTK", false)),
    shadowPatchName_(dict.lookup("shadowPatch")),
    shadowPatchIndex_(-1),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(false),
    dict_(dict),
    normalModelPtr_(NULL),
    frictionModelPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_
    (
        extendedGgiInterpolation::quickRejectNames_
        [
            dict.lookupOrDefault<word>("quickReject", "AABB")
        ]
    ),
    // regionOfInterest_
    // (
    //     dict.lookupOrDefault<boundBox>
    //     (
    //         "regionOfInterest",
    //         boundBox(vector::min, vector::max)
    //     )
    // ),
    curTimeIndex_(-1)
{
    Info<< "Creating " << solidContactFvPatchVectorField::typeName << " patch"
        << endl;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Shadow patch index
    polyPatchID shadow(shadowPatchName_, mesh.boundaryMesh());

    if (!shadow.active())
    {
        FatalErrorIn
        (
            "solidContactFvPatchScalarField::"
            "solidContactFvPatchScalarField(...)"
        )
            << "Shadow patch name " << shadowPatchName_ << " not found."
            << abort(FatalError);
    }

    shadowPatchIndex_ = shadow.index();

    // Master creates contact laws
    if (master_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));
    }

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
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
    master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModelPtr_(ptf.normalModelPtr_),
    frictionModelPtr_(ptf.frictionModelPtr_),
    zonePtr_(ptf.zonePtr_),
    shadowZonePtr_(ptf.shadowZonePtr_),
    zoneToZonePtr_(ptf.zoneToZonePtr_),
    quickReject_(ptf.quickReject_),
    //regionOfInterest_(ptf.regionOfInterest_),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Copy pointer objects

    if (ptf.zoneIndexPtr_)
    {
        zoneIndexPtr_ = new label(*ptf.zoneIndexPtr_);
    }

    if (ptf.shadowZoneIndexPtr_)
    {
        shadowZoneIndexPtr_ = new label(*ptf.shadowZoneIndexPtr_);
    }

    if (ptf.normalModelPtr_)
    {
        normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    }

    if (ptf.frictionModelPtr_)
    {
        frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new standAlonePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    master_(ptf.master_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModelPtr_(ptf.normalModelPtr_),
    frictionModelPtr_(ptf.frictionModelPtr_),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_(ptf.quickReject_),
    //regionOfInterest_(ptf.regionOfInterest_),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Copy pointer objects

    if (ptf.zoneIndexPtr_)
    {
        zoneIndexPtr_ = new label(*ptf.zoneIndexPtr_);
    }

    if (ptf.shadowZoneIndexPtr_)
    {
        shadowZoneIndexPtr_ = new label(*ptf.shadowZoneIndexPtr_);
    }

    if (ptf.normalModelPtr_)
    {
        normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    }

    if (ptf.frictionModelPtr_)
    {
        frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new standAlonePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solidContactFvPatchVectorField::~solidContactFvPatchVectorField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);

    // Force all data to be re-created when needed
    clearOut();
}


void Foam::solidContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


void Foam::solidContactFvPatchVectorField::clearOut()
{
    deleteDemandDrivenData(zoneIndexPtr_);
    deleteDemandDrivenData(shadowZoneIndexPtr_);
    deleteDemandDrivenData(normalModelPtr_);
    deleteDemandDrivenData(frictionModelPtr_);
    deleteDemandDrivenData(zonePtr_);
    deleteDemandDrivenData(shadowZonePtr_);
    deleteDemandDrivenData(zoneToZonePtr_);
}


void Foam::solidContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = this->db().time().timeIndex();

        // Let the contact models know that it is a new time-step, in case they
        // need to update anything
        normalModel().newTimeStep();
        frictionModel().newTimeStep();
    }

    // Move the master and slave zone to the deformed configuration
    moveZonesToDeformedConfiguration();

    // Delete the zone-to-zone interpolator weights as the zones have moved
    zoneToZone().movePoints(tensorField(0), tensorField(0), vectorField(0));

    // Calculate and apply contact forces
    if (master_)
    {
        vectorField shadowPatchFaceNormals
        (
            patch().boundaryMesh()[shadowPatchIndex()].size(),
            vector::zero
        );

        if (nonLinearGeometry())
        {
            // Use deformed configuration face unit normals
            shadowPatchFaceNormals =
                patchField
                (
                    shadowPatchIndex(),
                    shadowZoneIndex(),
                    shadowZone().faceNormals()
                );
        }
        else
        {
            // Use undeformed configuration face unit normals
            shadowPatchFaceNormals =
                patch().boundaryMesh()[shadowPatchIndex()].nf();
        }

        // Calculate normal contact forces
        normalModel().correct
        (
            shadowPatchFaceNormals,
            zoneToZone()
        );

        // Interpolate the master displacement increment to the slave patch as
        // it required by the friction model

        vectorField patchDD(patch().size(), vector::zero);
        vectorField shadowPatchDD
        (
            patch().boundaryMesh()[shadowPatchIndex()].size(),
            vector::zero
        );

        if (movingMesh())
        {
            // Updated Lagrangian, we will directly lookup the displacement
            // increment

            const volVectorField& DD = db().lookupObject<volVectorField>("DD");

            patchDD = DD.boundaryField()[patch().index()];
            shadowPatchDD = DD.boundaryField()[shadowPatchIndex()];
        }
        else
        {
            // We will lookup the total displacement and old total displacement

            const volVectorField& D = db().lookupObject<volVectorField>("D");

            patchDD =
                D.boundaryField()[patch().index()]
              - D.oldTime().boundaryField()[patch().index()];
            shadowPatchDD =
                D.boundaryField()[shadowPatchIndex()]
              - D.oldTime().boundaryField()[shadowPatchIndex()];
        }

        // Master zone DD
        const vectorField zoneDD =
            zoneField
            (
                zoneIndex(),
                patch().index(),
                patchDD
            );

        // Master patch DD interpolated to the slave patch
        const vectorField patchDDInterpToShadowPatch =
            patchField
            (
                shadowPatchIndex(),
                shadowZoneIndex(),
                zoneToZone().masterToSlave(zoneDD)()
            );

        // Calculate friction contact forces
        frictionModel().correct
        (
            normalModel().slavePressure(),
            shadowPatchFaceNormals,
            normalModel().areaInContact(),
            shadowPatchDD,
            patchDDInterpToShadowPatch
        );

        if (rigidMaster_)
        {
            // Set to master to traction free to mimic a rigid contact
            traction() = vector::zero;
        }
        else
        {
            // Interpolate slave traction to the master
            const vectorField slavePatchTraction =
               - frictionModel().slaveTractionForMaster()
               - normalModel().slavePressure();

            const vectorField slaveZoneTraction =
                zoneField
                (
                    shadowZoneIndex(),
                    shadowPatchIndex(),
                    slavePatchTraction
                );

            // We have two options for interpolating from the slave to the
            // master:
            // 1. face-to-face
            // 2. point-to-point

            const bool faceToFaceInterpolation_ = true;
            if (faceToFaceInterpolation_)
            {
                // Set the traction on the master patch
                traction() =
                    patchField
                    (
                        patch().index(),
                        zoneIndex(),
                        zoneToZone().slaveToMaster(slaveZoneTraction)()
                    );
            }
            else
            {
                // Interpolate the slave traction from the faces to the points
                const vectorField slaveZonePointTraction =
                    zoneFaceToPointInterpolate
                    (
                        shadowZoneIndex(), slaveZoneTraction
                    );

                // Interpolate the slave point tractions to the master points
                const vectorField masterZonePointTraction =
                    zoneToZone().slaveToMasterPointInterpolate
                    (
                        slaveZonePointTraction
                    );

                // Interpolate the master point tractions to the masater faces
                // and set the traction on the master patch
                traction() =
                    patchField
                    (
                        patch().index(),
                        zoneIndex(),
                        zonePointToFaceInterpolate
                        (
                            zoneIndex(), masterZonePointTraction
                        )()
                    );
            }
        }
    }
    else
    {
        // Set the traction on the slave patch
        traction() =
            frictionModel().slaveTraction() + normalModel().slavePressure();
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


// Foam::tmp<Foam::scalarField> Foam::solidContactFvPatchVectorField::Qc() const
// {
//     // Consider storing Qc instead of recalculating multiple times

//     if (!master())
//     {
//         FatalErrorIn
//         (
//             "Foam::tmp<Foam::scalarField> Foam::"
//             "solidContactFvPatchVectorField::Qc() const"
//         )   << "Only master can call Qc function!"
//             << abort(FatalError);
//     }

//     // For now, we assume traction is constant over time-step
//     // Todo: we should use trapezoidal rule
//     vectorField curTraction(patch().size(), vector::zero);

//     // sigma/sigmaCauchy is up-to-date as Qc is called after momentum loop
//     // has converged and sigma has been updated and mesh moved
//     if (db().foundObject<volSymmTensorField>("sigmaCauchy"))
//     {
//         const symmTensorField& sigma =
//             db().lookupObject<volSymmTensorField>
//             (
//                 "sigmaCauchy"
//             ).boundaryField()[patch().index()];

//         curTraction = patch().nf() & sigma;
//     }
//     else
//     {
//         const symmTensorField& sigma =
//             db().lookupObject<volSymmTensorField>
//             (
//                 "sigma"
//             ).boundaryField()[patch().index()];

//         curTraction = patch().nf() & sigma;
//     }

//     // Calculate slip

//     const vectorField slavePatchSlip = frictionModel().slip();

//     const vectorField slaveZoneSlip =
//         zoneField
//         (
//             shadowZoneIndex(),
//             shadowPatchIndex(),
//             slavePatchSlip
//         );

//     // Interpolate from slave to master

//     const vectorField masterZoneSlip =
//         zoneToZone().slaveToMaster(slaveZoneSlip);

//     const vectorField masterPatchSlip =
//         patchField
//         (
//             patch().index(),
//             zoneIndex(),
//             masterZoneSlip
//         );

//     const scalar deltaT = patch().boundaryMesh().mesh().time().deltaTValue();

//     // Increment of dissipated frictional energy for this timestep
//     // The dot product of the traction vectors and the slip vectors gives the
//     // dissipated frictional energy per unit area; which is always positive
//     return
//         tmp<scalarField>
//         (
//             new scalarField(mag(curTraction & (masterPatchSlip/deltaT)))
//         );
// }


void Foam::solidContactFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
    os.writeKeyword("shadowPatch")
        << patch().boundaryMesh().mesh().boundary()[shadowPatchIndex()].name()
        << token::END_STATEMENT << nl;
    // os.writeKeyword("regionOfInterest")
    //     << regionOfInterest_ << token::END_STATEMENT << nl;
    os.writeKeyword("writeZoneVTK")
        << writeZoneVTK_ << token::END_STATEMENT << nl;

    if (master_)
    {
        os.writeKeyword("rigidMaster") << rigidMaster_
            << token::END_STATEMENT << nl;

        os.writeKeyword("normalContactModel")
            << normalModel().type() << token::END_STATEMENT << nl;
        normalModel().writeDict(os);

        os.writeKeyword("frictionContactModel")
            << frictionModel().type() << token::END_STATEMENT << nl;
        frictionModel().writeDict(os);
    }

    if (writeZoneVTK_)
    {
        if
        (
            dimensionedInternalField().name() == "D"
         || dimensionedInternalField().name() == "DD"
        )
        {
            Info<< "Writing deformed zones to VTK" << endl;
            const word timeName =
                patch().boundaryMesh().mesh().time().timeName();
            zone().writeVTK("zone_" + timeName);
            shadowZone().writeVTK("shadowZone_" + timeName);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchVectorField, solidContactFvPatchVectorField);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace Foam

// ************************************************************************* //
