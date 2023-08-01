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

InClass
    solidContactPointPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidContactPointPatchVectorField.H"
#include "pointFields.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::solidContactPointPatchVectorField::moveZonesToDeformedConfiguration()
{
    // Only the master moves the zones
    if (!master_)
    {
        return;
    }

    // Method
    // We add the vertex/point displacements to the initial patch points
    // We need to take care in parallel, and also realise that the solidModel
    // might have a moving or stationary mesh

    // Lookup the point displacement field
    // Assemble the zone face displacement field to move the zones

    // First we will move the master zone

    // For a non-moving mesh, we will move the zones by the total
    // displacement, whereas for a moving mesh (updated Lagrangian), we will
    // move the zones by the displacement increment

    // Lookup the current total displacement field
    const pointVectorField& pointD =
        db().lookupObject<pointVectorField>("pointD");

    // Patch total displacement
    vectorField patchD(patchInternalField());

    if (movingMesh())
    {
        // Subtract the old-time displacement
        patchD -=
            pointD.oldTime().boundaryField()
            [
                patch().index()
            ].patchInternalField();
    }

    // The zone deformed points are the initial position plus the
    // displacement
    const pointField zoneNewPoints
    (
        zone().patchPointToGlobal(pointField(patch().localPoints() + patchD))
    );

    // Remove zone weights
    zone().movePoints(zoneNewPoints);

    // We need to use const_cast to move the standAlonePatch points as the
    // movePoints function only clears weights
    // Also, be careful to move the points as opposed to the localPoints
    const_cast<pointField&>(zone().globalPatch().points()) = zoneNewPoints;

    // Secondly we will move the shadow zones

    forAll(shadowPatchNames(), shadPatchI)
    {
        // Total displacement on the shadown patch
        vectorField shadowPatchD
        (
            pointD.boundaryField()
            [
                shadowPatchIndices()[shadPatchI]
            ].patchInternalField()
        );

        // For a non-moving mesh, we will move the zones by the total
        // displacement, whereas for a moving mesh (updated Lagrangian), we will
        // move the zones by the displacement increment

        if (movingMesh())
        {
            // Subtract the old-time displacement
            shadowPatchD -=
                pointD.oldTime().boundaryField()
                [
                   shadowPatchIndices()[shadPatchI]
                ].patchInternalField();
        }

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField shadowZoneNewPoints
        (
            shadowZones()[shadPatchI].patchPointToGlobal
            (
                pointField
                (
                    shadowZones()[shadPatchI].patch().localPoints()
                  + shadowPatchD
                )
            )
        );

        // Remove zone weights
        shadowZones()[shadPatchI].movePoints(shadowZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        const_cast<pointField&>
        (
            shadowZones()[shadPatchI].globalPatch().points()
        ) = shadowZoneNewPoints;
    }
}


void Foam::solidContactPointPatchVectorField::calcZone() const
{
    if (debug)
    {
        InfoIn("void Foam::solidContactPointPatchVectorField::calcZone() const")
            << patch().name() << " : making the zone" << endl;
    }

    if (!master_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactPointPatchVectorField::calcZone() const"
        )   << "Trying to create zone on a slave" << abort(FatalError);
    }

    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactPointPatchVectorField::calcZone() const"
        )   << "pointer already set" << abort(FatalError);
    }

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ = new globalPolyPatch
    (
        patch().name(),
        lookupSolidModel(db()).mesh()
    );
}


void Foam::solidContactPointPatchVectorField::calcShadowZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidContactPointPatchVectorField::calcShadowZones() const"
        )   << patch().name() << " : making the shadow zones" << endl;
    }

    if (!master_)
    {
        FatalErrorIn
        (
            "void Foam::solidContactPointPatchVectorField::"
            "calcShadowZones() const"
        )   << "Trying to create shadow zones on a slave" << abort(FatalError);
    }

    if (!shadowZones_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidContactPointPatchVectorField::"
            "calcShadowZones() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const wordList& shadPatchNames = shadowPatchNames();

    shadowZones_.setSize(shadPatchNames.size());

    forAll(shadowZones_, shadPatchI)
    {
        // Note: the main mesh will either be in the initial configuration or
        // the updated configuration
        shadowZones_.set
        (
            shadPatchI,
            new globalPolyPatch
            (
                shadPatchNames[shadPatchI],
                lookupSolidModel(db()).mesh()
            )
        );
    }
}


void Foam::solidContactPointPatchVectorField::calcZoneToZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidContactPointPatchVectorField::calcZoneToZones() const"
        )   << patch().name() << " : making the zoneToZone" << endl;
    }

    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::calcZoneToZones() const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    // Check master and slave patch
    const pointVectorField& field =
        db().lookupObject<pointVectorField>
        (
#ifdef OPENFOAMESIORFOUNDATION
            internalField().name()
#else
            dimensionedInternalField().name()
#endif
        );

    zoneToZones_.setSize(shadowPatchNames().size());

    forAll (zoneToZones_, shadPatchI)
    {
        const solidContactPointPatchVectorField& shadowPatchField =
            refCast<const solidContactPointPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadPatchI]]
            );

        if (master_)
        {
            if (shadowPatchField.master() == true)
            {
                FatalErrorIn
                (
                    "void solidContactFvPatchScalarField::"
                    "calcZoneToZones() const"
                )   << "There are two master patches!" << abort(FatalError);
            }
        }
        else
        {
            if (shadowPatchField.master() == false)
            {
                FatalErrorIn
                (
                    "void solidContactFvPatchScalarField::"
                    "calcZoneToZones() const"
                )   << "There is no master patch!" << abort(FatalError);
            }
        }

        if (master_)
        {
            // Create interpolation for patches
            zoneToZones_.set
            (
                shadPatchI,
#ifdef OPENFOAMESIORFOUNDATION
                new amiZoneInterpolation
                (
                    zone().globalPatch(),
                    shadowZones()[shadPatchI].globalPatch(),
                    faceAreaIntersect::tmMesh, // triMode
                    true,   // requireMatch
                    -1,     // lowWeightCorrection
                    false,  // reverseTarget
                    true    // use globalPolyPatch
                )
#else
                new newGgiStandAlonePatchInterpolation
                (
                    zone().globalPatch(),
                    shadowZones()[shadPatchI].globalPatch(),
                    tensorField(0),
                    tensorField(0),
                    vectorField(0), // Slave-to-master separation. Bug fix
                    true,           // global data
                    0,              // Master non-overlapping face tolerances
                    0,              // Slave non-overlapping face tolerances
                    // Do not rescale weighting factors, as it is wrong on
                    // partially covered faces
                    false,
                    quickReject_,
                    regionOfInterest_
                )
#endif
            );

#ifndef OPENFOAMESIORFOUNDATION
            // Check which point distance calculation method to use
            const Switch useNewPointDistanceMethod =
                dict_.lookupOrDefault<Switch>
                (
                    "useNewPointDistanceMethod", false
                );

            Info<< "    " << type() << ": " << patch().name() << nl
                << "        useNewPointDistanceMethod: "
                << useNewPointDistanceMethod
                << endl;

            zoneToZones_[shadPatchI].useNewPointDistanceMethod() =
                useNewPointDistanceMethod;

            // Check if the projectPointsToPatchBoundary switch is set
            const Switch projectPointsToPatchBoundary =
                dict_.lookupOrDefault<Switch>
                (
                    "projectPointsToPatchBoundary",
                    false
                );

            Info<< "        projectPointsToPatchBoundary: "
                << projectPointsToPatchBoundary
                << endl;

            zoneToZones_[shadPatchI].projectPointsToPatchBoundary() =
                projectPointsToPatchBoundary;

            if (dict_.found("checkPointDistanceOrientations"))
            {
                const Switch checkPointDistanceOrientations =
                    Switch(dict_.lookup("checkPointDistanceOrientations"));

                Info<< "        checkPointDistanceOrientations: "
                    << checkPointDistanceOrientations
                    << endl;

                zoneToZones_[shadPatchI].checkPointDistanceOrientations() =
                    checkPointDistanceOrientations;
            }

            // Check if the usePrevCandidateMasterNeighbors switch is set
            const Switch usePrevCandidateMasterNeighbors =
                dict_.lookupOrDefault<Switch>
                (
                    "usePrevCandidateMasterNeighbors",
                    false
                );

            Info<< "        usePrevCandidateMasterNeighbors: "
                << usePrevCandidateMasterNeighbors
                << endl;

            zoneToZones_[shadPatchI].usePrevCandidateMasterNeighbors() =
                usePrevCandidateMasterNeighbors;
#endif
        }
        else
        {
            FatalErrorIn
            (
                "void solidContactPointPatchVectorField::"
                "calcZoneToZones() const"
            )   << "Attempting to create GGIInterpolation on a shadow patch"
                << abort(FatalError);
        }
    }
}


// void Foam::solidContactPointPatchVectorField::calcContactPerShadow() const
// {
//     if (contactPerShadow_.size() > 0)
//     {
//         FatalErrorIn
//         (
//             "void thermalContactFvPatchScalarField::"
//             "calcContactPerShadow() const"
//         )   << "already calculated"
//             << abort(FatalError);
//     }

//     contactPerShadow_.setSize(shadowPatchNames().size());

//     forAll(contactPerShadow_, i)
//     {
//         contactPerShadow_.set
//         (
//             i,
//             new scalarField(patch().size(), 0.0)
//         );
//     }
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalPolyPatch&
Foam::solidContactPointPatchVectorField::zone() const
{
    if (master_)
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }
    else
    {
        const pointVectorField& field =
            db().lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                internalField().name()
#else
                dimensionedInternalField().name()
#endif
            );

        const solidContactPointPatchVectorField& shadowPatchField =
            refCast<const solidContactPointPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.zone();
    }
}


Foam::globalPolyPatch& Foam::solidContactPointPatchVectorField::zone()
{
    if (master_)
    {
        if (!zonePtr_)
        {
            calcZone();
        }

        return *zonePtr_;
    }
    else
    {
        const pointVectorField& field =
            db().lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                internalField().name()
#else
                dimensionedInternalField().name()
#endif
            );

        solidContactPointPatchVectorField& shadowPatchField =
            const_cast<solidContactPointPatchVectorField&>
            (
                refCast<const solidContactPointPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.zone();
    }
}


const Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidContactPointPatchVectorField::shadowZones() const
{
    if (master_)
    {
        if (shadowZones_.empty())
        {
            calcShadowZones();
        }

        return shadowZones_;
    }
    else
    {
        const pointVectorField& field =
            db().lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                internalField().name()
#else
                dimensionedInternalField().name()
#endif
            );

        const solidContactPointPatchVectorField& shadowPatchField =
            refCast<const solidContactPointPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.shadowZones();
    }
}


Foam::PtrList<Foam::globalPolyPatch>&
Foam::solidContactPointPatchVectorField::shadowZones()
{
    if (master_)
    {
        if (shadowZones_.empty())
        {
            calcShadowZones();
        }

        return shadowZones_;
    }
    else
    {
        const pointVectorField& field =
            db().lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                internalField().name()
#else
                dimensionedInternalField().name()
#endif
            );

        solidContactPointPatchVectorField& shadowPatchField =
            const_cast<solidContactPointPatchVectorField&>
            (
                refCast<const solidContactPointPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.shadowZones();
    }
}


#ifdef OPENFOAMESIORFOUNDATION
const Foam::PtrList<Foam::amiZoneInterpolation>&
#else
const Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
#endif
Foam::solidContactPointPatchVectorField::zoneToZones() const
{
    if (master_)
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
    }
    else
    {
        const pointVectorField& field =
            db().lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                internalField().name()
#else
                dimensionedInternalField().name()
#endif
            );

        const solidContactPointPatchVectorField& shadowPatchField =
            refCast<const solidContactPointPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[0]]
            );

        return shadowPatchField.zoneToZones();
    }
}


#ifdef OPENFOAMESIORFOUNDATION
Foam::PtrList<Foam::amiZoneInterpolation>&
#else
Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
#endif
Foam::solidContactPointPatchVectorField::zoneToZones()
{
    if (master_)
    {
        if (zoneToZones_.empty())
        {
            calcZoneToZones();
        }

        return zoneToZones_;
    }
    else
    {
        // We will const_cast the shadow patch so we can delete the weights when
        // the zones move
        const pointVectorField& field =
            db().lookupObject<pointVectorField>
            (
#ifdef OPENFOAMESIORFOUNDATION
                internalField().name()
#else
                dimensionedInternalField().name()
#endif
            );

        solidContactPointPatchVectorField& shadowPatchField =
            const_cast<solidContactPointPatchVectorField&>
            (
                refCast<const solidContactPointPatchVectorField>
                (
                    field.boundaryField()[shadowPatchIndices()[0]]
                )
            );

        return shadowPatchField.zoneToZones();
    }

}


#ifdef OPENFOAMESIORFOUNDATION
const Foam::amiZoneInterpolation&
#else
const Foam::newGgiStandAlonePatchInterpolation&
#endif
Foam::solidContactPointPatchVectorField::zoneToZoneForThisSlave() const
{
    if (master_)
    {
        FatalErrorIn
        (
            "const Foam::newGgiStandAlonePatchInterpolation&"
            " Foam::solidContactPointPatchVectorField::"
            "zoneToZoneForThisSlave() const"
        )   << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which zoneToZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
    label masterShadowID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterShadowID = shadPatchI;
            break;
        }
    }

    if (masterShadowID == -1)
    {
        FatalErrorIn
        (
            "const Foam::newGgiStandAlonePatchInterpolation&"
            " Foam::solidContactPointPatchVectorField::"
            "zoneToZoneForThisSlave() const"
        )   << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    // Return the zoneToZone between the master and the current patch
    return zoneToZones()[masterShadowID];
}


const Foam::globalPolyPatch&
Foam::solidContactPointPatchVectorField::zoneForThisSlave() const
{
    if (master_)
    {
        FatalErrorIn
        (
            "const Foam::globalPolyPatch&"
            " Foam::solidContactPointPatchVectorField::zoneForThisSlave() const"
        )   << "The master patch is not allowed to call this function"
            << abort(FatalError);
    }

    // The master may have multiple slaves so we need to find which shadowZone
    // corresponds to the current slave patch
    const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
    label masterShadowID = -1;
    forAll(shadPatchNames, shadPatchI)
    {
        if (shadPatchNames[shadPatchI] == patch().name())
        {
            masterShadowID = shadPatchI;
            break;
        }
    }

    if (masterShadowID == -1)
    {
        FatalErrorIn
        (
            "const Foam::globalPolyPatch&"
            " Foam::solidContactPointPatchVectorField::zoneForThisSlave() const"
        )   << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    // Return the zoneToZone between the master and the current patch
    return shadowZones()[masterShadowID];
}


// ************************************************************************* //
