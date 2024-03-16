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
    solidRigidContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidRigidContactFvPatchVectorField.H"
#include "pointFields.H"
#include "PrimitivePatchInterpolationTemplate.H"
#include "triSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::solidRigidContactFvPatchVectorField::moveZonesToDeformedConfiguration()
{
    // Method
    // We will interpolate the patch face displacements to the patch vertices
    // and then add these vertex/point displacements to the initial patch
    // points
    // We need to take care in parallel, and also realise that the solidModel
    // might have a moving or stationary mesh

    // Assemble the zone face displacement field to move the zones

    // First we will move the master zone

    vectorField zoneD(zone().globalPatch().size(), vector::zero);

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

        zoneD = zone().patchFaceToGlobal(patchDD);
    }
    else
    {
        // Non-moving mesh: we will move the zones by the total displacement

        // Lookup the current total displacement field
        const volVectorField& D = db().lookupObject<volVectorField>("D");

        // Take a reference to the patch face total displacement field
        const vectorField& patchD =
            D.boundaryField()[patch().index()];

        zoneD = zone().patchFaceToGlobal(patchD);
    }

    // Interpolate the zone face field to the zone points
    // Note: for efficiency, interpolation should be stored.  However, you
    // already thorw it away every time, so it will be made and discarded
    // HJ, 23/Jun/2017
    const pointField zonePointD =
        PrimitivePatchInterpolation<standAlonePatch>
        (
            zone().globalPatch()
        ).faceToPointInterpolate(zoneD);

    // The zone deformed points are the initial position plus the
    // displacement
    const pointField zoneNewPoints =
        zone().patchPointToGlobal
        (
            patch().patch().localPoints()
        )
      + zonePointD;

    // HJ, Comment:
    // This is violence - you should not really deform the patch that mirrors
    // the mesh.  However, since there is no point in keeping another copy
    // of the zone, I will leave the code as is.
    // HJ, 23/Jun/2017

    // Remove zone weights
    zone().movePoints(zoneNewPoints);

    // We need to use const_cast to move the standAlonePatch points as the
    // movePoints function only clears weights
    // Also, be careful to move the points as opposed to the localPoints
    const_cast<pointField&>(zone().globalPatch().points()) = zoneNewPoints;
}


void Foam::solidRigidContactFvPatchVectorField::moveTriSurfaces()
{
    forAll(shadowTriSurfNames(), triSurfI)
    {
        const pointField newPoints = shadowZonesMotion()[triSurfI].newPoints();

        // Remove zone weights
        shadowZones()[triSurfI].movePoints(newPoints);

        // Move the triSurface points
        const_cast<pointField&>
        (
            shadowZones()[triSurfI].points()
        ) = newPoints;
    }
}


void Foam::solidRigidContactFvPatchVectorField::calcZone() const
{
    if (debug)
    {
        InfoIn("void Foam::solidRigidContactFvPatchVectorField::"
        "calcZone() const")
            << patch().name() << " : making the zone" << endl;
    }

    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::calcZone() const"
        )   << "pointer already set" << abort(FatalError);
    }

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ = new globalPolyPatch
    (
        patch().name(),
        patch().boundaryMesh().mesh()
    );
}


void Foam::solidRigidContactFvPatchVectorField::calcShadowZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::"
            "calcShadowZones() const"
        )   << patch().name() << " : making the shadow zones" << endl;
    }

    if (!shadowZones_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::"
            "calcShadowZones() const"
        )   << "pointer already set" << abort(FatalError);
    }

    shadowZones_.setSize(shadowTriSurfNames().size());

    forAll(shadowZones_, triSurfI)
    {
        // Read the STL
        triSurface stl
        (
            fileName
            (
                "./constant/triSurfaces"/shadowTriSurfNames()[triSurfI] + ".stl"
            )
        );

        // Convert STL local faces to a faceList
        List<labelledTri> stlTriFaces = stl.localFaces();
        faceList stlFaces(stlTriFaces.size());
        forAll(stlFaces, faceI)
        {
            stlFaces[faceI] = face(stlTriFaces[faceI]);
        }

        // Create standAlonePatch from the stl
        shadowZones_.set
        (
            triSurfI,
            new standAlonePatch
            (
                stlFaces, stl.localPoints()
            )
        );
    }
}


void Foam::solidRigidContactFvPatchVectorField::calcShadowZonesMotion() const
{
    if (!shadowZonesMotion_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::"
            "calcShadowZonesMotion() const"
        )   << "pointer already set" << abort(FatalError);
    }

    shadowZonesMotion_.setSize(shadowTriSurfNames().size());

    forAll(shadowZonesMotion_, triSurfI)
    {
        shadowZonesMotion_.set
        (
            triSurfI,
            new triSurfaceMotion
            (
                db().time(),
                dict_.subDict("triSurfaces").subDict
                (
                    shadowTriSurfNames()[triSurfI]
                ),
                shadowZones()[triSurfI]
            )
        );
    }
}


void Foam::solidRigidContactFvPatchVectorField::calcZoneToZones() const
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::"
            "calcZoneToZones() const"
        )   << patch().name() << " : making the zoneToZone" << endl;
    }

    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorIn
        (
            "void solidRigidContactFvPatchScalarField::calcZoneToZones() const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    zoneToZones_.setSize(shadowZones().size());

    forAll(zoneToZones_, triSurfI)
    {
        // Create interpolation for patches
        zoneToZones_.set
        (
            triSurfI,
            new newGgiStandAlonePatchInterpolation
            (
                zone().globalPatch(),
                shadowZones()[triSurfI],
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation
                true,           // global data
                0,              // Master non-overlapping face tolerances
                0,              // Slave non-overlapping face tolerances
                true,           // Rescale weighting factors
                false,
                quickReject_,
                regionOfInterest_
            )
        );

        // Check which point distance calculation method to use
        const Switch useNewPointDistanceMethod =
            dict_.lookupOrDefault<Switch>("useNewPointDistanceMethod", false);

        Info<< "    " << type() << ": " << patch().name() << nl
            << "        useNewPointDistanceMethod: "
            << useNewPointDistanceMethod
            << endl;

        zoneToZones_[triSurfI].useNewPointDistanceMethod() =
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

        zoneToZones_[triSurfI].projectPointsToPatchBoundary() =
            projectPointsToPatchBoundary;

        if (dict_.found("checkPointDistanceOrientations"))
        {
            const Switch checkPointDistanceOrientations =
                Switch(dict_.lookup("checkPointDistanceOrientations"));

            Info<< "        checkPointDistanceOrientations: "
                << checkPointDistanceOrientations
                << endl;

            zoneToZones_[triSurfI].checkPointDistanceOrientations() =
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

        zoneToZones_[triSurfI].usePrevCandidateMasterNeighbors() =
            usePrevCandidateMasterNeighbors;
    }
}


void Foam::solidRigidContactFvPatchVectorField::calcContactPerShadow() const
{
    if (contactPerShadow_.size() > 0)
    {
        FatalErrorIn
        (
            "void thermalContactFvPatchScalarField::"
            "calcContactPerShadow() const"
        )   << "already calculated"
            << abort(FatalError);
    }

    contactPerShadow_.setSize(shadowTriSurfNames().size());

    forAll(contactPerShadow_, i)
    {
        contactPerShadow_.set
        (
            i,
            new scalarField(patch().size(), 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::globalPolyPatch&
Foam::solidRigidContactFvPatchVectorField::zone() const
{
    if (!zonePtr_)
    {
        calcZone();
    }

    return *zonePtr_;
}


Foam::globalPolyPatch& Foam::solidRigidContactFvPatchVectorField::zone()
{
    if (!zonePtr_)
    {
        calcZone();
    }

    return *zonePtr_;
}


const Foam::PtrList<Foam::standAlonePatch>&
Foam::solidRigidContactFvPatchVectorField::shadowZones() const
{
    if (shadowZones_.empty())
    {
        calcShadowZones();
    }

    return shadowZones_;
}


Foam::PtrList<Foam::standAlonePatch>&
Foam::solidRigidContactFvPatchVectorField::shadowZones()
{
    if (shadowZones_.empty())
    {
        calcShadowZones();
    }

    return shadowZones_;
}


const Foam::PtrList<Foam::triSurfaceMotion>&
Foam::solidRigidContactFvPatchVectorField::shadowZonesMotion() const
{
    if (shadowZonesMotion_.empty())
    {
        calcShadowZonesMotion();
    }

    return shadowZonesMotion_;
}


Foam::PtrList<Foam::triSurfaceMotion>&
Foam::solidRigidContactFvPatchVectorField::shadowZonesMotion()
{
    if (shadowZonesMotion_.empty())
    {
        calcShadowZonesMotion();
    }

    return shadowZonesMotion_;
}


const Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solidRigidContactFvPatchVectorField::zoneToZones() const
{
    if (zoneToZones_.empty())
    {
        calcZoneToZones();
    }

    return zoneToZones_;
}


Foam::PtrList<Foam::newGgiStandAlonePatchInterpolation>&
Foam::solidRigidContactFvPatchVectorField::zoneToZones()
{
    if (zoneToZones_.empty())
    {
        calcZoneToZones();
    }

    return zoneToZones_;
}


// const Foam::newGgiStandAlonePatchInterpolation&
// Foam::solidRigidContactFvPatchVectorField::zoneToZoneForThisSlave() const
// {
//     if (master_)
//     {
//         FatalErrorIn
//         (
//             "const Foam::newGgiStandAlonePatchInterpolation&"
//             " Foam::solidRigidContactFvPatchVectorField::"
//             "zoneToZoneForThisSlave() const"
//         )   << "The master patch is not allowed to call this function"
//             << abort(FatalError);
//     }

//     // The master may have multiple slaves so we need to find which zoneToZone
//     // corresponds to the current slave patch
//     const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
//     label masterShadowID = -1;
//     forAll(shadPatchNames, shadPatchI)
//     {
//         if (shadPatchNames[shadPatchI] == patch().name())
//         {
//             masterShadowID = shadPatchI;
//             break;
//         }
//     }

//     if (masterShadowID == -1)
//     {
//         FatalErrorIn
//         (
//             "const Foam::newGgiStandAlonePatchInterpolation&"
//             " Foam::solidRigidContactFvPatchVectorField::"
//             "zoneToZoneForThisSlave() const"
//         )   << "Something went wrong when looking for the shadowPatch"
//             << abort(FatalError);
//     }

//     // Return the zoneToZone between the master and the current patch
//     return zoneToZones()[masterShadowID];
// }


// const Foam::globalPolyPatch&
// Foam::solidRigidContactFvPatchVectorField::zoneForThisSlave() const
// {
//     if (master_)
//     {
//         FatalErrorIn
//         (
//             "const Foam::globalPolyPatch&"
//             " Foam::solidRigidContactFvPatchVectorField::zoneForThisSlave() const"
//         )   << "The master patch is not allowed to call this function"
//             << abort(FatalError);
//     }

//     // The master may have multiple slaves so we need to find which shadowZone
//     // corresponds to the current slave patch
//     const wordList& shadPatchNames = shadowPatchField().shadowPatchNames();
//     label masterShadowID = -1;
//     forAll(shadPatchNames, shadPatchI)
//     {
//         if (shadPatchNames[shadPatchI] == patch().name())
//         {
//             masterShadowID = shadPatchI;
//             break;
//         }
//     }

//     if (masterShadowID == -1)
//     {
//         FatalErrorIn
//         (
//             "const Foam::globalPolyPatch&"
//             " Foam::solidRigidContactFvPatchVectorField::zoneForThisSlave() const"
//         )   << "Something went wrong when looking for the shadowPatch"
//             << abort(FatalError);
//     }

//     // Return the zoneToZone between the master and the current patch
//     return shadowZones()[masterShadowID];
// }


// ************************************************************************* //
