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
//#include "standardPenalty.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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

    // If the displacement increment field (DD) exists then we will assume that
    // the solidModel is using an incremental approach
    // TODO: we can add a function to solidModel called incremental() and
    // movingMesh()
    if (db().foundObject<volVectorField>("DD"))
    {
        // Linear geometry incremental models
        FatalErrorIn
        (
            "solidContactFvPatchVectorField::moveZonesToDeformedConfiguration()"
        ) << "wip: on my to-do list" << abort(FatalError);
    }
    else if (db().foundObject<volVectorField>("F"))
    {
        // Nonlinear geometry models
        FatalErrorIn
        (
            "solidContactFvPatchVectorField::moveZonesToDeformedConfiguration()"
        ) << "wip: on my to-do list" << abort(FatalError);
    }
    else
    {
        // Lookup the current total displacement field
        const volVectorField& D = db().lookupObject<volVectorField>("D");

        // Take a reference to the patch face total displacement field
        const vectorField& patchD = D.boundaryField()[patch().index()];
        const vectorField& shadowPatchD = D.boundaryField()[shadowPatchIndex()];

        // Assemble the zone face total displacement field
        const vectorField zoneD =
            zoneField(zoneIndex(), patch().index(), patchD);
        const vectorField shadowZoneD =
            zoneField(shadowZoneIndex(), shadowPatchIndex(), shadowPatchD);

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
        zone().movePoints(zoneNewPoints);
        shadowZone().movePoints(shadowZoneNewPoints);
    }
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
            patch().index(),        // master
            shadowPatchIndex(),     // slave
            zoneIndex(),            // master face zone ID
            shadowZoneIndex(),      // slave face zone ID
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
            patch().index(), // master
            shadowPatchIndex(), // slave
            zoneIndex(), // master face zone ID
            shadowZoneIndex() // slave face zone ID
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


// void Foam::solidContactFvPatchVectorField::calcAllPointsDeformed() const
// {
//     if (allPointsDeformedPtr_)
//     {
//         FatalErrorIn
//         (
//             "void Foam::solidContactFvPatchVectorField::"
//             "calcAllPointsDeformed() const"
//         )   << "allPointsDeformedPtr_ already set" << abort(FatalError);
//     }

//     const fvMesh& mesh = patch().boundaryMesh().mesh();

//     allPointsDeformedPtr_ = new vectorField(mesh.allPoints());
// }


// Foam::vectorField&
// Foam::solidContactFvPatchVectorField::allPointsDeformed()
// {
//     if (master())
//     {
//         if (!allPointsDeformedPtr_)
//         {
//             calcAllPointsDeformed();
//         }

//         return *allPointsDeformedPtr_;
//     }

//     // slave
//     const volVectorField& field =
//         db().lookupObject<volVectorField>
//         (
//             this->dimensionedInternalField().name()
//         );

//     solidContactFvPatchVectorField& shadowPatchField =
//         const_cast<solidContactFvPatchVectorField&>
//         (
//             refCast<const solidContactFvPatchVectorField>
//             (
//                 field.boundaryField()[shadowPatchIndex_]
//             )
//         );

//     return shadowPatchField.allPointsDeformed();
// }


// const Foam::vectorField&
// Foam::solidContactFvPatchVectorField::allPointsDeformed() const
// {
//     if (master())
//     {
//         if (!allPointsDeformedPtr_)
//         {
//             calcAllPointsDeformed();
//         }

//         return *allPointsDeformedPtr_;
//     }

//     // slave
//     const volVectorField& field =
//         db().lookupObject<volVectorField>
//         (
//             this->dimensionedInternalField().name()
//         );

//     const solidContactFvPatchVectorField& shadowPatchField =
//         refCast<const solidContactFvPatchVectorField>
//         (
//             field.boundaryField()[shadowPatchIndex_]
//         );

//     return shadowPatchField.allPointsDeformed();
// }


// void Foam::solidContactFvPatchVectorField::updateAllPointsDeformed()
// {
//     if (master())
//     {
//         vectorField& p = allPointsDeformed();

//         // lookup displacement field
//         const pointVectorField& field =
//             db().lookupObject<pointVectorField>
//             (
//                 "point" + this->dimensionedInternalField().name()
//             );

//         // If pointDU is not updated in the outer iterations of the solver then
//         // we must update it here
//         if (!solverUpdatesPointDU_)
//         {
//             constCastUpdatePointDU();
//         }

//         // Total field is used for non-moving mesh incremental solvers
//         // We will check if the deformation gradient F is defined
//         bool incrementalNonMovingMesh = false;
//         if
//         (
//             fieldName_ == "DD"
//          && !db().lookupObject<volTensorField>("F")
//         )
//         {
//             incrementalNonMovingMesh = true;
//         }

//         const pointVectorField& totalField =
//             db().lookupObject<pointVectorField>
//             (
//                 "pointD"
//             );

//         // pointDU does not include unused points from the global face zones
//         const label nUsedPoints = field.size();

//         // Move zone points

//         const vectorField& oldZoneP = oldZonePoints();
//         const vectorField& oldShadowZoneP = oldShadowZonePoints();

//         // Zone may not have been created so we will use the zones in the mesh
//         //const labelList& zoneMeshPoints = zone().meshPoints();
//         //const labelList& shadowZoneMeshPoints = shadowZone().meshPoints();
//         const fvMesh& mesh = patch().boundaryMesh().mesh();
//         const labelList& zoneMeshPoints =
//             mesh.faceZones()[zoneIndex()]().meshPoints();
//         const labelList& shadowZoneMeshPoints =
//             mesh.faceZones()[shadowZoneIndex()]().meshPoints();

//         forAll(zoneMeshPoints, pI)
//         {
//             const label pointID = zoneMeshPoints[pI];

//             p[pointID] = oldZoneP[pI];

//             if (pointID < nUsedPoints)
//             {
//                 p[pointID] += field[pointID];

//                 if (incrementalNonMovingMesh)
//                 {
//                     p[pointID] += totalField[pointID];
//                 }
//             }
//         }

//         forAll(shadowZoneMeshPoints, pI)
//         {
//             const label pointID = shadowZoneMeshPoints[pI];

//             p[pointID] = oldShadowZoneP[pI];

//             if (pointID < nUsedPoints)
//             {
//                 p[pointID] += field[pointID];

//                 if (incrementalNonMovingMesh)
//                 {
//                     p[pointID] += totalField[pointID];
//                 }
//             }
//         }

//         // Set values for unused points

//         const vectorField& pointFieldI = field.internalField();
//         const vectorField& pointTotalFieldI = field.internalField();

//         const labelList& gFaceZones = globalFaceZones();

//         forAll(gFaceZones, zoneI)
//         {
//             const label curZoneID = gFaceZones[zoneI];

//             const labelList& curMap =
//                 globalToLocalFaceZonePointMap()[zoneI];

//             const labelList& curZoneMeshPoints =
//                 mesh.faceZones()[curZoneID]().meshPoints();

//             vectorField curGlobalZonePointDispl
//                 (
//                     curZoneMeshPoints.size(),
//                     vector::zero
//                 );

//             //-Inter-proc points are shared by multiple procs
//             // pointNumProc is the number of procs which a point lies on
//             scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

//             forAll(curGlobalZonePointDispl, globalPointI)
//             {
//                 label localPoint = curMap[globalPointI];

//                 if(curZoneMeshPoints[localPoint] < mesh.nPoints())
//                 {
//                     label procPoint = curZoneMeshPoints[localPoint];

//                     curGlobalZonePointDispl[globalPointI] =
//                         pointFieldI[procPoint];

//                     if (incrementalNonMovingMesh)
//                     {
//                         curGlobalZonePointDispl[globalPointI] +=
//                             pointTotalFieldI[procPoint];
//                     }

//                     pointNumProcs[globalPointI] = 1;
//                 }
//             }

//             if (Pstream::parRun())
//             {
//                 reduce(curGlobalZonePointDispl, sumOp<vectorField>());
//                 reduce(pointNumProcs, sumOp<scalarField>());

//                 // nNow average the displacement between all procs
//                 curGlobalZonePointDispl /= pointNumProcs;
//             }

//             // The curZonePointsDisplGlobal now contains the correct face zone
//             // displacement in a global master processor order, now convert
//             // them back into the local proc order

//             vectorField curZonePointDispl
//             (
//                 curZoneMeshPoints.size(),
//                 vector::zero
//             );

//             forAll(curGlobalZonePointDispl, globalPointI)
//             {
//                 label localPoint = curMap[globalPointI];

//                 curZonePointDispl[localPoint] =
//                     curGlobalZonePointDispl[globalPointI];
//             }

//             forAll(curZonePointDispl, pointI)
//             {
//                 // Unused points
//                 if (curZoneMeshPoints[pointI] >= mesh.nPoints())
//                 {
//                     // Pout<< "correcting motion for point "
//                     //     << curZoneMeshPoints[pointI] << endl;
//                     p[curZoneMeshPoints[pointI]] += curZonePointDispl[pointI];
//                 }
//             }
//         }
//     }
//     else
//     {
//         const volVectorField& field =
//             db().objectRegistry::lookupObject<volVectorField>
//             (
//                 this->dimensionedInternalField().name()
//             );

//         // We will const cast the shadow patch to update the master points

//         solidContactFvPatchVectorField& shadowPatchField =
//             const_cast<solidContactFvPatchVectorField&>
//             (
//                 refCast<const solidContactFvPatchVectorField>
//                 (
//                     field.boundaryField()[shadowPatchIndex_]
//                 )
//             );

//         return shadowPatchField.updateAllPointsDeformed();
//     }
// }


// void Foam::solidContactFvPatchVectorField::constCastUpdatePointDU() const
// {
//     if
//     (
//         fieldName_ != "DU"
//      || nonLinear_ != nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
//     )
//     {
//         FatalErrorIn("solidContactFvPatchVectorField::constCastUpdatePointDU")
//             << "Only implemented for DU and updateLagrangianKirchhoff"
//             << abort(FatalError);
//     }

//     // Method taken from old function
//     // solidContactFvPatchVectorField::moveFaceZonePatches()

//     // Method: we get the total displacement field for the global
//     // face zone patches. We then interpolate these face values
//     // to the vertices. And we move the vertices by these
//     // interpolated displacements, so the global face zone patches
//     // should be in the same deformed position on all procs.
//     // Note: we use const_cast here to update pointDU

//     // Update face zone patch interpolators
//     //masterFaceZonePatchInterpolatorPtr_->movePoints();
//     //slaveFaceZonePatchInterpolatorPtr_->movePoints();

//     // Lookup displacement field and const_cast as we will update it
//     pointVectorField& pointDU =
//         const_cast<pointVectorField&>
//         (
//             db().lookupObject<pointVectorField>
//             (
//                 "point" + this->dimensionedInternalField().name()
//             )
//         );

//     const fvMesh& mesh = patch().boundaryMesh().mesh();

//     // Get local total displacement fields
//     const volVectorField& dispField =
//         db().lookupObject<volVectorField>("DU");

//     const vectorField globalMasterU =
//         zoneField
//         (
//             zoneIndex(),
//             patch().index(),
//             dispField.boundaryField()[patch().index()]
//         );

//     const vectorField globalSlaveU =
//         zoneField
//         (
//             shadowZoneIndex(),
//             shadowPatchIndex(),
//             dispField.boundaryField()[shadowPatchIndex()]
//         );

//     // Interpolate displacement from face centre to vertices

//     vectorField globalMasterPointDU =
//         zoneFaceToPointInterpolate
//         (
//             zoneIndex(), globalMasterU
//         );

//     vectorField globalSlavePointDU =
//         zoneFaceToPointInterpolate
//         (
//             shadowZoneIndex(), globalSlaveU
//         );

//     // Remove displacement is empty direction for 2-D
//     if (mesh.nSolutionD() == 2)
//     {
//         const Vector<label>& solD = mesh.solutionD();

//         forAll(solD, dirI)
//         {
//             if (solD[dirI] == -1)
//             {
//                 forAll(globalMasterPointDU, pointI)
//                 {
//                     globalMasterPointDU[pointI][dirI] = 0.0;
//                 }

//                 forAll(globalSlavePointDU, pointI)
//                 {
//                     globalSlavePointDU[pointI][dirI] = 0.0;
//                 }
//             }
//         }
//     }

//     // Overwrite pointDU boundary points using const_cast

//     const label nLocalPoints = mesh.nPoints();

//     {
//         const labelList& mp = zone().meshPoints();
//         forAll(mp, pI)
//         {
//             const label pointID = mp[pI];

//             if (pointID < nLocalPoints)
//             {
//                 pointDU[pointID] = globalMasterPointDU[pI];
//             }
//         }
//     }
//     {
//         const labelList& mp = shadowZone().meshPoints();
//         forAll(mp, pI)
//         {
//             const label pointID = mp[pI];

//             if (pointID < nLocalPoints)
//             {
//                 pointDU[pointID] = globalSlavePointDU[pI];
//             }
//         }
//     }
// }


// const Foam::vectorField&
// Foam::solidContactFvPatchVectorField::oldZonePoints() const
// {
//     if (master())
//     {
//         if (!oldZonePointsPtr_)
//         {
//             calcOldZonePoints();
//         }

//         return *oldZonePointsPtr_;
//     }

//     // slave

//     const volVectorField& field =
//         db().lookupObject<volVectorField>
//         (
//             this->dimensionedInternalField().name()
//         );

//     const solidContactFvPatchVectorField& shadowPatchField =
//         refCast<const solidContactFvPatchVectorField>
//         (
//             field.boundaryField()[shadowPatchIndex_]
//         );

//     return shadowPatchField.oldZonePoints();
// }


// void Foam::solidContactFvPatchVectorField::calcOldZonePoints() const
// {
//     if (oldZonePointsPtr_)
//     {
//         FatalErrorIn
//         (
//             "void Foam::solidContactFvPatchVectorField::"
//             "calcOldZonePoints() const"
//         )   << "pointer already set"
//             << abort(FatalError);
//     }

//     if
//     (
//         nonLinear_ != nonLinearGeometry::OFF
//         && nonLinear_ != nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
//     )
//     {
//         FatalErrorIn
//         (
//             "void vectorField& Foam::solidContactFvPatchVectorField::"
//             "oldZonePoints() const"
//         )   << "only implemented for nonLinear updatedLagrangianKirchhoff"
//             << abort(FatalError);
//     }

//     const fvMesh& mesh = patch().boundaryMesh().mesh();

//     oldZonePointsPtr_ =
//         new vectorField(mesh.faceZones()[zoneIndex()]().localPoints());
// }


// const Foam::vectorField&
// Foam::solidContactFvPatchVectorField::oldShadowZonePoints() const
// {
//     if (master())
//     {
//         if (!oldShadowZonePointsPtr_)
//         {
//             calcOldShadowZonePoints();
//         }

//         return *oldShadowZonePointsPtr_;
//     }

//     // slave

//     const volVectorField& field =
//         db().lookupObject<volVectorField>
//         (
//             this->dimensionedInternalField().name()
//         );

//     const solidContactFvPatchVectorField& shadowPatchField =
//         refCast<const solidContactFvPatchVectorField>
//         (
//             field.boundaryField()[shadowPatchIndex_]
//         );

//     return shadowPatchField.oldZonePoints();
// }


// void Foam::solidContactFvPatchVectorField::calcOldShadowZonePoints() const
// {
//     if (oldShadowZonePointsPtr_)
//     {
//         FatalErrorIn
//         (
//             "void Foam::solidContactFvPatchVectorField::"
//             "calcOldShadowZonePoints() const"
//         )   << "pointer already set"
//             << abort(FatalError);
//     }

//     if
//     (
//         nonLinear_ != nonLinearGeometry::OFF
//         && nonLinear_ != nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
//     )
//     {
//         FatalErrorIn
//         (
//             "void vectorField& Foam::solidContactFvPatchVectorField::"
//             "oldShadowZonePoints() const"
//         )   << "only implemented for nonLinear updatedLagrangianKirchhoff"
//             << abort(FatalError);
//     }

//     const fvMesh& mesh = patch().boundaryMesh().mesh();

//     oldShadowZonePointsPtr_ =
//         new vectorField(mesh.faceZones()[shadowZoneIndex()]().localPoints());
// }


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
        new primitiveFacePatch
        (
            // faceList(mesh.faceZones()[zoneIndex()].size()),
            // allPointsDeformed()
            mesh.faceZones()[zoneIndex()]().localFaces(),
            mesh.faceZones()[zoneIndex()]().localPoints()
        );

    // primitiveFacePatch& patch = *zonePtr_;

    // const faceList& f = mesh.allFaces();

    // const labelList& addr = mesh.faceZones()[zoneIndex()];
    // const boolList& flip = mesh.faceZones()[zoneIndex()].flipMap();

    // forAll (addr, faceI)
    // {
    //     if (flip[faceI])
    //     {
    //         patch[faceI] = f[addr[faceI]].reverseFace();
    //     }
    //     else
    //     {
    //         patch[faceI] = f[addr[faceI]];
    //     }
    // }
}


const Foam::primitiveFacePatch&
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


Foam::primitiveFacePatch& Foam::solidContactFvPatchVectorField::zone()
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
        new primitiveFacePatch
        (
            // faceList(mesh.faceZones()[zoneIndex()].size()),
            // allPointsDeformed()
            mesh.faceZones()[shadowZoneIndex()]().localFaces(),
            mesh.faceZones()[shadowZoneIndex()]().localPoints()
        );

    // primitiveFacePatch& patch = *shadowZonePtr_;

    // const faceList& f = mesh.allFaces();

    // const labelList& addr = mesh.faceZones()[shadowZoneIndex()];
    // const boolList& flip = mesh.faceZones()[shadowZoneIndex()].flipMap();

    // forAll (addr, faceI)
    // {
    //     if (flip[faceI])
    //     {
    //         patch[faceI] = f[addr[faceI]].reverseFace();
    //     }
    //     else
    //     {
    //         patch[faceI] = f[addr[faceI]];
    //     }
    // }
}


const Foam::primitiveFacePatch&
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


Foam::primitiveFacePatch& Foam::solidContactFvPatchVectorField::shadowZone()
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
            new extendedGgiZoneInterpolation
            (
                zone(),
                shadowZone(),
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation. Bug fix
                0,              // Master non-overlapping face tolerances
                0,              // Slave non-overlapping face tolerances
                true,           // Rescale weighting factors
                quickReject_,
                regionOfInterest_
            );

        Info<< "Region of interest: " << regionOfInterest_ << endl;

        // Testing: increase BB span
        //debug::tolerances().set("GGIFaceBoundBoxExtendSpanFraction", 10);
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


const Foam::extendedGgiZoneInterpolation&
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


Foam::extendedGgiZoneInterpolation&
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


// void Foam::solidContactFvPatchVectorField::calcSlaveFaceNormals() const
// {
//     // Create slave face normals
//     if (slaveFaceNormalsPtr_)
//     {
//         FatalErrorIn
//         (
//             "void solidContactFvPatchScalarField::calcSlaveFaceNormals() const"
//         )   << "slaveFaceNormals pointer already calculated"
//             << abort(FatalError);
//     }

//     slaveFaceNormalsPtr_ =
//         new vectorField
//         (
//             patch().boundaryMesh().mesh().boundary()[shadowPatchIndex()].size(),
//             vector::zero
//         );
// }


// Foam::vectorField&
// Foam::solidContactFvPatchVectorField::slaveFaceNormals()
// {
//     if (master())
//     {
//         if (!slaveFaceNormalsPtr_)
//         {
//             calcSlaveFaceNormals();
//         }

//         return *slaveFaceNormalsPtr_;
//     }

//     // Slave

//     const volVectorField& field =
//         db().lookupObject<volVectorField>
//         (
//             this->dimensionedInternalField().name()
//         );

//  // We will const_cast the shadow patch so we can delete the weights when the
//     // zones move
//     solidContactFvPatchVectorField& shadowPatchField =
//         const_cast<solidContactFvPatchVectorField&>
//         (
//             refCast<const solidContactFvPatchVectorField>
//             (
//                 field.boundaryField()[shadowPatchIndex_]
//             )
//         );

//     return shadowPatchField.slaveFaceNormals();
// }


void Foam::solidContactFvPatchVectorField::calcGlobalFaceZones() const
{
    // Create slave face normals
    if (globalFaceZonesPtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::"
            "calcGlobalFaceZones() const"
        )   << "pointer already calculated"
            << abort(FatalError);
    }

    globalFaceZonesPtr_ = new labelList(0);

    labelList& globalFaceZones = *globalFaceZonesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (Pstream::parRun())
    {
        SLList<label> globalFaceZonesSet;

        // Lookup globalFaceZones from decomposeParDict
        // To be fixed for FSI cases
        IOdictionary decompDict
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().time().system(),
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (decompDict.found("globalFaceZones"))
        {
            wordList globalFaceZoneNames(decompDict.lookup("globalFaceZones"));

            const faceZoneMesh& faceZones = mesh.faceZones();

            forAll(globalFaceZoneNames, nameI)
            {
                const label zoneID =
                    faceZones.findZoneID(globalFaceZoneNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn
                    (
                        "solidContactFvPatchVectorField::"
                        "solidContactFvPatchVectorField"
                    )   << "cannot find globalFaceZone:"
                        << " " << globalFaceZoneNames[nameI]
                        << abort(FatalError);
                }

                globalFaceZonesSet.insert(zoneID);
            }
        }

        globalFaceZones = labelList(globalFaceZonesSet);
    }
}


const Foam::labelList&
Foam::solidContactFvPatchVectorField::globalFaceZones() const
{
    if (master())
    {
        if (!globalFaceZonesPtr_)
        {
            calcGlobalFaceZones();
        }

        return *globalFaceZonesPtr_;
    }

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.globalFaceZones();
}


void Foam::solidContactFvPatchVectorField::
calcGlobalToLocalFaceZonePointMap() const
{
    // Create slave face normals
    if (globalToLocalFaceZonePointMapPtr_)
    {
        FatalErrorIn
        (
            "void solidContactFvPatchScalarField::"
            "calcGlobalToLocalFaceZonePointMap() const"
        )   << "slaveFaceNormals pointer already calculated"
            << abort(FatalError);
    }

    const labelList& gFaceZones = globalFaceZones();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    fileName mapName = "globalToLocalFaceZonePointMapSolidContact";

    word timeName = mesh.time().timeName();

    bool mapExists = false;

    {
        IOobject mapHeader
        (
            mapName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (mapHeader.headerOk())
        {
            mapExists = true;
        }
        else
        {
            // Check previous time-step
            instantList times = mesh.time().times();

            if (times.size() > 1)
            {
                word prevTimeName = times[times.size() - 1].name();

                Pout<< "Reading face map from " << prevTimeName << endl;

                IOobject mapPrevHeader
                (
                    mapName,
                    prevTimeName,
                    mesh,
                    IOobject::MUST_READ
                );

                if (mapPrevHeader.headerOk())
                {
                    mapExists = true;
                    timeName = prevTimeName;
                }
            }
        }
    }

    globalToLocalFaceZonePointMapPtr_ =
        new IOList<labelList>
        (
            IOobject
            (
                mapName,
                timeName,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            gFaceZones.size()
        );

    IOList<labelList>& globalToLocalFaceZonePointMap =
        *globalToLocalFaceZonePointMapPtr_;

    // Calculate map if it has not been read
    if (!mapExists && Pstream::parRun())
    {
        forAll(gFaceZones, zoneI)
        {
            label curZoneID = gFaceZones[zoneI];

            Info<< "Creating faceMap for globalFaceZones "
                << mesh.faceZones()[curZoneID].name()<< endl;

            labelList curMap(mesh.faceZones()[curZoneID]().nPoints(), -1);

            vectorField fzGlobalPoints =
                mesh.faceZones()[curZoneID]().localPoints();

            // Set all slave points to zero because only the master order is
            // used
            if (!Pstream::master())
            {
                fzGlobalPoints = vector::zero;
            }

            //- pass points to all procs
            reduce(fzGlobalPoints, sumOp<vectorField>());

            // Now every proc has the master's list of FZ points
            // every proc must now find the mapping from their local FZ points
            // to the global FZ points

            const vectorField& fzLocalPoints =
                mesh.faceZones()[curZoneID]().localPoints();

            const edgeList& fzLocalEdges =
                mesh.faceZones()[curZoneID]().edges();

            const labelListList& fzPointEdges =
                mesh.faceZones()[curZoneID]().pointEdges();

            scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

            forAll(minEdgeLength, pI)
            {
                const labelList& curPointEdges = fzPointEdges[pI];

                forAll(curPointEdges, eI)
                {
                    scalar Le =
                        fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
                    if (Le < minEdgeLength[pI])
                    {
                        minEdgeLength[pI] = Le;
                    }
                }
            }

            forAll(fzGlobalPoints, globalPointI)
            {
                //scalar minDist = GREAT;
                bool pointFound = false;

                forAll(fzLocalPoints, procPointI)
                {
                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[procPointI]
                            - fzGlobalPoints[globalPointI]
                        );

                    if (curDist < 1e-4*minEdgeLength[procPointI])
                    {
                        curMap[globalPointI] = procPointI;

                        if (pointFound)
                        {
                            FatalError
                                << "findGlobalFaceZones: point found twice!"
                                << abort(FatalError);
                        }
                        pointFound = true;
                    }
                }
            }

            forAll(curMap, globalPointI)
            {
                if (curMap[globalPointI] == -1)
                {
                    FatalErrorIn
                    (
                        "solidContactFvPatchVectorField::"
                        "calcGlobalToLocalFaceZonePointMap()"
                    )   << "local to global face zone point map is not correct"
                        << " for zone " << zoneI
                        << abort(FatalError);
                }
            }

            globalToLocalFaceZonePointMap[zoneI] = curMap;
        }
    }
    else
    {
        Info<< "globalToLocalFaceZonePointMap read from file" << endl;
    }
}


const Foam::IOList<Foam::labelList>&
Foam::solidContactFvPatchVectorField::globalToLocalFaceZonePointMap() const
{
    if (master())
    {
        if (!globalToLocalFaceZonePointMapPtr_)
        {
            calcGlobalToLocalFaceZonePointMap();
        }

        return *globalToLocalFaceZonePointMapPtr_;
    }

    // Slave

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            this->dimensionedInternalField().name()
        );

    // We will const_cast the shadow patch so we can delete the weights when the
    // zones move
    const solidContactFvPatchVectorField& shadowPatchField =
        refCast<const solidContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndex_]
        );

    return shadowPatchField.globalToLocalFaceZonePointMap();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    master_(false),
    shadowPatchName_("undefined"),
    shadowPatchIndex_(-1),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(false),
    dict_(NULL),
    // normalModelPtr_(NULL),
    // frictionModelPtr_(NULL),
    //allPointsDeformedPtr_(NULL),
    // oldZonePointsPtr_(NULL),
    // oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_(Foam::extendedGgiInterpolation::AABB),
    regionOfInterest_(vector::min, vector::max),
    curTimeIndex_(-1),
    // slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL)
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
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    // normalModelPtr_(NULL),
    // frictionModelPtr_(NULL),
    //allPointsDeformedPtr_(NULL),
    // oldZonePointsPtr_(NULL),
    // oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_(ptf.quickReject_),
    regionOfInterest_(ptf.regionOfInterest_),
    curTimeIndex_(ptf.curTimeIndex_),
    // slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL)
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

    // if (ptf.normalModelPtr_)
    // {
    //     normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    // }

    // if (ptf.frictionModelPtr_)
    // {
    //     frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    // }

    // if (ptf.allPointsDeformedPtr_)
    // {
    //     allPointsDeformedPtr_ = new vectorField(*ptf.allPointsDeformedPtr_);
    // }

    // if (ptf.oldZonePointsPtr_)
    // {
    //     oldZonePointsPtr_ = new pointField(*ptf.oldZonePointsPtr_);
    // }

    // if (ptf.oldShadowZonePointsPtr_)
    // {
    //     oldShadowZonePointsPtr_ = new pointField(*ptf.oldShadowZonePointsPtr_);
    // }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new primitiveFacePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new primitiveFacePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed

    // if (ptf.slaveFaceNormalsPtr_)
    // {
    //     slaveFaceNormalsPtr_ = new vectorField(*ptf.slaveFaceNormalsPtr_);
    // }

    if (ptf.globalFaceZonesPtr_)
    {
        globalFaceZonesPtr_ = new labelList(*ptf.globalFaceZonesPtr_);
    }

    if (ptf.globalToLocalFaceZonePointMapPtr_)
    {
        globalToLocalFaceZonePointMapPtr_ =
            new IOList<labelList>(*ptf.globalToLocalFaceZonePointMapPtr_);
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    master_(dict.lookup("master")),
    shadowPatchName_(dict.lookup("shadowPatch")),
    shadowPatchIndex_(-1),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(false),
    dict_(dict),
    // normalModelPtr_(NULL),
    // frictionModelPtr_(NULL),
    //allPointsDeformedPtr_(NULL),
    // oldZonePointsPtr_(NULL),
    // oldShadowZonePointsPtr_(NULL),
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
    regionOfInterest_
    (
        dict.lookupOrDefault<boundBox>
        (
            "regionOfInterest",
            boundBox(vector::min, vector::max)
        )
    ),
    curTimeIndex_(-1),
    //slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL)
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
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    // normalModelPtr_(ptf.normalModelPtr_),
    // frictionModelPtr_(ptf.frictionModelPtr_),
    //allPointsDeformedPtr_(ptf.allPointsDeformedPtr_),
    // oldZonePointsPtr_(ptf.oldZonePointsPtr_),
    // oldShadowZonePointsPtr_(ptf.oldShadowZonePointsPtr_),
    zonePtr_(ptf.zonePtr_),
    shadowZonePtr_(ptf.shadowZonePtr_),
    zoneToZonePtr_(ptf.zoneToZonePtr_),
    quickReject_(ptf.quickReject_),
    regionOfInterest_(ptf.regionOfInterest_),
    curTimeIndex_(ptf.curTimeIndex_),
    //slaveFaceNormalsPtr_(ptf.slaveFaceNormalsPtr_),
    globalFaceZonesPtr_(ptf.globalFaceZonesPtr_),
    globalToLocalFaceZonePointMapPtr_(ptf.globalToLocalFaceZonePointMapPtr_)
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

    // if (ptf.normalModelPtr_)
    // {
    //     normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    // }

    // if (ptf.frictionModelPtr_)
    // {
    //     frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    // }

    // // if (ptf.allPointsDeformedPtr_)
    // // {
    // //     allPointsDeformedPtr_ = new vectorField(*ptf.allPointsDeformedPtr_);
    // // }

    // if (ptf.oldZonePointsPtr_)
    // {
    //     oldZonePointsPtr_ = new pointField(*ptf.oldZonePointsPtr_);
    // }

    // if (ptf.oldShadowZonePointsPtr_)
    // {
    //     oldShadowZonePointsPtr_ = new pointField(*ptf.oldShadowZonePointsPtr_);
    // }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new primitiveFacePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new primitiveFacePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed

    // if (ptf.slaveFaceNormalsPtr_)
    // {
    //     slaveFaceNormalsPtr_ = new vectorField(*ptf.slaveFaceNormalsPtr_);
    // }

    if (ptf.globalFaceZonesPtr_)
    {
        globalFaceZonesPtr_ = new labelList(*ptf.globalFaceZonesPtr_);
    }

    if (ptf.globalToLocalFaceZonePointMapPtr_)
    {
        globalToLocalFaceZonePointMapPtr_ =
            new IOList<labelList>(*ptf.globalToLocalFaceZonePointMapPtr_);
    }
}


Foam::solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndexPtr_(NULL),
    shadowZoneIndexPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    //normalModelPtr_(ptf.normalModelPtr_),
    //frictionModelPtr_(ptf.frictionModelPtr_),
    //allPointsDeformedPtr_(NULL),
    //oldZonePointsPtr_(NULL),
    //oldShadowZonePointsPtr_(NULL),
    zonePtr_(NULL),
    shadowZonePtr_(NULL),
    zoneToZonePtr_(NULL),
    quickReject_(ptf.quickReject_),
    regionOfInterest_(ptf.regionOfInterest_),
    curTimeIndex_(ptf.curTimeIndex_),
    //slaveFaceNormalsPtr_(NULL),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL)
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

    // if (ptf.normalModelPtr_)
    // {
    //     normalModelPtr_ = ptf.normalModelPtr_->clone().ptr();
    // }

    // if (ptf.frictionModelPtr_)
    // {
    //     frictionModelPtr_ = ptf.frictionModelPtr_->clone().ptr();
    // }

    // if (ptf.allPointsDeformedPtr_)
    // {
    //     allPointsDeformedPtr_ = new vectorField(*ptf.allPointsDeformedPtr_);
    // }

    // if (ptf.oldZonePointsPtr_)
    // {
    //     oldZonePointsPtr_ = new pointField(*ptf.oldZonePointsPtr_);
    // }

    // if (ptf.oldShadowZonePointsPtr_)
    // {
    //     oldShadowZonePointsPtr_ = new pointField(*ptf.oldShadowZonePointsPtr_);
    // }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new primitiveFacePatch(*ptf.zonePtr_);
    }

    if (ptf.shadowZonePtr_)
    {
        shadowZonePtr_ = new primitiveFacePatch(*ptf.shadowZonePtr_);
    }

    // We will not copy zoneToZonePtr_; it will have to be re-created when
    // needed

    // if (ptf.slaveFaceNormalsPtr_)
    // {
    //     slaveFaceNormalsPtr_ = new vectorField(*ptf.slaveFaceNormalsPtr_);
    // }

    if (ptf.globalFaceZonesPtr_)
    {
        globalFaceZonesPtr_ = new labelList(*ptf.globalFaceZonesPtr_);
    }

    if (ptf.globalToLocalFaceZonePointMapPtr_)
    {
        globalToLocalFaceZonePointMapPtr_ =
            new IOList<labelList>(*ptf.globalToLocalFaceZonePointMapPtr_);
    }
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
    //deleteDemandDrivenData(normalModelPtr_);
    //deleteDemandDrivenData(frictionModelPtr_);
    //deleteDemandDrivenData(allPointsDeformedPtr_);
    //deleteDemandDrivenData(oldZonePointsPtr_);
    //deleteDemandDrivenData(oldShadowZonePointsPtr_);
    deleteDemandDrivenData(zonePtr_);
    deleteDemandDrivenData(shadowZonePtr_);
    deleteDemandDrivenData(zoneToZonePtr_);
    //deleteDemandDrivenData(slaveFaceNormalsPtr_);
    deleteDemandDrivenData(globalFaceZonesPtr_);
    deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
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

        // PC: we need to remove all references to nonLinearGeometry
        // if
        // (
        //     master_
        //  && (
        //         nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN
        //      || nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
        //      || nonLinear_ == nonLinearGeometry::DEFORMED_LAGRANGIAN
        //     )
        // )
        // {
        //     deleteDemandDrivenData(oldZonePointsPtr_);
        //     deleteDemandDrivenData(oldShadowZonePointsPtr_);
        // }

        // Update penalty factor
        // PC: we should define a "newTimeStep" function and then the contact
        // law can do whatever it likes
        // if (normalModel().type() == standardPenalty::typeName && master_)
        // {
        //     standardPenalty& sp =
        //         refCast<standardPenalty>(normalModel());

        //     sp.updatePenaltyFactor();
        // }
        // Let the contact models know that it is a new time-step, in case they
        // need to update anything
        // TO BE DEFINED in normal and friction model
        //normalModel().newTimeStep();
        //frictionModel().newTimeStep();
    }

    // Move the master and slave zone to the deformed configuration
    //zone().clearGeom();
    //shadowZone().clearGeom();
    moveZonesToDeformedConfiguration();

    // Delete the zone-to-zone interpolator weights as the zones have moved
    zoneToZone().movePoints(tensorField(0), tensorField(0), vectorField(0));

    // Update the contact zones to the deformed position
    // This allows us to calculate the distances and also the weights for
    // interpolation
    //updateAllPointsDeformed();

    // Calculate and apply contact forces
    if (master_)
    {
        // Calculate the slave patch face unit normals as they are units by both
        // the normal and friction models
        const vectorField slavePatchFaceNormals =
            patchField
            (
                shadowPatchIndex(),
                shadowZoneIndex(),
                shadowZone().faceNormals()
            );

        // Calculate normal contact forces
        normalModel().correct
        (
            slavePatchFaceNormals,
            zoneToZone()
        );

        // Interpolate the master displacement increment to the slave patch as
        // it required by the friction model
        // Can I tidy this up? Maybe a function to pass this straight to the
        // friction model
        // Maybe just pass the slip!?
        if (db().foundObject<volVectorField>("DD"))
        {
            FatalErrorIn("solidContact::updateCoeffs()")
                << "wip: DD not ready yet!" << abort(FatalError);
        }
        const volVectorField& D = db().lookupObject<volVectorField>("D");
        const vectorField masterPatchDD =
            D.boundaryField()[patch().index()]
          - D.oldTime().boundaryField()[patch().index()];
        const vectorField slavePatchDD =
            D.boundaryField()[shadowPatchIndex()]
          - D.oldTime().boundaryField()[shadowPatchIndex()];

        // Master zone DD
        const vectorField masterZoneDD =
            zoneField
            (
                zoneIndex(),
                patch().index(),
                masterPatchDD
            );

        // Master patch DD interpolated to the slave patch
        const vectorField masterPatchDDInterpToSlavePatch =
            patchField
            (
                shadowPatchIndex(),
                shadowZoneIndex(),
                zoneToZone().masterToSlave(masterZoneDD)()
            );

        // Calculate friction contact forces
        frictionModel().correct
        (
            normalModel().slavePressure(),
            slavePatchFaceNormals,
            normalModel().areaInContact(),
            slavePatchDD,
            masterPatchDDInterpToSlavePatch
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
    os.writeKeyword("regionOfInterest")
        << regionOfInterest_ << token::END_STATEMENT << nl;

    if (master_)
    {
        os.writeKeyword("rigidMaster") << rigidMaster_
            << token::END_STATEMENT << nl;

        // os.writeKeyword("normalContactModel")
        //     << normalModel().type() << token::END_STATEMENT << nl;
        // normalModel().writeDict(os);

        // os.writeKeyword("frictionContactModel")
        //     << frictionModel().type() << token::END_STATEMENT << nl;
        // frictionModel().writeDict(os);
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
