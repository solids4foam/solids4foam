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
    standardPenalty

\*---------------------------------------------------------------------------*/

#include "standardPenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(standardPenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, standardPenalty, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


void standardPenalty::calcPenaltyFactor() const
{
    // Set penalty factor

    // Approximate penaltyFactor from the mechanical properties
    // This can then be scaled using the penaltyScale

    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();

    // Lookup implicit stiffness = 2*mu + lambda, approximately equal to the
    // bulk modulus
    const volScalarField& impK = mesh_.lookupObject<volScalarField>("impK");

    // Avarage contact patch bulk modulus
    scalar masterK = gAverage(impK.boundaryField()[masterPatchIndex]);
    scalar slaveK = gAverage(impK.boundaryField()[slavePatchIndex]);
    scalar bulkModulus = 0.5*(masterK + slaveK);

    // Average contact patch face area
    scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

    // Average contact patch cell volume
    scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    const volScalarField::DimensionedInternalField& V = mesh_.V();
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[masterPatchIndex].faceCells();
        forAll(mesh_.boundary()[masterPatchIndex], facei)
        {
            masterV[facei] = V[faceCells[facei]];
        }
    }
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[slavePatchIndex].faceCells();
        forAll(mesh_.boundary()[slavePatchIndex], facei)
        {
            slaveV[facei] = V[faceCells[facei]];
        }
    }
    scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

    // Approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    penaltyFactorPtr_ =
        new scalar(penaltyScale_*bulkModulus*faceArea/cellVolume);

    Info<< "    normal penalty factor: " << *penaltyFactorPtr_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPenalty::standardPenalty
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
:
    normalContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZoneID,
        slaveFaceZoneID,
        masterFaceZonePatch,
        slaveFaceZonePatch
    ),
    normalContactModelDict_(dict.subDict(name + "NormalModelDict")),
    mesh_(patch.boundaryMesh().mesh()),
    writeDebugFile_
    (
        normalContactModelDict_.lookupOrDefault<Switch>("writeDebugFile", false)
    ),
    slavePressure_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    globalSlavePointPenetration_(slaveFaceZonePatch.points().size(), 0.0),
    slavePointPenetration_
    (
        mesh_.boundaryMesh()[slavePatchID].localPoints().size(), 0.0
    ),
    areaInContact_(slavePressure_.size(), 0.0),
    penaltyFactorPtr_(NULL),
    penaltyScale_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("penaltyScale", 1.0)
    ),
    relaxFac_
    (
        normalContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.02
        )
    ),
    totalSlavePointPressure_(slavePointPenetration_.size(), 0.0),
    contactIterNum_(0),
    infoFreq_(normalContactModelDict_.lookupOrDefault<int>("infoFrequency", 1)),
    contactFilePtr_(NULL)
{
    if (Pstream::master() && writeDebugFile_)
    {
        const word masterName = mesh_.boundary()[masterPatchID].name();
        const word slaveName = mesh_.boundary()[slavePatchID].name();

        const fileName contactFileDir = "contact";

        mkDir(contactFileDir);

        contactFilePtr_ =
            new OFstream
            (
                contactFileDir/"normalContact_"+masterName+"_"+slaveName+".txt"
            );

        OFstream& contactFile = *contactFilePtr_;

        int width = 20;
        contactFile
            << "Time";
        contactFile.width(width);
        contactFile
        << "ContactIter";
        contactFile.width(width);
        contactFile
            << "slaveContactVer";
        contactFile.width(width);
        contactFile
            << "penetration";
        contactFile.width(width);
        contactFile
            << "maxSlaveTrac";
        contactFile.width(width);
        contactFile
            << "numMissedPoints" << endl;
    }
}


standardPenalty::standardPenalty(const standardPenalty& nm)
:
    normalContactModel(nm),
    normalContactModelDict_(nm.normalContactModelDict_),
    mesh_(nm.mesh_),
    writeDebugFile_(nm.writeDebugFile_),
    slavePressure_(nm.slavePressure_),
    globalSlavePointPenetration_(nm.globalSlavePointPenetration_),
    slavePointPenetration_(nm.slavePointPenetration_),
    areaInContact_(nm.areaInContact_),
    penaltyFactorPtr_(NULL),
    penaltyScale_(nm.penaltyScale_),
    relaxFac_(nm.relaxFac_),
    contactIterNum_(nm.contactIterNum_),
    infoFreq_(nm.infoFreq_),
    contactFilePtr_(NULL)
{
    if (nm.penaltyFactorPtr_)
    {
        penaltyFactorPtr_ = new scalar(*nm.penaltyFactorPtr_);
    }

    if (nm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*nm.contactFilePtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //


void standardPenalty::correct
(
    const vectorField& slavePatchFaceNormals,
    const extendedGgiStandAlonePatchInterpolation& zoneToZone
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();

    // Calculate the point distances for the slave zone
    globalSlavePointPenetration_ =
        zoneToZone.slavePointDistanceToIntersection();

    // Transfer the global zone distances to patch patch distances
    forAll(slavePointPenetration_, pointI)
    {
        // The local point values are kept at the start of the global field
        slavePointPenetration_[pointI] = globalSlavePointPenetration_[pointI];
    }

    // Calculate area in contact for slave patch

    const faceList& slavePatchLocalFaces =
        mesh.boundaryMesh()[slavePatchIndex].localFaces();

    const pointField& slavePatchLocalPoints =
        mesh.boundaryMesh()[slavePatchIndex].localPoints();

    forAll(slavePatchLocalFaces, faceI)
    {
        areaInContact_[faceI] =
            slavePatchLocalFaces[faceI].areaInContact
            (
                slavePatchLocalPoints,
                slavePointPenetration_
            );
    }

    // Calculate traction increments

    int numSlaveContactPoints = 0;
    const scalar minSlavePointPenetration = gMin(globalSlavePointPenetration_);
    const scalar penaltyFac = penaltyFactor();

    forAll(totalSlavePointPressure_, pointI)
    {
        // If a point has penetrated (i.e. if the penetration is negative),
        if (slavePointPenetration_[pointI] < 0.0)
        {
            // Count points in contact
            numSlaveContactPoints++;

            // The force is linearly proportional the penetration, like a spring
            totalSlavePointPressure_[pointI] =
                penaltyFac*slavePointPenetration_[pointI];
        }
        else
        {
            totalSlavePointPressure_[pointI] = 0.0;
        }
    }


    // Interpolate point traction to faces

    // Create local patch interpolation: No need to interpolate using the entire
    // face zone patch
    // This is using the original mesh weights; it would be better to use the
    // deformed patch, and it might help
    // This could be changed to used the current weights
    // WarningIn("standardPenalty::correct(...)")
    //    << "Philip: check pointToFace interpolator" << endl;
    primitivePatchInterpolation localSlaveInterpolator
    (
        mesh.boundaryMesh()[slavePatchIndex]
    );

    vectorField newSlavePressure =
        localSlaveInterpolator.pointToFaceInterpolate<scalar>
        (
            totalSlavePointPressure_
        )*slavePatchFaceNormals;

    // Under-relax pressure
    slavePressure_ =
        relaxFac_*newSlavePressure + (1.0 - relaxFac_)*slavePressure_;

    scalar maxMagSlaveTraction = 0.0;
    if (slavePressure_.size() > 0)
    {
        maxMagSlaveTraction = max(mag(slavePressure_));
    }
    reduce(maxMagSlaveTraction, maxOp<scalar>());
    reduce(numSlaveContactPoints, sumOp<int>());

    // Write to contact file
    if
    (
        Pstream::master()
     && (contactIterNum_++ %  infoFreq_ == 0)
     && writeDebugFile_
    )
    {
        OFstream& contactFile = *contactFilePtr_;
        int width = 20;
        contactFile
            << mesh.time().value();
        contactFile.width(width);
        contactFile
            << contactIterNum_;
        contactFile.width(width);
        contactFile
            << numSlaveContactPoints;
        contactFile.width(width);
        contactFile
            << minSlavePointPenetration;
        contactFile.width(width);
        contactFile
            << maxMagSlaveTraction;
        contactFile.width(width);
        contactFile
            << endl;
    }
}


scalar standardPenalty::penaltyFactor() const
{
    if (!penaltyFactorPtr_)
    {
        calcPenaltyFactor();
    }

    return *penaltyFactorPtr_;
}


void standardPenalty::writeDict(Ostream& os) const
{
    word keyword(name() + "NormalModelDict");
    os.writeKeyword(keyword)
        << normalContactModelDict_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
