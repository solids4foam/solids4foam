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

\*---------------------------------------------------------------------------*/

#include "standardPenaltyFriction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(standardPenaltyFriction, 0);
    addToRunTimeSelectionTable
    (
        frictionContactModel,
        standardPenaltyFriction,
        dictionary
    );
}


// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

void Foam::standardPenaltyFriction::calcFrictionPenaltyFactor()
{
    // Set penalty factor using a similar method to the normal
    // contact where we approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale

    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();

    // Lookup implicit stiffness = 2*mu + lambda, approximately equal to the
    // bulk modulus
    const volScalarField& impK = mesh_.lookupObject<volScalarField>("impK");

    // Avarage contact patch bulk modulus
    scalar masterK = gAverage(impK.boundaryField()[masterPatchIndex]);
    scalar slaveK = gAverage(impK.boundaryField()[slavePatchIndex]);

    // avarage contact patch shear modulus
    scalar modulus = 0.5*(masterK + slaveK);

    // average contact patch face area
    scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

    // average contact patch cell volume
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

    // approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    frictionPenaltyFactorPtr_ =
        new scalar(frictionPenaltyScale_*modulus*faceArea/cellVolume);

    Info<< "    friction penalty factor: " << *frictionPenaltyFactorPtr_
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::standardPenaltyFriction::standardPenaltyFriction
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID
)
:
    frictionContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZoneID,
        slaveFaceZoneID
    ),
    frictionContactModelDict_(dict.subDict(name + "FrictionModelDict")),
    frictionLawPtr_(NULL),
    mesh_(patch.boundaryMesh().mesh()),
    writeDebugFile_
    (
        frictionContactModelDict_.lookupOrDefault<Switch>
        (
            "writeDebugFile",
            false
        )
    ),
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    prevSlaveTraction_(slaveTraction_),
    slip_(slaveTraction_),
    frictionPenaltyFactorPtr_(NULL),
    frictionPenaltyScale_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("penaltyScale", 1.0)
    ),
    relaxFac_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.1
        )
    ),
    contactIterNum_(0),
    infoFreq_
    (
        frictionContactModelDict_.lookupOrDefault<int>("infoFrequency", 1)
    ),
    contactFilePtr_(NULL)
{
    // Create friction law
    frictionLawPtr_ =
        frictionLaw::New
        (
            frictionContactModelDict_.lookup("frictionLaw"),
            *this,
            frictionContactModelDict_
        ).ptr();

    // master proc open contact info file
    if (Pstream::master() && writeDebugFile_)
    {
        const word masterName = mesh_.boundary()[masterPatchID].name();
        const word slaveName = mesh_.boundary()[slavePatchID].name();

        const fileName contactFileDir = "contact";

        mkDir(contactFileDir);

        contactFilePtr_ =
            new OFstream
            (
                contactFileDir/
                "frictionContact_" + masterName + "_" + slaveName + ".txt"
            );

        OFstream& contactFile = *contactFilePtr_;

        int width = 20;
        contactFile
            << "time";
        contactFile.width(width);
        contactFile
            << "iterNum";
        contactFile.width(width);
        contactFile
            << "penaltyScale";
        contactFile.width(width);
        contactFile
            << "slipFaces";
        contactFile.width(width);
        contactFile
            << "stickFaces";
        contactFile.width(width);
        contactFile
            << "maxMagSlaveTraction";
        contactFile.width(width);
        contactFile
            << "maxMagSlip" << endl;
    }
}


// Construct as a copy
Foam::standardPenaltyFriction::standardPenaltyFriction
(
    const standardPenaltyFriction& fm
)
:
    frictionContactModel(fm),
    frictionContactModelDict_(fm.frictionContactModelDict_),
    frictionLawPtr_(fm.frictionLawPtr_->clone().ptr()),
    mesh_(fm.mesh_),
    writeDebugFile_(fm.writeDebugFile_),
    slaveTraction_(fm.slaveTraction_),
    prevSlaveTraction_(fm.prevSlaveTraction_),
    slip_(fm.slip_),
    frictionPenaltyFactorPtr_(NULL),
    frictionPenaltyScale_(fm.frictionPenaltyScale_),
    relaxFac_(fm.relaxFac_),
    contactIterNum_(fm.contactIterNum_),
    infoFreq_(fm.infoFreq_),
    contactFilePtr_(NULL)
{
    if (fm.frictionPenaltyFactorPtr_)
    {
        frictionPenaltyFactorPtr_ = new scalar(*fm.frictionPenaltyFactorPtr_);
    }

    if (fm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*fm.contactFilePtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standardPenaltyFriction::correct
(
    const vectorField& slavePressure,
    const vectorField& slaveFaceNormals,
    const scalarField& areaInContact,
    const vectorField& slaveDD,
    const vectorField& masterDDInterpToSlave
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();

    // Calculate slave shear traction increments
    const scalarField magSlavePressure = mag(slavePressure);
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    const scalarField oldStickSlip = stickSlip;
    scalar frictionPenaltyFac = frictionPenaltyFactor();
    scalar maxMagSlip = 0.0;
    scalarField slipTraction(magSlavePressure.size(), 0.0);

    forAll(magSlavePressure, faceI)
    {
        if (areaInContact[faceI] > SMALL)
        {
            // Compute slip as the we need the difference of DD between the
            // master and slave
            slip_[faceI] = slaveDD[faceI] - masterDDInterpToSlave[faceI];

            // The shear traction direction is got by removing the normal
            // component of the DD
            //     (I - sqr(n)) removes the normal
            //    sqr(n) would remove the shear
            slip_[faceI] = (I - sqr(slaveFaceNormals[faceI])) & slip_[faceI];

            slaveTraction_[faceI] = -frictionPenaltyFac*slip_[faceI];

            const scalar magSlip = mag(slip_[faceI]);
            maxMagSlip = max(maxMagSlip, magSlip);

            const scalar deltaT = mesh.time().deltaTValue();
            const vector slaveVelocity = slaveDD[faceI]/deltaT;
            const vector masterVelocity = masterDDInterpToSlave[faceI]/deltaT;

            // Traction to cause slipping i.e. the maximum shear traction the
            // face can hold for the given pressure, velocities, temperature,
            // etc.
            // Note: the actual friction law is implemented through the run-time
            // selectable frictionLaw
            slipTraction[faceI] =
                frictionLawPtr_->slipTraction
                (
                    magSlavePressure[faceI],    // Contact pressure
                    slip_[faceI],               // Slip vector
                    slaveVelocity,              // Velocity of slave face
                    masterVelocity,             // Velocity of master face
                    slavePatchIndex,            // Slave patch index
                    faceI                       // Local slave face ID
                );

            if ((mag(slaveTraction_[faceI]) - slipTraction[faceI]) > SMALL)
            {
                // Analogous to plasticity
                // slip is a combination of elastic slip and plastic slip
                // elastic slip should be zero but is finite due to penalty
                // stiffness. plastic slip is the permanent deformation
                slaveTraction_[faceI] =
                    slipTraction[faceI]*(-slip_[faceI]/magSlip);

                numSlipFaces++;
                stickSlip[faceI] = 1;
            }
            else
            {
                numStickFaces++;
                stickSlip[faceI] = 2;
            }
        }
        else
        {
            // No friction if pressure is negative or zero or face is not in
            // contact
            slaveTraction_[faceI] = vector::zero;
            slipTraction[faceI] = 0.0;
            stickSlip[faceI] = 0;
        }
    }

    // StickSlip field is just for visualisation but we will under-relax it to
    // allow us to see if a face is jumping between stick and slip
    stickSlip = relaxFac_*stickSlip + (1.0 - relaxFac_)*oldStickSlip;

    // Under-relax traction
    slaveTraction_ =
        relaxFac_*slaveTraction_ + (1.0 - relaxFac_)*prevSlaveTraction_;

    // Update the previous traction
    prevSlaveTraction_ = slaveTraction_;

    scalar maxMagSlaveTraction = 0.0;
    if (slaveTraction_.size() > 0)
    {
        maxMagSlaveTraction = max(mag(slaveTraction_));
    }
    reduce(maxMagSlaveTraction, maxOp<scalar>());
    reduce(numSlipFaces, sumOp<int>());
    reduce(numStickFaces, sumOp<int>());

    // Writes to contact info file
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
            << frictionPenaltyScale_;
        contactFile.width(width);
        contactFile
            << numSlipFaces;
        contactFile.width(width);
        contactFile
            << numStickFaces;
        contactFile.width(width);
        contactFile
            << maxMagSlaveTraction;
        contactFile.width(width);
        contactFile
            << maxMagSlip << endl;
    }
}


void Foam::standardPenaltyFriction::writeDict(Ostream& os) const
{
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
        << frictionContactModelDict_;
}


// ************************************************************************* //
