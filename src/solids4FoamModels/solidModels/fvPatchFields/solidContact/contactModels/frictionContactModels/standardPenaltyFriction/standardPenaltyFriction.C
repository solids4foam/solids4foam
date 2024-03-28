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

#include "standardPenaltyFriction.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

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

#ifdef OPENFOAM_NOT_EXTEND
    typedef labelUList unallocLabelList;
#endif
}


// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

Foam::scalar Foam::standardPenaltyFriction::frictionPenaltyFactor()
{
    if (frictionPenaltyFactor_ < -SMALL)
    {
        calcFrictionPenaltyFactor();
    }

    return frictionPenaltyFactor_;
}


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

    // Note: for solidRigidContact, only the master index is set
    if (masterPatchIndex > -1 && slavePatchIndex > -1)
    {
        // Avarage contact patch bulk modulus
        const scalar masterK = gAverage(impK.boundaryField()[masterPatchIndex]);
        const scalar slaveK = gAverage(impK.boundaryField()[slavePatchIndex]);

        // avarage contact patch shear modulus
        const scalar modulus = 0.5*(masterK + slaveK);

        // average contact patch face area
        const scalar masterMagSf =
            gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
        const scalar slaveMagSf =
            gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
        const scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

        // average contact patch cell volume
        scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
        scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);

#ifdef OPENFOAM_NOT_EXTEND
    const DimensionedField<scalar, volMesh>& V = mesh_.V();
#else
    const volScalarField::DimensionedInternalField& V = mesh_.V();
#endif

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
        const scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

        // approximate penalty factor based on Hallquist et al.
        // we approximate penalty factor for traction instead of force
        frictionPenaltyFactor_ =
            frictionPenaltyScale_*modulus*faceArea/cellVolume;
    }
    else if (slavePatchIndex > -1)
    {
        // Avarage contact patch bulk modulus
        const scalar modulus = gAverage(impK.boundaryField()[slavePatchIndex]);

        // average contact patch face area
        const scalar faceArea =
            gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);

        // average contact patch cell volume
        scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);

#ifdef OPENFOAM_NOT_EXTEND
    const DimensionedField<scalar, volMesh>& V = mesh_.V();
#else
    const volScalarField::DimensionedInternalField& V = mesh_.V();
#endif

        {
            const unallocLabelList& faceCells =
                mesh_.boundary()[slavePatchIndex].faceCells();
            forAll(mesh_.boundary()[slavePatchIndex], facei)
            {
                slaveV[facei] = V[faceCells[facei]];
            }
        }
        const scalar cellVolume = gAverage(slaveV);

        // approximate penalty factor based on Hallquist et al.
        // we approximate penalty factor for traction instead of force
        frictionPenaltyFactor_ =
            frictionPenaltyScale_*modulus*faceArea/cellVolume;
    }
    else
    {
        FatalErrorIn("void standardPenaltyFriction::calcPenaltyFactor() const")
            << "This is unexpected! Neither the master nor slave index are set!"
            << abort(FatalError);
    }

    Info<< "    friction penalty factor: " << frictionPenaltyFactor_
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
    const label slavePatchID
)
:
    frictionContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID
    ),
    frictionContactModelDict_(dict.subDict(name + "FrictionModelDict")),
    frictionLawPtr_(),
    mesh_(patch.boundaryMesh().mesh()),
    tractionVolField_
    (
        IOobject
        (
            "shearTraction_" + mesh_.boundaryMesh()[slavePatchID].name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimPressure, vector::zero)
    ),
    slipOnSlave_
    (
        tractionVolField_.boundaryField()[slavePatchID].size(),
        vector::zero
    ),
    slipOnMaster_
    (
        tractionVolField_.boundaryField()[masterPatchID].size(),
        vector::zero
    ),
    frictionPenaltyFactor_(-1),
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
    contactIterNum_(0)
{
    // Create friction law
    frictionLawPtr_.set
    (
        frictionLaw::New
        (
            word(frictionContactModelDict_.lookup("frictionLaw")),
            *this,
            frictionContactModelDict_
        ).ptr()
    );
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
    tractionVolField_(fm.tractionVolField_),
    slipOnSlave_(fm.slipOnSlave_),
    slipOnMaster_(fm.slipOnMaster_),
    frictionPenaltyFactor_(fm.frictionPenaltyFactor_),
    frictionPenaltyScale_(fm.frictionPenaltyScale_),
    relaxFac_(fm.relaxFac_),
    contactIterNum_(fm.contactIterNum_)
{}


// * * * * * * * * * * * * * * * *  Destructor  * * * * * * * * * * * * * * //

Foam::standardPenaltyFriction::~standardPenaltyFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standardPenaltyFriction::correct
(
    const vectorField& slavePressure,
    const vectorField& slaveFaceNormals,
    const scalarField& slavePatchAreaInContact,
    const vectorField& slaveDD,
    const vectorField& masterDDInterpToSlave
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();

    // Calculate slave shear traction increments
    const scalarField magSlavePressure(mag(slavePressure));
    // label numSlipFaces = 0;
    // label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    const scalarField oldStickSlip = stickSlip;
    const scalar frictionPenaltyFac = frictionPenaltyFactor();
    scalar maxMagSlip = 0.0;
    scalarField slipTraction(magSlavePressure.size(), 0.0);
    vectorField newSlaveTraction(slipTraction.size(), vector::zero);

    forAll(magSlavePressure, faceI)
    {
        if (slavePatchAreaInContact[faceI] > SMALL)
        {
            // Compute slip as the we need the difference of DD between the
            // master and slave
            slipOnSlave_[faceI] = slaveDD[faceI] - masterDDInterpToSlave[faceI];

            // The shear traction direction is got by removing the normal
            // component of the DD
            //     (I - sqr(n)) removes the normal
            //    sqr(n) would remove the shear
            slipOnSlave_[faceI] =
                (I - sqr(slaveFaceNormals[faceI])) & slipOnSlave_[faceI];

            newSlaveTraction[faceI] = -frictionPenaltyFac*slipOnSlave_[faceI];

            const scalar magSlip = mag(slipOnSlave_[faceI]);
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
                frictionLawPtr_().slipTraction
                (
                    magSlavePressure[faceI],    // Contact pressure
                    slipOnSlave_[faceI],               // Slip vector
                    slaveVelocity,              // Velocity of slave face
                    masterVelocity,             // Velocity of master face
                    slavePatchIndex,            // Slave patch index
                    faceI                       // Local slave face ID
                );

            if ((mag(newSlaveTraction[faceI]) - slipTraction[faceI]) > SMALL)
            {
                // Analogous to plasticity
                // slip is a combination of elastic slip and plastic slip
                // elastic slip should be zero but is finite due to penalty
                // stiffness. plastic slip is the permanent deformation
                newSlaveTraction[faceI] =
                    slipTraction[faceI]*(-slipOnSlave_[faceI]/magSlip);

                // numSlipFaces++;
                stickSlip[faceI] = 1;
            }
            else
            {
                // numStickFaces++;
                stickSlip[faceI] = 2;
            }
        }
        else
        {
            // No friction if pressure is negative or zero or face is not in
            // contact
            newSlaveTraction[faceI] = vector::zero;
            slipTraction[faceI] = 0.0;
            stickSlip[faceI] = 0;
        }
    }

    // StickSlip field is just for visualisation but we will under-relax it to
    // allow us to see if a face is jumping between stick and slip
    stickSlip = relaxFac_*stickSlip + (1.0 - relaxFac_)*oldStickSlip;

    // Under-relax traction
    slaveTraction() =
        relaxFac_*newSlaveTraction + (1.0 - relaxFac_)*slaveTraction();
}


void Foam::standardPenaltyFriction::correct
(
    const vectorField& patchPressure,
    const vectorField& patchFaceNormals,
    const scalarField& patchAreaInContact,
    const vectorField& DD,
    const vectorField& shadowDDInterpToPatch,
    const bool master
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;
    label patchIndex = -1;

    if (master)
    {
        patchIndex = masterPatchID();
    }
    else
    {
        patchIndex = slavePatchID();
    }

    const scalarField magPressure(mag(patchPressure));

    const scalar frictionPenaltyFac = frictionPenaltyFactor();

    scalarField slipTraction(magPressure.size(), 0.0);
    vectorField newTraction(slipTraction.size(), vector::zero);
    scalar maxMagSlip = 0.0;


    forAll(magPressure, faceI)
    {
        if (patchAreaInContact[faceI] > SMALL)
        {
            vector& faceSlip = slip(master)[faceI];

            // Compute slip as the we need the difference of DD between patches
            faceSlip = DD[faceI] - shadowDDInterpToPatch[faceI];

            // The shear traction direction is got by removing the normal
            // component of the DD
            faceSlip = (I - sqr(patchFaceNormals[faceI])) & faceSlip;

            newTraction[faceI] = -frictionPenaltyFac*faceSlip;

            const scalar magSlip = mag(faceSlip);
            maxMagSlip = max(maxMagSlip, magSlip);

            const scalar deltaT = mesh.time().deltaTValue();
            const vector faceVelocity = DD[faceI]/deltaT;
            const vector shadowFaceVelocity =
                shadowDDInterpToPatch[faceI]/deltaT;

            // Traction to cause slipping i.e. the maximum shear traction the
            // face can hold for the given pressure, velocities, temperature,
            // etc.
            // Note: the actual friction law is implemented through the run-time
            // selectable frictionLaw
            slipTraction[faceI] =
                frictionLawPtr_().slipTraction
                (
                    magPressure[faceI],         // Contact pressure
                    faceSlip,                   // Slip vector
                    faceVelocity,               // Face velocity
                    shadowFaceVelocity,         // Shadow face velocity
                    patchIndex,                 // Patch index
                    faceI                       // Local face ID
                );

            if ((mag(newTraction[faceI]) - slipTraction[faceI]) > SMALL)
            {
                // Analogous to plasticity
                // slip is a combination of elastic slip and plastic slip
                // elastic slip should be zero but is finite due to penalty
                // stiffness. plastic slip is the permanent deformation
                newTraction[faceI] =
                    slipTraction[faceI]*(-faceSlip/magSlip);

                // stickSlip field is not updated. Should be added.
            }
        }
        else
        {
            // No friction if pressure is negative or zero or face is not in
            // contact
            newTraction[faceI] = vector::zero;
            slipTraction[faceI] = 0.0;
        }
    }

    vectorField& patchTraction =
        tractionVolField_.boundaryField()[patchIndex];

    // Under-relax traction
    patchTraction =
        relaxFac_*newTraction + (1.0 - relaxFac_)*patchTraction;
}


void Foam::standardPenaltyFriction::autoMap(const fvPatchFieldMapper& m)
{
    frictionContactModel::autoMap(m);

    if (debug)
    {
        InfoIn
        (
            "void standardPenaltyFriction::autoMap(const fvPatchFieldMapper& m)"
        )   << "autoMap" << endl;
    }

#ifdef OPENFOAM_ORG
    m(slipOnSlave_, slipOnSlave_);
    m(slipOnMaster_, slipOnMaster_);
#else
    slipOnSlave_.autoMap(m);
    slipOnMaster_.autoMap(m);
#endif

    // The internal fields for the volFields should always be zero
    // We will reset them as they may not be zero after field advection
#ifdef OPENFOAM_NOT_EXTEND
    tractionVolField_.primitiveFieldRef() = vector::zero;
#else
    tractionVolField_.internalField() = vector::zero;
#endif
}


void Foam::standardPenaltyFriction::writeDict(Ostream& os) const
{
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
        << frictionContactModelDict_;
}


// ************************************************************************* //
