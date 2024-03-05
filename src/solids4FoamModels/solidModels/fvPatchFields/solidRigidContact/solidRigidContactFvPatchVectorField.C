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
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "lookupSolidModel.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidRigidContactFvPatchVectorField::movingMesh() const
{
    // Check if the solid model moves the mesh
    return lookupSolidModel(patch().boundaryMesh().mesh()).movingMesh();
}


void Foam::solidRigidContactFvPatchVectorField::makeShadowTriSurfNames
(
    const dictionary& dict
) const
{
    // Count the number of STLs
    const dictionary& triSurfDict = dict.subDict("triSurfaces");

    const wordList triSurfDictToc = triSurfDict.toc();

    if (triSurfDictToc.size() == 0)
    {
        FatalErrorIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::"
            "makeShadowTriSurfNames"
        )   << "No triSurfaces are defined!" << abort(FatalError);
    }

    shadowTriSurfNames_.setSize(triSurfDictToc.size());
    shadowTriSurfNames_ = triSurfDictToc;
}


void Foam::solidRigidContactFvPatchVectorField::makeNormalModels
(
    const dictionary& dict
) const
{
    normalModels_.setSize(shadowZones().size());

    forAll(normalModels_, shadPatchI)
    {
        // Check if only one shadow patch is defined
        const dictionary* contactDictPtr = NULL;
        if (normalModels_.size() == 1)
        {
            if (dict.found("normalContactModel") )
            {
                contactDictPtr = &dict;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + shadowTriSurfNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + shadowTriSurfNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;

        // Create contact model
        normalModels_.set
        (
            shadPatchI,
            normalContactModel::New
            (
                word(contactDict.lookup("normalContactModel")),
                // Pass any patch as it is just used to get access to the mesh
                patch(),
                contactDict,
                // We switch the master and slave for the contact model
                -1,                 // master in contact model
                patch().index(),    // slave in contact model
                shadowZones()[shadPatchI],
                zone().globalPatch()
            ).ptr()
        );
       }

    // Initialise penalty scales to -1
    Info<< "    Initialising stored previous normalPenaltyFactors" << endl;
    normalPenaltyFactors_.setSize(normalModels_.size(), -1);
}


void Foam::solidRigidContactFvPatchVectorField::makeFrictionModels
(
    const dictionary& dict
) const
{
    frictionModels_.setSize(shadowTriSurfNames().size());

    forAll (frictionModels_, triSurfI)
    {
        // Check if only one shadow patch is defined
        const dictionary* contactDictPtr = NULL;
        if (normalModels_.size() == 1)
        {
            if (dict.found("normalContactModel") )
            {
                contactDictPtr = &dict;
            }
            else
            {
                contactDictPtr =
                    &dict.subDict
                    (
                        patch().name() + "_to_"
                      + shadowTriSurfNames()[triSurfI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + shadowTriSurfNames()[triSurfI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;

        // Create contact model
        frictionModels_.set
        (
            triSurfI,
            frictionContactModel::New
            (
                word(contactDict.lookup("frictionContactModel")),
                patch(),
                contactDict,
                -1,
                patch().index()
            ).ptr()
        );
    }
}


void Foam::solidRigidContactFvPatchVectorField::clearOut()
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::clearOut()"
        )   << patch().name() << " : clearOut" << endl;
    }

    deleteDemandDrivenData(zonePtr_);
    // No need to clear the shadowZones as they are unchanging triSurfaces
    //shadowZones_.clear();
    zoneToZones_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidRigidContactFvPatchVectorField::solidRigidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    dict_(),
    writeZoneVTK_(false),
    writePointDistanceFields_(false),
    shadowTriSurfNames_(),
    normalModels_(),
    frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    shadowZones_(),
    shadowZonesMotion_(),
    zoneToZones_(),
    quickReject_(Foam::newGgiInterpolation::AABB),
    regionOfInterestTopCorner_(vector::max),
    regionOfInterestBottomCorner_(vector::min),
    regionOfInterest_(vector::min, vector::max),
    contact_(0),
    contactPerShadow_(),
    curTimeIndex_(-1)
{}


Foam::solidRigidContactFvPatchVectorField::solidRigidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    dict_(dict),
    writeZoneVTK_(dict.lookupOrDefault<Switch>("writeZoneVTK", false)),
    writePointDistanceFields_
    (
        dict.lookupOrDefault<Switch>("writePointDistanceFields", false)
    ),
    shadowTriSurfNames_(),
    normalModels_(),
    frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    shadowZones_(),
    shadowZonesMotion_(),
    zoneToZones_(),
    quickReject_
    (
        newGgiInterpolation::quickRejectNames_
        [
            dict.lookupOrDefault<word>("quickReject", "AABB")
        ]
    ),
    regionOfInterestTopCorner_
    (
        dict.lookupOrDefault<vector>
        (
            "regionOfInterestTopCorner",
            vector::max
        )
    ),
    regionOfInterestBottomCorner_
    (
        dict.lookupOrDefault<vector>
        (
            "regionOfInterestBottomCorner",
            vector::min
        )
    ),
    regionOfInterest_
    (
//        boundBox
//        (
//            regionOfInterestBottomCorner_,
//            regionOfInterestTopCorner_
//        )
        dict.lookupOrDefault<boundBox>
        (
            "regionOfInterest",
            boundBox(vector::min, vector::max)
        )
    ),
    contact_(patch().size(), 0.0),
    contactPerShadow_(),
    curTimeIndex_(-1)
{
    Info<< "Creating " << solidRigidContactFvPatchVectorField::typeName
        << " patch" << endl;

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


Foam::solidRigidContactFvPatchVectorField::solidRigidContactFvPatchVectorField
(
    const solidRigidContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    dict_(ptf.dict_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowTriSurfNames_(ptf.shadowTriSurfNames_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    shadowZonesMotion_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerShadow_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects: they will be re-created.
}

#ifndef OPENFOAM_ORG
Foam::solidRigidContactFvPatchVectorField::solidRigidContactFvPatchVectorField
(
    const solidRigidContactFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
    dict_(ptf.dict_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowTriSurfNames_(ptf.shadowTriSurfNames_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    shadowZonesMotion_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerShadow_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects
}
#endif

Foam::solidRigidContactFvPatchVectorField::solidRigidContactFvPatchVectorField
(
    const solidRigidContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    dict_(ptf.dict_),
    writeZoneVTK_(ptf.writeZoneVTK_),
    writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowTriSurfNames_(ptf.shadowTriSurfNames_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    shadowZonesMotion_(),
    zoneToZones_(),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    contact_(ptf.contact_),
    contactPerShadow_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solidRigidContactFvPatchVectorField::
~solidRigidContactFvPatchVectorField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidRigidContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::autoMap\n"
            "(\n"
            "    const fvPatchFieldMapper& m\n"
            ")"
        )   << patch().name() << " ::autoMap" << endl;
    }

    solidTractionFvPatchVectorField::autoMap(m);

    if (shadowTriSurfNames_.size() > 0)
    {
        contact_.autoMap(m);

        forAll(contactPerShadow_, shadI)
        {
            contactPerShadow_[shadI].autoMap(m);
        }

        // Let the contact models know about the mapping
        forAll(normalModels_, modelI)
        {
            normalModels_[modelI].autoMap(m);
            frictionModels_[modelI].autoMap(m);
        }
    }

    // Force all data to be re-created when needed
    clearOut();

    // Reset normal pelanty factors to reinitialise normal models
    normalPenaltyFactors_ =  -1;
}


void Foam::solidRigidContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);

    const solidRigidContactFvPatchVectorField& dmptf =
        refCast<const solidRigidContactFvPatchVectorField>(ptf);

    if (shadowTriSurfNames_.size() > 0)
    {
        contact_.rmap(dmptf.contact_, addr);

        forAll(contactPerShadow_, shadI)
        {
            contactPerShadow_[shadI].rmap(dmptf.contactPerShadow_[shadI], addr);
        }
    }
}


const Foam::wordList&
Foam::solidRigidContactFvPatchVectorField::shadowTriSurfNames() const
{
    if (shadowTriSurfNames_.size() == 0)
    {
        makeShadowTriSurfNames(dict_);
    }

    return shadowTriSurfNames_;
}


Foam::PtrList<Foam::normalContactModel>&
Foam::solidRigidContactFvPatchVectorField::normalModels()
{
    if (normalModels_.size() == 0)
    {
        makeNormalModels(dict_);
    }

    return normalModels_;
}


const Foam::PtrList<Foam::normalContactModel>&
Foam::solidRigidContactFvPatchVectorField::normalModels() const
{
    if (normalModels_.size() == 0)
    {
        makeNormalModels(dict_);
    }

    return normalModels_;
}


Foam::PtrList<Foam::frictionContactModel>&
Foam::solidRigidContactFvPatchVectorField::frictionModels()
{
    if (frictionModels_.size() == 0)
    {
        makeFrictionModels(dict_);
    }

    return frictionModels_;
}


const Foam::PtrList<Foam::frictionContactModel>&
Foam::solidRigidContactFvPatchVectorField::frictionModels() const
{
    if (frictionModels_.size() == 0)
    {
        makeFrictionModels(dict_);
    }

    return frictionModels_;
}


void Foam::solidRigidContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = db().time().timeIndex();

        // Let the contact models know that it is a new time-step, in case
        // they need to update anything
        forAll(shadowTriSurfNames(), triSurfI)
        {
            normalModels()[triSurfI].newTimeStep();
            frictionModels()[triSurfI].newTimeStep();

            // Force N^2 contact search at least once per time-step
            zoneToZones()[triSurfI].clearPrevCandidateMasterNeighbors();
        }

        // Move rigid triSurfaces to their new position
        moveTriSurfaces();
    }

    // Move the master and slave zone to the deformed configuration
    moveZonesToDeformedConfiguration();

    // Delete the zone-to-zone interpolator weights as the zones have moved
    const wordList& shadTriSurfNames = shadowTriSurfNames();
    forAll(shadTriSurfNames, triSurfI)
    {
        zoneToZones()[triSurfI].movePoints
        (
            tensorField(0), tensorField(0), vectorField(0)
        );
    }

    // Calculate and accumulate contact traction for all triSurfaces
    // Note: the contact models think the current patch is the slave

    traction() = vector::zero;

    forAll(shadTriSurfNames, triSurfI)
    {
        // Get the patch face deformed unit normals
        const vectorField patchFaceNormals =
            zone().globalFaceToPatch(zone().globalPatch().faceNormals());

        // Get the master patch displacement increment field

        vectorField patchDD(patch().size(), vector::zero);

        if (movingMesh())
        {
            // Updated Lagrangian, we will directly lookup the displacement
            // increment

            const volVectorField& DD =
                db().lookupObject<volVectorField>("DD");

            patchDD = DD.boundaryField()[patch().index()];
        }
        else
        {
            // We will lookup the total displacement and old total
            // displacement

            const volVectorField& D =
                db().lookupObject<volVectorField>("D");

            patchDD =
                D.boundaryField()[patch().index()]
              - D.oldTime().boundaryField()[patch().index()];
        }

        // Calculate the displacement increments on the triSurfaces and
        // interpolate it to the master surface
        const vectorField interpShadowPatchDD =
            zone().globalFaceToPatch
            (
                zoneToZones()[triSurfI].slaveToMaster
                (
                    shadowZonesMotion()[triSurfI].faceDisplacementIncrement()
                )()
            );

        // Calculate normal contact forces
        normalModels()[triSurfI].correct
        (
            patchFaceNormals,
            zone().globalPointToPatch
            (
                zoneToZones()[triSurfI].masterPointDistanceToIntersection()
            ),
            patchDD,
            interpShadowPatchDD
        );

        // Calculate friction contact forces
        frictionModels()[triSurfI].correct
        (
            normalModels()[triSurfI].slavePressure(),
            patchFaceNormals,
            normalModels()[triSurfI].areaInContact(),
            patchDD,
            interpShadowPatchDD
        );

        // Calculate the traction contribution for this triSurface
        const vectorField tractionForThisTriSurf =
            normalModels()[triSurfI].slavePressure()
          + frictionModels()[triSurfI].slaveTraction();

        // Add traction contribution from this triSurface
        traction() += tractionForThisTriSurf;

        // Update contactPerShadow field
        // Note: this is used by thermalContact to know which faces
        // are in contact
        const scalarField magTraction = mag(tractionForThisTriSurf);
        const scalar tol = 1e-6*gMax(magTraction);
        scalarField& contactForThisShadow = contactPerShadow()[triSurfI];
        forAll(contactForThisShadow, faceI)
        {
            if (magTraction[faceI] > tol)
            {
                contactForThisShadow[faceI] = 1.0;
            }
            else
            {
                contactForThisShadow[faceI] = 0.0;
            }
        }
    }

    // Accumulate the contact indicator field
    contact_ = 0.0;
    PtrList<scalarField>& contactPerShadow = this->contactPerShadow();
    forAll(contactPerShadow, shadI)
    {
        contact_ += contactPerShadow[shadI];
    }

    // Scale any face in contact with more than one shadow
    if (gMax(contact_) > (1.0 + SMALL))
    {
        forAll(contact_, faceI)
        {
            if (contact_[faceI] > (1.0 + SMALL))
            {
                // Update the contact weights corresponding to each shadow
                scalar sumContact = 0.0;
                forAll(contactPerShadow, shadI)
                {
                    contactPerShadow[shadI][faceI] /= contact_[faceI];
                    sumContact += contactPerShadow[shadI][faceI];
                }

                if (sumContact > (1.0 + SMALL))
                {
                    FatalErrorIn
                    (
                        "void solidRigidContactFvPatchVectorField::"
                        "updateCoeffs()"
                    )   << "There is a problem normalising the contact field"
                        << ", sumContact is: " << sumContact
                        << abort(FatalError);
                }

                // Reset accumulated contact face value to 1.0
                contact_[faceI] = 1.0;
            }
        }
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::solidRigidContactFvPatchVectorField::frictionHeatRate() const
{
    // Consider storing frictionHeatRate instead of recalculating multiple times

    // For now, we assume traction is constant over time-step
    // Todo: we should use trapezoidal rule
    vectorField curTraction(patch().size(), vector::zero);

    tmp<scalarField> tfrictionHeatRate
    (
        new scalarField(curTraction.size(), 0.0)
    );
    scalarField& frictionHeatRate = tfrictionHeatRate();

    forAll(shadowTriSurfNames(), triSurfI)
    {
        // Calculate slip

        const vectorField masterPatchSlip = frictionModels()[triSurfI].slip();

        const scalar deltaT =
            patch().boundaryMesh().mesh().time().deltaTValue();

        // Accumulate frictionHeatRate for each shadow patch

        // Rate of dissipated frictional energy for this timestep
        // The dot product of the traction vectors and the slip vectors gives
        // the dissipated frictional energy per unit area; which is always
        // positive
        frictionHeatRate += mag(traction() & (masterPatchSlip/deltaT));
    }

    return tfrictionHeatRate;
}


Foam::PtrList<Foam::scalarField>&
Foam::solidRigidContactFvPatchVectorField::contactPerShadow()
{
    if (contactPerShadow_.size() == 0)
    {
        calcContactPerShadow();
    }

    return contactPerShadow_;
}


const Foam::PtrList<Foam::scalarField>&
Foam::solidRigidContactFvPatchVectorField::contactPerShadow() const
{
    if (contactPerShadow_.size() == 0)
    {
        calcContactPerShadow();
    }

    return contactPerShadow_;
}


void Foam::solidRigidContactFvPatchVectorField::write(Ostream& os) const
{
    // If the shadowPatchIndices pointer is not set then we will assume that the
    // contact models were not created and nothing has changed; so we will just
    // output the input dict unchanged
    if (shadowTriSurfNames_.size() == 0)
    {
        // Overwrite fields in the dict
        dictionary& dict = const_cast<dictionary&>(dict_);

        dict.remove("gradient");
        dict.remove("value");
        dict.remove("traction");
        dict.remove("pressure");

        //dict.add("gradient", gradient());
        const vectorField& patchValue = *this;
        //dict.add("value", patchValue);
        //dict.add("traction", traction());
        //dict.add("pressure", pressure());

        // Write the dictionary
        dict_.write(os, false);

        gradient().writeEntry("gradient", os);
        patchValue.writeEntry("value", os);
        traction().writeEntry("traction", os);
        pressure().writeEntry("pressure", os);

        return;
    }

    solidTractionFvPatchVectorField::write(os);

    // Write triSurfaces subDict
    dict_.subDict("triSurfaces").write(os, false);

    os.writeKeyword("regionOfInterest")
        << regionOfInterest_ << token::END_STATEMENT << nl;
    os.writeKeyword("regionOfInterestTopCorner")
        << regionOfInterestTopCorner_ << token::END_STATEMENT << nl;
    os.writeKeyword("regionOfInterestBottomCorner")
        << regionOfInterestBottomCorner_ << token::END_STATEMENT << nl;
    os.writeKeyword("writeZoneVTK")
        << writeZoneVTK_ << token::END_STATEMENT << nl;
    os.writeKeyword("writePointDistanceFields")
        << writePointDistanceFields_ << token::END_STATEMENT << nl;
    os.writeKeyword("useNewPointDistanceMethod")
        << dict_.lookupOrDefault<Switch>
        (
            "useNewPointDistanceMethod", false
        )
        << token::END_STATEMENT << nl;

    os.writeKeyword("projectPointsToPatchBoundary")
        << dict_.lookupOrDefault<Switch>
        (
            "projectPointsToPatchBoundary", false
        )
        << token::END_STATEMENT << nl;

    os.writeKeyword("checkPointDistanceOrientations")
        << dict_.lookupOrDefault<Switch>
        (
            "checkPointDistanceOrientations", false
        )
        << token::END_STATEMENT << nl;

    os.writeKeyword("usePrevCandidateMasterNeighbors")
        << dict_.lookupOrDefault<Switch>
        (
            "usePrevCandidateMasterNeighbors", false
        )
        << token::END_STATEMENT << nl;

    if (shadowTriSurfNames_.size() == 1)
    {
        os.writeKeyword("normalContactModel")
            << normalModels()[0].type()
                << token::END_STATEMENT << nl;
        normalModels()[0].writeDict(os);

        os.writeKeyword("frictionContactModel")
            << frictionModels()[0].type()
            << token::END_STATEMENT << nl;
        frictionModels()[0].writeDict(os);
    }
    else
    {
        forAll(shadowTriSurfNames_, triSurfI)
        {
            os  << patch().name() << "_to_"
                << shadowTriSurfNames_[triSurfI] << "_dict" << nl
                << '{' << endl;

            os.writeKeyword("normalContactModel")
                << normalModels()[triSurfI].type()
                << token::END_STATEMENT << nl;
            normalModels()[triSurfI].writeDict(os);

            os.writeKeyword("frictionContactModel")
                << frictionModels()[triSurfI].type()
                << token::END_STATEMENT << nl;
            frictionModels()[triSurfI].writeDict(os);

            os  << '}' << endl;
        }
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

            zone().globalPatch().writeVTK("zone_" + timeName);
        }
    }

    // Write out deformed triSurfaces
    const fileName triSurfDir(db().time().path()/"deformedTriSurfaces");
    if (!isDir(triSurfDir))
    {
        mkDir(triSurfDir);
    }

    const word timeIndex = name(db().time().timeIndex());

    forAll(shadowZones(), triSurfI)
    {
        shadowZones()[triSurfI].writeVTK
        (
            fileName(triSurfDir/shadowTriSurfNames()[triSurfI] + timeIndex)
        );
    }

    // Write out point distance fields for master and slave
    if (writePointDistanceFields_)
    {
        if (normalModels().size() != 1)
        {
            FatalErrorIn
            (
                "void Foam::solidRigidContactFvPatchVectorField::"
                "write(Ostream& os) const"
            )   << "The 'writePointDistanceFields' is currently only "
                << "implemented for one-to-one contact"
                << abort(FatalError);
        }

        // Take a reference to the mesh for convenience
        const polyMesh& mesh = patch().patch().boundaryMesh().mesh();

        // Create the point mesh, which is needed for the point field
        const pointMesh& pMesh = pointMesh::New(mesh);

        // Create the point distance fields

        pointScalarField dist
        (
            IOobject
            (
                "pointDistance",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

        pointVectorField distVecs
        (
            IOobject
            (
                "pointDistanceVectors",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimless, vector::zero)
        );

        // Lookup the master point distance to intersection
        const scalarField pd =
            zone().globalPointToPatch
            (
                zoneToZones()[0].masterPointDistanceToIntersection()
            );
        const vectorField pdVecs =
            zone().globalPointToPatch
            (
                zoneToZones()[0].masterPointDistanceVectorsToIntersection()
            );

        // Transfer the patch point distances into the dist point field
        const labelList& meshPoints = patch().patch().meshPoints();
        forAll(pd, pI)
        {
            const label pointID = meshPoints[pI];
            dist[pointID] = pd[pI];
            distVecs[pointID] = pdVecs[pI];
        }

        // Write the field
        InfoIn
        (
            "void Foam::solidRigidContactFvPatchVectorField::"
            "write(Ostream& os) const"
        )   << "Writing point distance fields: " << dist.name()
            << " and " << distVecs.name() << endl;
        dist.write();
        distVecs.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        solidRigidContactFvPatchVectorField
    );
}


// ************************************************************************* //
