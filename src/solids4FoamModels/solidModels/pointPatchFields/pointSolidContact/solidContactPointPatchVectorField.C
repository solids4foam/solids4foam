/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "solidContactPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#include "lookupSolidModel.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidContactPointPatchVectorField::movingMesh() const
{
    // Check if the solid model moves the mesh
    return lookupSolidModel(db()).movingMesh();
}


void Foam::solidContactPointPatchVectorField::makeShadowPatchNames
(
    const dictionary& dict
) const
{
    if (master_)
    {
        // Check if only one shadow patch is specified
        if (dict.found("shadowPatch"))
        {
            Info<< "Reading individual shadowPatch" << endl;

            // Just one shadow patch
            shadowPatchNames_.setSize(1);
            shadowPatchNames_[0] = word(dict.lookup("shadowPatch"));
        }
        else if (dict.found("shadowPatches"))
        {
            Info<< "Reading list of shadowPatches" << endl;

            // Shadow patches defined as a list
            shadowPatchNames_ = wordList(dict.lookup("shadowPatches"));
        }
        else
        {
            FatalErrorIn
            (
                "void solidContactPointPatchVectorField::"
                " makeShadowPatchNames() const"
            )   << "'shadowPatch' OR 'shadowPatches' should be defined"
                << abort(FatalError);
        }

        // It is an error to defined both shadowPatch and shadowPatches
        if (dict.found("shadowPatch") && dict.found("shadowPatches"))
        {
            FatalErrorIn
            (
                "void solidContactPointPatchVectorField::"
                " makeShadowPatchNames() const"
            )   << "'shadowPatch' OR 'shadowPatches' should be defined: "
                << "not both!" << abort(FatalError);
        }
    }
    else
    {
        // If this is not the master then we will assume there is only one
        // shadow i.e. the master is the shadow
        shadowPatchNames_.setSize(1);
        shadowPatchNames_[0] = word(dict.lookup("shadowPatch"));
    }
}


void Foam::solidContactPointPatchVectorField::calcShadowPatchIndices() const
{
    if (shadowPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "void solidContactPointPatchVectorField::"
            " calcShadowPatchIndices() const"
        )   << "pointer already set" << abort(FatalError);
    }

    shadowPatchIndicesPtr_ = new labelList(shadowPatchNames().size(), -1);
    labelList& shadowPatchIndices = *shadowPatchIndicesPtr_;

    forAll(shadowPatchIndices, shadPatchI)
    {
        shadowPatchIndices[shadPatchI] =
            patch().boundaryMesh().findPatchID
            (
                shadowPatchNames_[shadPatchI]
            );

        if (shadowPatchIndices[shadPatchI] == -1)
        {
            FatalErrorIn
            (
                "void solidContactPointPatchVectorField::"
                " calcShadowPatchIndices() const"
            )   << "shadowPatch " << shadowPatchNames_[shadPatchI]
                << " not found!" << abort(FatalError);
        }
    }
}


void Foam::solidContactPointPatchVectorField::makeNormalModels
(
    const dictionary& dict
) const
{
    normalModels_.setSize(shadowPatchNames().size());

    forAll (normalModels_, shadPatchI)
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
                      + shadowPatchNames()[shadPatchI] + "_dict"
                    );
            }
        }
        else
        {
            contactDictPtr =
                &dict.subDict
                (
                    patch().name() + "_to_"
                  + shadowPatchNames()[shadPatchI] + "_dict"
                );
        }
        const dictionary& contactDict = *contactDictPtr;

        // Create contact model
        normalModels_.set
        (
            shadPatchI,
            pointNormalContactModel::New
            (
                word(contactDict.lookup("normalContactModel")),
                lookupSolidModel(db()).mesh(),
                contactDict,
                patch().index(),                  // master
                shadowPatchIndices()[shadPatchI], // slave
                zone().globalPatch(),
                shadowZones()[shadPatchI].globalPatch()
            ).ptr()
        );
    }

    // Initialise penalty scales to -1
    Info<< "    Initialising stored previous normalPenaltyFactors" << endl;
    normalPenaltyFactors_.setSize(normalModels_.size(), -1);
}


// void Foam::solidContactPointPatchVectorField::makeFrictionModels
// (
//     const dictionary& dict
// ) const
// {
//     frictionModels_.setSize(shadowPatchNames().size());

//     forAll (frictionModels_, shadPatchI)
//     {
//         // Check if only one shadow patch is defined
//         const dictionary* contactDictPtr = NULL;
//         if (frictionModels_.size() == 1)
//         {
//             if (dict.found("frictionContactModel") )
//             {
//                 contactDictPtr = &dict;
//             }
//             else
//             {
//                 contactDictPtr =
//                     &dict.subDict
//                     (
//                         patch().name() + "_to_"
//                       + shadowPatchNames()[shadPatchI] + "_dict"
//                     );
//             }
//         }
//         else
//         {
//             contactDictPtr =
//                 &dict.subDict
//                 (
//                     patch().name() + "_to_"
//                   + shadowPatchNames()[shadPatchI] + "_dict"
//                 );
//         }
//         const dictionary& contactDict = *contactDictPtr;

//         // Create contact model
//         frictionModels_.set
//         (
//             shadPatchI,
//             frictionContactModel::New
//             (
//                 word(contactDict.lookup("frictionContactModel")),
//                 patch(),
//                 contactDict,
//                 patch().index(),                 // master
//                 shadowPatchIndices()[shadPatchI] // slave
//             ).ptr()
//         );
//     }
// }


void Foam::solidContactPointPatchVectorField::clearOut()
{
    if (debug)
    {
        InfoIn
        (
            "void Foam::solidContactPointPatchVectorField::clearOut()"
        )   << patch().name() << " : clearOut" << endl;
    }

    deleteDemandDrivenData(shadowPatchIndicesPtr_);
    deleteDemandDrivenData(zonePtr_);
    shadowZones_.clear();
    zoneToZones_.clear();
    // scaleTractionFieldPtr_.clear();
}

// Foam::scalarField
// Foam::solidContactPointPatchVectorField::scaleTractionField() const
// {
//     if (scaleTractionFieldPtr_.empty())
//     {
//         makeScaleTractionField();
//     }

//     return scaleTractionFieldPtr_();
// }


// void Foam::solidContactPointPatchVectorField::makeScaleTractionField() const
// {
//     if (scaleTractionFieldPtr_.valid())
//     {
//         FatalErrorIn(type() + "::makeScaleTractionField()")
//             << "Pointer already set!" << abort(FatalError);
//     }

//     scaleTractionFieldPtr_.set(new scalarField(patch().size(), 1.0));
//     scalarField& scaleTractionField = scaleTractionFieldPtr_();

//     // Find all faces on the patch that are adjacent to faces on the
//     // downstream patch

//     const fvMesh& mesh = patch().boundaryMesh().mesh();
//     const scalar scaleFactor =
//         readScalar(dict_.lookup("downstreamScaleFactor"));

//     // Downstream patch name
//     const word patchName = word(dict_.lookup("downstreamPatchName"));
//     const label patchID = mesh.boundaryMesh().findPatchID(patchName);
//     if (patchID == -1)
//     {
//         FatalErrorIn(type() + "::makeScaleTractionField()")
//             << "Cannot find patch " << patchName << abort(FatalError);
//     }

//     const unallocLabelList& faceCells = patch().faceCells();
//     const cellList& cells = mesh.cells();

//     forAll(scaleTractionField, fI)
//     {
//         const label cellID = faceCells[fI];
//         const cell& curCell = cells[cellID];

//         forAll(curCell, cfI)
//         {
//             const label cellFaceID = curCell[cfI];

//             if (!mesh.isInternalFace(cellFaceID))
//             {
//                 if (mesh.boundaryMesh().whichPatch(cellFaceID) == patchID)
//                 {
//                     scaleTractionField[fI] = scaleFactor;
//                     break;
//                 }
//             }
//         }
//     }
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidContactPointPatchVectorField::solidContactPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(p, iF),
    dict_(),
    master_(false),
    // writePointDistanceFields_(false),
    shadowPatchNames_(),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(false),
    normalModels_(),
    // frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
#ifdef FOAMEXTEND
    writeZoneVTK_(false),
    quickReject_(Foam::newGgiInterpolation::AABB),
    regionOfInterestTopCorner_(vector::max),
    regionOfInterestBottomCorner_(vector::min),
    regionOfInterest_(vector::min, vector::max),
#endif
    // contact_(0),
    // contactPerShadow_(),
    // scaleFaceTractionsNearDownstreamPatch_(false),
    // scaleTractionFieldPtr_(),
    curTimeIndex_(-1)
{}


solidContactPointPatchVectorField::solidContactPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    solidTractionPointPatchVectorField(p, iF),
    dict_(dict),
    master_(dict.lookup("master")),
    // writePointDistanceFields_
    // (
    //     dict.lookupOrDefault<Switch>("writePointDistanceFields", false)
    // ),
    shadowPatchNames_(),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(false),
    normalModels_(),
    // frictionModels_(),
    normalPenaltyFactors_(),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
#ifdef FOAMEXTEND
    writeZoneVTK_(dict.lookupOrDefault<Switch>("writeZoneVTK", false)),
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
        boundBox
        (
            regionOfInterestBottomCorner_,
            regionOfInterestTopCorner_
        )
//        dict.lookupOrDefault<boundBox>
//        (
//            "regionOfInterest",
//            boundBox(vector::min, vector::max)
//        )
    ),
#endif
    // contact_(patch().size(), 0.0),
    // contactPerShadow_(),
    // scaleFaceTractionsNearDownstreamPatch_
    // (
    //     dict.lookupOrDefault<Switch>
    //     (
    //         "scaleFaceTractionsNearDownstreamPatch",
    //         Switch(false)
    //     )
    // ),
    // scaleTractionFieldPtr_(),
    curTimeIndex_(-1)
{
    if (debug)
    {
        Info<< "Creating " << solidContactPointPatchVectorField::typeName
            << " patch" << endl;
    }

    // Master creates contact laws
    if (master_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));

        // if (debug)
        // {
        //     Info<< "    writePointDistanceFields: " << writePointDistanceFields_
        //         << endl;
        // }

        // if (scaleFaceTractionsNearDownstreamPatch_)
        // {
        //     WarningIn(type() + "::" + type())
        //         << "scaleFaceTractionsNearDownstreamPatch can only be applied on"
        //         << "the slave patch: this option will be ignored for the master "
        //         << "patch!" << endl;
        // }
    }

    if (dict.found("value"))
    {
        solidTractionPointPatchVectorField::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        solidTractionPointPatchVectorField::operator==(vector::zero);
    }

    // Todo: use read if present for traction and pressure
    traction() = vector::zero;
    pressure() = 0.0;
}


solidContactPointPatchVectorField::solidContactPointPatchVectorField
(
    const solidContactPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    solidTractionPointPatchVectorField(p, iF),
    dict_(ptf.dict_),
    master_(ptf.master_),
    // writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowPatchNames_(ptf.shadowPatchNames_),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    // frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
#ifdef FOAMEXTEND
    writeZoneVTK_(ptf.writeZoneVTK_),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
#endif
    // contact_(ptf.contact_),
    // contactPerShadow_(),
    // scaleFaceTractionsNearDownstreamPatch_
    // (
    //     ptf.scaleFaceTractionsNearDownstreamPatch_
    // ),
    // scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects: they will be re-created.
}


#ifndef OPENFOAMFOUNDATION
solidContactPointPatchVectorField::solidContactPointPatchVectorField
(
    const solidContactPointPatchVectorField& ptf
)
:
    solidTractionPointPatchVectorField(ptf),
    dict_(ptf.dict_),
    master_(ptf.master_),
    // writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowPatchNames_(ptf.shadowPatchNames_),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    // frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
    #ifdef FOAMEXTEND
    writeZoneVTK_(ptf.writeZoneVTK_),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
    #endif
    // contact_(ptf.contact_),
    // contactPerShadow_(),
    // scaleFaceTractionsNearDownstreamPatch_
    // (
    //     ptf.scaleFaceTractionsNearDownstreamPatch_
    // ),
    // scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects
}
#endif


solidContactPointPatchVectorField::solidContactPointPatchVectorField
(
    const solidContactPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(ptf, iF),
    dict_(ptf.dict_),
    master_(ptf.master_),
    // writePointDistanceFields_(ptf.writePointDistanceFields_),
    shadowPatchNames_(ptf.shadowPatchNames_),
    shadowPatchIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    normalModels_(ptf.normalModels_),
    // frictionModels_(ptf.frictionModels_),
    normalPenaltyFactors_(ptf.normalPenaltyFactors_.size(), -1),
    zonePtr_(NULL),
    shadowZones_(),
    zoneToZones_(),
#ifdef FOAMEXTEND
    writeZoneVTK_(ptf.writeZoneVTK_),
    quickReject_(ptf.quickReject_),
    regionOfInterestTopCorner_(ptf.regionOfInterestTopCorner_),
    regionOfInterestBottomCorner_(ptf.regionOfInterestBottomCorner_),
    regionOfInterest_(ptf.regionOfInterest_),
#endif
    // contact_(ptf.contact_),
    // contactPerShadow_(),
    // scaleFaceTractionsNearDownstreamPatch_
    // (
    //     ptf.scaleFaceTractionsNearDownstreamPatch_
    // ),
    // scaleTractionFieldPtr_(),
    curTimeIndex_(ptf.curTimeIndex_)
{
    // Do not copy pointer objects
}

// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solidContactPointPatchVectorField::
~solidContactPointPatchVectorField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void solidContactPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
// #ifdef OPENFOAMFOUNDATION
//     m(contact_, contact_);
// #else
//     contact_.autoMap(m);
// #endif

//     scaleTractionFieldPtr_.clear();

//     if (contactPerShadow_.size())
//     {
//         forAll(contactPerShadow_, shadI)
//         {
// #ifdef OPENFOAMFOUNDATION
//             m(contactPerShadow_[shadI], contactPerShadow_[shadI]);
// #else
//             contactPerShadow_[shadI].autoMap(m);
// #endif
//         }
//     }

    // if (shadowPatchNames_.size() > 0)
    // {
    //     // Let the contact models know about the mapping
    //     // Be careful, we must pass slave
    //     // FIX PC 21-Sep-17: move this check inside if (shadowPatchNames ... )
    //     if (!master_)
    //     {
    //         normalModelForThisSlave().autoMap(m);
    //         // frictionModelForThisSlave().autoMap(m);
    //     }
    // }

    // Force all data to be re-created when needed
    clearOut();

    // Reset normal pelanty factors to reinitialise normal models
    normalPenaltyFactors_ =  -1;
}


// Grab the values using rmap
void solidContactPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    // solidTractionPointPatchVectorField::rmap(ptf, addr);

    // const solidContactPointPatchVectorField& tiptf =
    //   refCast<const solidContactPointPatchVectorField>(ptf);

    // // PC, I'm not sure if this "if" check should be here...
    // if (shadowPatchNames_.size() > 0)
    // {
    //     contact_.rmap(dmptf.contact_, addr);

    //     if (dmptf.contactPerShadow_.size())
    //     {
    //         // Force contactPerShadow to be initialised
    //         contactPerShadow();

    //         forAll(contactPerShadow_, shadI)
    //         {
    //             contactPerShadow_[shadI].rmap
    //             (
    //                 dmptf.contactPerShadow_[shadI], addr
    //             );
    //         }
    //     }
    // }

    // scaleTractionFieldPtr_.clear();
}


const Foam::wordList&
Foam::solidContactPointPatchVectorField::shadowPatchNames() const
{
    if (shadowPatchNames_.size() == 0)
    {
        makeShadowPatchNames(dict_);
    }

    return shadowPatchNames_;
}


const Foam::labelList&
Foam::solidContactPointPatchVectorField::shadowPatchIndices() const
{
    if (!shadowPatchIndicesPtr_)
    {
        calcShadowPatchIndices();
    }

    return *shadowPatchIndicesPtr_;
}


const Foam::solidContactPointPatchVectorField&
Foam::solidContactPointPatchVectorField::shadowPatchField() const
{
    if (shadowPatchIndices().size() != 1)
    {
        FatalErrorIn
        (
            "const Foam::solidContactPointPatchVectorField&\n"
            "Foam::solidContactPointPatchVectorField::shadowPatchField() const"
        )   << "This function can only be called for a patch with 1 shadow "
            << "patch; this patch has " << shadowPatchIndices().size()
            << " shadow patches!" << abort(FatalError);
    }

    const pointVectorField& field =
#ifdef OPENFOAMESIORFOUNDATION
        db().lookupObject<pointVectorField>(internalField().name());
#else
        db().lookupObject<pointVectorField>(dimensionedInternalField().name());
#endif

    return
        refCast<const solidContactPointPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[0]]
        );
}


Foam::PtrList<Foam::pointNormalContactModel>&
Foam::solidContactPointPatchVectorField::normalModels()
{
    if (master_)
    {
        if (normalModels_.size() == 0)
        {
            makeNormalModels(dict_);
        }

        return normalModels_;
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

        return shadowPatchField.normalModels();
    }
}


const Foam::PtrList<Foam::pointNormalContactModel>&
Foam::solidContactPointPatchVectorField::normalModels() const
{
    if (master_)
    {
        if (normalModels_.size() == 0)
        {
            makeNormalModels(dict_);
        }

        return normalModels_;
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

        return shadowPatchField.normalModels();
    }
}


// Foam::PtrList<Foam::frictionContactModel>&
// Foam::solidContactPointPatchVectorField::frictionModels()
// {
//     if (master_)
//     {
//         if (frictionModels_.size() == 0)
//         {
//             makeFrictionModels(dict_);
//         }

//         return frictionModels_;
//     }
//     else
//     {
//         const volVectorField& field =
//             db().lookupObject<volVectorField>
//             (
// #ifdef OPENFOAMESIORFOUNDATION
//                 internalField().name()
// #else
//                 dimensionedInternalField().name()
// #endif
//             );

//         solidContactPointPatchVectorField& shadowPatchField =
//             const_cast<solidContactPointPatchVectorField&>
//             (
//                 refCast<const solidContactPointPatchVectorField>
//                 (
//                     field.boundaryField()[shadowPatchIndices()[0]]
//                 )
//             );

//         return shadowPatchField.frictionModels();
//     }
// }


// const Foam::PtrList<Foam::frictionContactModel>&
// Foam::solidContactPointPatchVectorField::frictionModels() const
// {
//     if (master_)
//     {
//         if (frictionModels_.size() == 0)
//         {
//             makeFrictionModels(dict_);
//         }

//         return frictionModels_;
//     }
//     else
//     {
//         const volVectorField& field =
//             db().lookupObject<volVectorField>
//             (
// #ifdef OPENFOAMESIORFOUNDATION
//                 internalField().name()
// #else
//                 dimensionedInternalField().name()
// #endif
//             );

//         solidContactPointPatchVectorField& shadowPatchField =
//             const_cast<solidContactPointPatchVectorField&>
//             (
//                 refCast<const solidContactPointPatchVectorField>
//                 (
//                     field.boundaryField()[shadowPatchIndices()[0]]
//                 )
//             );

//         return shadowPatchField.frictionModels();
//     }
// }


Foam::pointNormalContactModel&
Foam::solidContactPointPatchVectorField::normalModelForThisSlave()
{
    if (master_)
    {
        FatalErrorIn
        (
            "pointNormalContactModel&"
            "solidContactPointPatchVectorField::normalModelForThisSlave()"
        )   << "The master is not allowed to called this fucntion!"
            << abort(FatalError);
    }

    // Lookup the master patch corresponding to the current slave patch
    const pointVectorField& field =
        db().lookupObject<pointVectorField>
        (
#ifdef OPENFOAMESIORFOUNDATION
            internalField().name()
#else
            dimensionedInternalField().name()
#endif
        );

    if (returnReduce(shadowPatchIndices().size() == 0, maxOp<bool>()))
    {
        FatalError
            << "shadowPatchIndices().size() == 0" << exit(FatalError);
    }

    const solidContactPointPatchVectorField& masterPatchField =
        refCast<const solidContactPointPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[0]]
        );

    // The master may have multiple slaves so we need to find which model
    // corresponds to the current slave patch
    const wordList& shadPatchNames = masterPatchField.shadowPatchNames();
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
            "void solidContactPointPatchVectorField::"
            "normalModelForThisSlave()"
        )   << "Something went wrong when looking for the shadowPatch"
            << abort(FatalError);
    }

    return normalModels()[masterShadowID];
}


// Foam::frictionContactModel&
// Foam::solidContactPointPatchVectorField::frictionModelForThisSlave()
// {
//     if (master_)
//     {
//         FatalErrorIn
//         (
//             "frictionContactModel& "
//             "solidContactPointPatchVectorField::frictionModelForThisSlave()"
//         )   << "The master is not allowed to called this function!"
//             << abort(FatalError);
//     }

//     // Lookup the master patch corresponding to the current slave patch
//     const volVectorField& field =
//         db().lookupObject<volVectorField>
//         (
// #ifdef OPENFOAMESIORFOUNDATION
//             internalField().name()
// #else
//             dimensionedInternalField().name()
// #endif
//         );

//     const solidContactPointPatchVectorField& masterPatchField =
//         refCast<const solidContactPointPatchVectorField>
//         (
//             field.boundaryField()[shadowPatchIndices()[0]]
//         );

//     // The master may have multiple slaves so we need to find which model
//     // corresponds to the current slave patch
//     const wordList& shadPatchNames = masterPatchField.shadowPatchNames();
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
//             "void solidContactPointPatchVectorField::"
//             "frictionModelForThisSlave()"
//         )   << "Something went wrong when looking for the shadowPatch"
//             << abort(FatalError);
//     }

//     return frictionModels()[masterShadowID];
// }


void solidContactPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = this->db().time().timeIndex();

//         if (master_)
//         {
//             // Let the contact models know that it is a new time-step, in case
//             // they need to update anything
//             forAll(shadowPatchNames(), shadPatchI)
//             {
//                 normalModels()[shadPatchI].newTimeStep();
//                 frictionModels()[shadPatchI].newTimeStep();

// #ifdef FOAMEXTEND
//                 // Force N^2 contact search at least once per time-step
//                 zoneToZones()[shadPatchI].clearPrevCandidateMasterNeighbors();
// #endif
//             }
//         }
    }

    // Move the master and slave zone to the deformed configuration
    moveZonesToDeformedConfiguration();

    // Delete the zone-to-zone interpolator weights as the zones have moved
    const wordList& shadPatchNames = shadowPatchNames();
    forAll(shadPatchNames, shadPatchI)
    {
#ifdef OPENFOAMESIORFOUNDATION
        // Clear the interpolators each time
        // Ther is no need to do this in the loop but we will keep it here for
        // tidiness
        zoneToZones_.clear();
#else
        zoneToZones()[shadPatchI].movePoints
        (
            tensorField(0), tensorField(0), vectorField(0)
        );
#endif
    }

    // Calculate and apply contact forces
    if (master_)
    {
        // Reset the traction to zero as we will accumulate it over all the
        // shadow patches
        traction() = vector::zero;

        forAll(shadPatchNames, shadPatchI)
        {
            // Calculate the slave patch point unit normals as they are used by
            // both the normal and friction models
            const vectorField shadowPatchPointNormals
            (
                shadowZones()[shadPatchI].globalPointToPatch
                (
                    shadowZones()[shadPatchI].globalPatch().pointNormals()
                )
            );

            // Interpolate the master displacement increment to the slave patch
            // as it is required by specific normal and friction contact models

            // vectorField patchDD(patch().size(), vector::zero);
            // vectorField shadowPatchDD
            // (
            //     patch().boundaryMesh()[shadowPatchIndices()[shadPatchI]].size(),
            //     vector::zero
            // );

            // We will lookup the total displacement and old total
            // displacement

            // const pointVectorField& pointD =
            //     db().lookupObject<pointVectorField>("pointD");

            // patchDD =
            //     pointD.boundaryField()[patch().index()]
            //   - pointD.oldTime().boundaryField()[patch().index()];
            // shadowPatchDD =
            //     pointD.boundaryField()[shadowPatchIndices()[shadPatchI]]
            //   - pointD.oldTime().boundaryField()
            //     [
            //         shadowPatchIndices()[shadPatchI]
            //     ];

            // Master zone DD
            // const vectorField zoneDD(zone().patchPointToGlobal(patchDD));

//             // Master patch DD interpolated to the slave patch
//             const vectorField patchDDInterpToShadowPatch
//              (
//                 shadowZones()[shadPatchI].globalPointToPatch
//                 (
// #ifdef OPENFOAMESIORFOUNDATION
//                     fixThis() //zoneToZones()[shadPatchI].interpolateToTarget(zoneDD)()
//                     // this is for faces: we need point-to-point interpolate
// #else
//                     fixThis() //zoneToZones()[shadPatchI].masterToSlave(zoneDD)()
// #endif
//                 )
//              );

            // Calculate normal contact forces
            // shadowPatchDD is the DD on the shadow patch, whereas
            // patchDDInterpToShadowPatch is the master patch DD interpolated to
            // the shadow; and the difference between these two is the slip (and
            // also the normal component of DD)
            normalModels()[shadPatchI].correct
            (
                shadowPatchPointNormals,
                shadowZones()[shadPatchI].globalPointToPatch
                (
#ifdef OPENFOAMESIORFOUNDATION
                    zoneToZones()[shadPatchI].targetPointDistanceToIntersection()
#else
                    zoneToZones()[shadPatchI].slavePointDistanceToIntersection()
#endif
                ) //,
                // zoneToZones()[shadPatchI],
                //shadowPatchDD,
                //patchDDInterpToShadowPatch
            );

            // Calculate friction contact forces: to be implemented
            // frictionModels()[shadPatchI].correct
            // (
            //     normalModels()[shadPatchI].slavePressure(),
            //     shadowPatchFaceNormals,
            //     normalModels()[shadPatchI].areaInContact(),
            //     shadowPatchDD,
            //     patchDDInterpToShadowPatch
            // );

            if (rigidMaster_)
            {
                // Set to master to traction free to mimic a rigid contact
                traction() = vector::zero;

                // Set contact indicator field
                // contactPerShadow()[shadPatchI] = 0.0;
            }
            else
            {
                // Interpolate slave traction to the master
                const vectorField slavePatchTraction
                (
                   // - frictionModels()[shadPatchI].slaveTractionForMaster()
                   - normalModels()[shadPatchI].slavePressure()
                );

                const vectorField slaveZoneTraction
                (
                    shadowZones()[shadPatchI].patchPointToGlobal
                    (
                        slavePatchTraction
                    )
                );

                // Calculate traction for this contact
                vectorField tractionForThisShadow
                (
                    zone().globalPointToPatch
                    (
#ifdef OPENFOAMESIORFOUNDATION
                        zoneToZones()[shadPatchI].interpolateToSourcePoints
#else
                        fixThis() //zoneToZones()[shadPatchI].slaveToMaster
#endif
                        (
                            slaveZoneTraction
                        )()
                    )
                );

                // Accumulate the traction on the master patch
                traction() += tractionForThisShadow;

                // // Update contactPerShadow field
                // // Note: this is used by thermalContact to know which faces
                // // are in contact
                // const scalarField magTraction(mag(tractionForThisShadow));
                // const scalar tol = 1e-6*gMax(magTraction);
                // scalarField& contactForThisShadow =
                //     contactPerShadow()[shadPatchI];
                // forAll(contactForThisShadow, faceI)
                // {
                //     if (magTraction[faceI] > tol)
                //     {
                //         contactForThisShadow[faceI] = 1.0;
                //     }
                //     else
                //     {
                //         contactForThisShadow[faceI] = 0.0;
                //     }
                // }
            }
        }
    }
    else
    {
        // Set the traction on the slave patch
        // The master stores the friction and normal models, so we need to find
        // which models correspond to the current shadow
        traction() =
            normalModelForThisSlave().slavePressure();
          //   frictionModelForThisSlave().slaveTraction()
          // + normalModelForThisSlave().slavePressure();

        // // TESTING - START
        // // Scale traction vectors on faces, which share an edge with the
        // // downstream patch
        // // This is an attempt to fix an issue where the first row of faces
        // // deform unphysically when being drawn into the die
        // if (scaleFaceTractionsNearDownstreamPatch_)
        // {
        //     traction() *= scaleTractionField();
        // }
        // // TESTING - END

        // // Update contactPerShadow field
        // // Note: this is used by thermalContact to know which faces
        // // are in contact
        // const scalarField magTraction(mag(traction()));
        // const scalar tol = 1e-6*gMax(magTraction);
        // scalarField& contactForThisShadow = contactPerShadow()[0];
        // forAll(contactForThisShadow, faceI)
        // {
        //     if (magTraction[faceI] > tol)
        //     {
        //         contactForThisShadow[faceI] = 1.0;
        //     }
        //     else
        //     {
        //         contactForThisShadow[faceI] = 0.0;
        //     }
        // }
    }

    // Accumulate the contact indicator field
    // contact_ = 0.0;
    // PtrList<scalarField>& contactPerShadow = this->contactPerShadow();
    // forAll(contactPerShadow, shadI)
    // {
    //     contact_ += contactPerShadow[shadI];
    // }

    // // Scale any face in contact with more than one shadow
    // if (gMax(contact_) > (1.0 + SMALL))
    // {
    //     forAll(contact_, faceI)
    //     {
    //         if (contact_[faceI] > (1.0 + SMALL))
    //         {
    //             // Update the contact weights corresponding to each shadow
    //             scalar sumContact = 0.0;
    //             forAll(contactPerShadow, shadI)
    //             {
    //                 contactPerShadow[shadI][faceI] /= contact_[faceI];
    //                 sumContact += contactPerShadow[shadI][faceI];
    //             }

    //             if (sumContact > (1.0 + SMALL))
    //             {
    //                 FatalErrorIn
    //                 (
    //                     "void solidContactPointPatchVectorField::"
    //                     "updateCoeffs()"
    //                 )   << "There is a problem normalising the contact field"
    //                     << ", sumContact is: " << sumContact
    //                     << abort(FatalError);
    //             }

    //             // Reset accumulated contact face value to 1.0
    //             contact_[faceI] = 1.0;
    //         }
    //     }
    // }

    solidTractionPointPatchVectorField::initEvaluate(commsType);
}


// Write
void solidContactPointPatchVectorField::write(Ostream& os) const
{
    // If the shadowPatchIndices pointer is not set then we will assume that the
    // contact models were not created and nothing has changed; so we will just
    // output the input dict unchanged
    if (shadowPatchNames_.size() == 0)
    {
        // Overwrite fields in the dict
        dictionary& dict = const_cast<dictionary&>(dict_);

        // dict.remove("gradient");
        dict.remove("value");
        dict.remove("traction");
        dict.remove("pressure");

        //dict.add("gradient", gradient());
        const vectorField patchValue(patchInternalField());
        dict.add("value", patchValue);
        //dict.add("traction", traction());
        //dict.add("pressure", pressure());

        // Write the dictionary
        dict_.write(os, false);

#ifdef OPENFOAMFOUNDATION
        // writeEntry(os, "gradient", gradient());
        writeEntry(os, "value", patchValue);
        writeEntry(os, "traction", traction());
        writeEntry(os, "pressure", pressure());
#else
        // gradient().writeEntry("gradient", os);
        patchValue.writeEntry("value", os);
        traction().writeEntry("traction", os);
        pressure().writeEntry("pressure", os);
#endif

        return;
    }

    solidTractionPointPatchVectorField::write(os);

    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;

    const wordList& shadPatchNames = shadowPatchNames();
    if (shadPatchNames.size() == 1)
    {
        os.writeKeyword("shadowPatch")
            << shadPatchNames[0] << token::END_STATEMENT << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "shadowPatches", shadowPatchNames());
#else
        shadowPatchNames().writeEntry("shadowPatches", os);
#endif
    }

// #ifdef FOAMEXTEND
//     os.writeKeyword("regionOfInterest")
//         << regionOfInterest_ << token::END_STATEMENT << nl;
//     os.writeKeyword("regionOfInterestTopCorner")
//         << regionOfInterestTopCorner_ << token::END_STATEMENT << nl;
//     os.writeKeyword("regionOfInterestBottomCorner")
//         << regionOfInterestBottomCorner_ << token::END_STATEMENT << nl;
//     os.writeKeyword("writeZoneVTK")
//         << writeZoneVTK_ << token::END_STATEMENT << nl;
// #endif

    // os.writeKeyword("writePointDistanceFields")
    //     << writePointDistanceFields_ << token::END_STATEMENT << nl;
    // os.writeKeyword("scaleFaceTractionsNearDownstreamPatch")
    //     << scaleFaceTractionsNearDownstreamPatch_ << token::END_STATEMENT << nl;
    // if (scaleFaceTractionsNearDownstreamPatch_)
    // {
    //     os.writeKeyword("downstreamScaleFactor")
    //         << readScalar(dict_.lookup("downstreamScaleFactor"))
    //         << token::END_STATEMENT << nl;
    //     os.writeKeyword("downstreamPatchName")
    //         << word(dict_.lookup("downstreamPatchName")) << token::END_STATEMENT
    //         << nl;
    // }

    if (master_)
    {
        os.writeKeyword("rigidMaster") << rigidMaster_
            << token::END_STATEMENT << nl;

        if (shadowPatchNames_.size() == 1)
        {
            os.writeKeyword("normalContactModel")
                << normalModels()[0].type()
                << token::END_STATEMENT << nl;
            normalModels()[0].writeDict(os);

            // os.writeKeyword("frictionContactModel")
            //     << frictionModels()[0].type()
            //     << token::END_STATEMENT << nl;
            // frictionModels()[0].writeDict(os);

            // os.writeKeyword("useNewPointDistanceMethod")
            //     << dict_.lookupOrDefault<Switch>
            //     (
            //         "useNewPointDistanceMethod", false
            //     )
            //     << token::END_STATEMENT << nl;

            // os.writeKeyword("projectPointsToPatchBoundary")
            //     << dict_.lookupOrDefault<Switch>
            //     (
            //         "projectPointsToPatchBoundary", false
            //     )
            //     << token::END_STATEMENT << nl;

            // os.writeKeyword("checkPointDistanceOrientations")
            //     << dict_.lookupOrDefault<Switch>
            //     (
            //         "checkPointDistanceOrientations", false
            //     )
            //     << token::END_STATEMENT << nl;

            // os.writeKeyword("usePrevCandidateMasterNeighbors")
            //     << dict_.lookupOrDefault<Switch>
            //     (
            //         "usePrevCandidateMasterNeighbors", false
            //     )
            //     << token::END_STATEMENT << nl;
        }
        else
        {
            forAll(shadowPatchNames_, shadPatchI)
            {
                os  << patch().name() << "_to_"
                    << shadowPatchNames_[shadPatchI] << "_dict" << nl
                    << '{' << endl;

                os.writeKeyword("normalContactModel")
                    << normalModels()[shadPatchI].type()
                    << token::END_STATEMENT << nl;
                normalModels()[shadPatchI].writeDict(os);

                // os.writeKeyword("frictionContactModel")
                //     << frictionModels()[shadPatchI].type()
                //     << token::END_STATEMENT << nl;
                // frictionModels()[shadPatchI].writeDict(os);

                os  << '}' << endl;
            }
        }
    }

// #ifdef FOAMEXTEND
//     if (writeZoneVTK_)
//     {
//         if
//         (
//             dimensionedInternalField().name() == "D"
//          || dimensionedInternalField().name() == "DD"
//         )
//         {
//             Info<< "Writing deformed zones to VTK" << endl;
//             const word timeName =
//                 patch().boundaryMesh().mesh().time().timeName();

//             zone().globalPatch().writeVTK("zone_" + timeName);

//             forAll(shadowZones(), shadI)
//             {
//                 shadowZones()[shadI].globalPatch().writeVTK
//                 (
//                     "shadowZone_" + timeName
//                 );
//             }
//         }
//     }
// #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    solidContactPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
