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

Class
    solidContactPointPatchVectorField

Description
    The point displacement version of solidContactFvPatchField.

    Used to enforce a traction condition for point based solid mechanics
    solvers.

SourceFiles
    solidContactPointPatchVectorField.C
    solidContactPointPatchVectorFieldCalcContact.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef solidContactPointPatchVectorField_H
#define solidContactPointPatchVectorField_H

#include "solidTractionPointPatchVectorField.H"
#include "pointNormalContactModel.H"
// #include "frictionContactModel.H"
#include "globalPolyPatch.H"
#include "Switch.H"
#include "boundBox.H"
#include "standAlonePatch.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "amiZoneInterpolation.H"
#else
    #include "newGgiInterpolation.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// #ifdef OPENFOAMESIORFOUNDATION
//     typedef pointPatchFieldMapper PointPatchFieldMapper;
// #endif

/*---------------------------------------------------------------------------*\
            Class solidContactPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidContactPointPatchVectorField
:
    public solidTractionPointPatchVectorField
{

    // Private data

        //- Store a copy of the patch dictionary
        //  This allows us to use lazy evaluation when creating the contact
        //  models
        const dictionary dict_;

        //- Is this patch the master
        const Switch master_;

        // //- Write the contact point distances fields for viewing
        // const Switch writePointDistanceFields_;

        //- Shadow patch names
        mutable wordList shadowPatchNames_;

        //- Shadow patch indices
        mutable labelList* shadowPatchIndicesPtr_;

        //- Is the contact boundary active
        Switch rigidMaster_;

        //- Normal contact model pointers
        mutable PtrList<pointNormalContactModel> normalModels_;

        // //- Friction contact model pointer
        // mutable PtrList<frictionContactModel> frictionModels_;

        //- Store current penalty scale factors (if using penalty contact model)
        mutable scalarField normalPenaltyFactors_;

        //- Master face zone kept in deformed configuration
        mutable globalPolyPatch* zonePtr_;

        //- Shadow face zones kept in deformed configuration
        //  Stored by the master patch
        mutable PtrList<globalPolyPatch> shadowZones_;

        //- Zone-to-zone interpolations
// #ifdef OPENFOAMESIORFOUNDATION
        mutable PtrList<amiZoneInterpolation> zoneToZones_;
// #else
//         mutable PtrList<newGgiStandAlonePatchInterpolation> zoneToZones_;

//         //- Write the deformed contact zones to VTK surface files
//         const Switch writeZoneVTK_;

//         //- GGI quick reject algorithm
//         newGgiInterpolation::quickReject quickReject_;

//         //- Top corner of the region where te contact is expected to be
//         // This aims to save time in the GGI distance calculation
//         const vector regionOfInterestTopCorner_;

//         //- Bottom corner of the region where te contact is expected to be
//         // This aims to save time in the GGI distance calculation
//         const vector regionOfInterestBottomCorner_;

//         //- Region where the contact is expected to be
//         //  This aims to save time in the GGI distance calculations
//         const boundBox regionOfInterest_;
// #endif

        // //- Contact field for all the contacts, stored on the current patch
        // //  1 means in contact
        // //  0 means not in contact
        // scalarField contact_;

        // //- Contact field corresponding to each shadow patch, stored on the
        // //  current patch
        // //  1 means in contact
        // //  0 means not in contact
        // mutable PtrList<scalarField> contactPerShadow_;

        // //- Flag to enable approach for scaling traction near the downstream
        // //  patch to fix an issue with unphysical behaviour
        // Switch scaleFaceTractionsNearDownstreamPatch_;

        // //- Field to scale traction of slave patch
        // mutable autoPtr<scalarField> scaleTractionFieldPtr_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Make the shadow patch names
        void makeShadowPatchNames(const dictionary& dict) const;

        //- Calculate the shadow patch indices
        void calcShadowPatchIndices() const;

        //- Make normal contact models
        void makeNormalModels(const dictionary& dict) const;

        // //- Make friction contact models
        // void makeFrictionModels(const dictionary& dict) const;

        // //- Clear out demand driven data
        // //virtual void clearOut();

        //- Create the zone
        void calcZone() const;

        //- Create the shadow zones
        void calcShadowZones() const;

        //- Create the zone-to-zone interpolators
        void calcZoneToZones() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();

        // // Set the contactPerShadow field
        // void calcContactPerShadow() const;

        // //- Return field to scale the traction field on the slave patch
        // scalarField scaleTractionField() const;

        // //- Make the scaleTractionField
        // void makeScaleTractionField() const;


public:

    //- Runtime type information
    TypeName("pointSolidContact");


    // Constructors

        //- Construct from patch and internal field
        solidContactPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidContactPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        solidContactPointPatchVectorField
        (
            const solidContactPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

#ifndef OPENFOAMFOUNDATION
        //- Construct as copy
        solidContactPointPatchVectorField
        (
            const solidContactPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new solidContactPointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidContactPointPatchVectorField
        (
            const solidContactPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new solidContactPointPatchVectorField(*this, iF)
            );
        }

    //- Destructor
    virtual ~solidContactPointPatchVectorField();


    // Member functions

        // Access

            //- Return if this patch is the slave
            Switch master() const
            {
                return master_;
            }

            //- Return the shadow patch names
            const wordList& shadowPatchNames() const;

            //- Return the shadow patch indices
            const labelList& shadowPatchIndices() const;

            //- Return the shadow patch
            //  This function is not allowed if shadowPatchIndices.size() > 1
            const solidContactPointPatchVectorField& shadowPatchField() const;

            //- Return reference to the normal contact models
            PtrList<pointNormalContactModel>& normalModels();

            //- Return const reference to the normal contact models
            const PtrList<pointNormalContactModel>& normalModels() const;

            // //- Return reference to the friction contact models
            // PtrList<frictionContactModel>& frictionModels();

            // //- Return const reference to the friction contact models
            // const PtrList<frictionContactModel>& frictionModels() const;

            //- Return normal contact model for the slave patch calling this
            //  function
            pointNormalContactModel& normalModelForThisSlave();

            // //- Return friction contact model for the slave patch calling this
            // //  function
            // frictionContactModel& frictionModelForThisSlave();

            //- Return const reference to deformed zone
            const globalPolyPatch& zone() const;

            //- Return reference to deformed zone
            globalPolyPatch& zone();

            //- Return const reference to deformed shadow zones
            const PtrList<globalPolyPatch>& shadowZones() const;

            //- Return reference to deformed shadow zones
            PtrList<globalPolyPatch>& shadowZones();

#ifdef OPENFOAMESIORFOUNDATION
            //- Return const reference to patch-to-patch interpolation
            const PtrList<amiZoneInterpolation>& zoneToZones() const;

            //- Return reference to patch-to-patch interpolation
            PtrList<amiZoneInterpolation>& zoneToZones();
#else
            fixThis();

            // //- Return const reference to patch-to-patch interpolation
            // const PtrList<newGgiStandAlonePatchInterpolation>&
            //     zoneToZones() const;

            // //- Return reference to patch-to-patch interpolation
            // PtrList<newGgiStandAlonePatchInterpolation>& zoneToZones();
#endif

            //- Return const reference to zone-to-zone interpolation for the
            //  current slave patch; an error is thrown if the master calls this
            //  function
#ifdef OPENFOAMESIORFOUNDATION
            const amiZoneInterpolation& zoneToZoneForThisSlave() const;
#else
            const newGgiStandAlonePatchInterpolation&
            zoneToZoneForThisSlave() const;
#endif

            //- Return the global patch for the current slave patch
            //  An error is thrown if the master calls this function
            const globalPolyPatch& zoneForThisSlave() const;

            // //- Return reference to contact field
            // const scalarField& contact() const
            // {
            //     return contact_;
            // }

            // //- Return reference to contactPerShadow fields
            // PtrList<scalarField>& contactPerShadow();

            // //- Return const reference to contactPerShadow fields
            // const PtrList<scalarField>& contactPerShadow() const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto
            //  this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        //- Insert boundary value into the internal field
        //  Fixed values need to be inserted first. In the second go,
        //  there's nothing to do
        virtual void initEvaluate
        (
#ifdef OPENFOAMESIORFOUNDATION
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );

        //- Clear out demand driven data
        virtual void clearOut();

        // //- Return the rate of dissipated energy due to friction for the
        // //  master patch
        // virtual tmp<scalarField> frictionHeatRate() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
