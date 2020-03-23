/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    solidRigidContactFvPatchVectorField

Description
    Contact boundary conditions for finite volume solidMechanics solvers.

    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the normalContactModel and
    frictionContactModel.

    The boundary condition keeps a copy of the master and slave patch and
    keeps them in the deformed configuration; this allows the calculation of
    contact distances and interpolation between the globalPatchZones.

    The distance calculations and interpolations are performed by the GGI class.

    More details in:

    P. Cardiff, A. Karać, A. Ivanković: Development of a Finite Volume contact
    solver based on the penalty method. Computational Materials Science,
    03/2014, 64:283–284. DOI: 10.1016/j.commatsci.2012.03.011.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

    P. Cardiff, Development of the Finite Volume Method for Hip Joint Stress
    Analysis, PhD thesis, University College Dublin, 2012.


SourceFiles
    solidRigidContactFvPatchVectorField.C
    solidRigidContactFvPatchVectorFieldCalcContact.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Rewrite by Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#ifndef solidRigidContactFvPatchVectorField_H
#define solidRigidContactFvPatchVectorField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"
#include "globalPolyPatch.H"
#include "newGgiInterpolation.H"
#include "Switch.H"
#include "boundBox.H"
#include "standAlonePatch.H"
#include "triSurfaceMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class solidRigidContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidRigidContactFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    // Private data

        //- Store a copy of the patch dictionary
        //  This allows us to use lazy evaluation when creating the contact
        //  models
        const dictionary dict_;

        //- Write the deformed contact zones to VTK surface files
        const Switch writeZoneVTK_;

        //- Write the contact point distances fields for viewing
        const Switch writePointDistanceFields_;

        //- Shadow triSurface names
        mutable wordList shadowTriSurfNames_;

        //- Normal contact model pointers
        mutable PtrList<normalContactModel> normalModels_;

        //- Friction contact model pointer
        mutable PtrList<frictionContactModel> frictionModels_;

        //- Store current penalty scale factors (if using penalty contact model)
        mutable scalarField normalPenaltyFactors_;

        //- Master face zone kept in deformed configuration
        mutable globalPolyPatch* zonePtr_;

        //- Shadow face zones kept in deformed configuration
        //  Stored by the master patch
        //  These are STL surfaces
        mutable PtrList<standAlonePatch> shadowZones_;

        //- List of motions for the shadow zones
        mutable PtrList<triSurfaceMotion> shadowZonesMotion_;

        //- Zone-to-zone interpolations
        mutable PtrList<newGgiStandAlonePatchInterpolation> zoneToZones_;

        //- GGI quick reject algorithm
        newGgiInterpolation::quickReject quickReject_;

        //- Top corner of the region where te contact is expected to be
        // This aims to save time in the GGI distance calculation
        const vector regionOfInterestTopCorner_;

        //- Bottom corner of the region where te contact is expected to be
        // This aims to save time in the GGI distance calculation
        const vector regionOfInterestBottomCorner_;

        //- Region where the contact is expected to be
        //  This aims to save time in the GGI distance calculations
        const boundBox regionOfInterest_;

        //- Contact field for all the contacts, stored on the current patch
        //  1 means in contact
        //  0 means not in contact
        scalarField contact_;

        //- Contact field corresponding to each shadow patch, stored on the
        //  current patch
        //  1 means in contact
        //  0 means not in contact
        mutable PtrList<scalarField> contactPerShadow_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Make the shadow triSurface names
        void makeShadowTriSurfNames(const dictionary& dict) const;

        //- Calculate the shadow patch indices
        // void calcShadowPatchIndices() const;

        //- Make normal contact models
        void makeNormalModels(const dictionary& dict) const;

        //- Make friction contact models
        void makeFrictionModels(const dictionary& dict) const;

        //- Clear out demand driven data
        virtual void clearOut();

        //- Create the zone
        void calcZone() const;

        //- Create the shadow zones
        void calcShadowZones() const;

        //- Create the shadow zones motion
        void calcShadowZonesMotion() const;

        //- Create the zone-to-zone interpolators
        void calcZoneToZones() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();

        //- Move the rigid triSurfaces to their deformed position
        void moveTriSurfaces();

        // Set the contactPerShadow field
        void calcContactPerShadow() const;


public:

    //- Runtime type information
    TypeName("solidRigidContact");


    // Constructors

        //- Construct from patch and internal field
        solidRigidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidRigidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidRigidContactFvPatchVectorField onto
        //  a new patch
        solidRigidContactFvPatchVectorField
        (
            const solidRigidContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidRigidContactFvPatchVectorField
        (
            const solidRigidContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidRigidContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidRigidContactFvPatchVectorField
        (
            const solidRigidContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidRigidContactFvPatchVectorField(*this, iF)
            );
        }

    //- Destructor
    virtual ~solidRigidContactFvPatchVectorField();


    // Member functions

        // Access

            //- Return the shadow triSurface names
            const wordList& shadowTriSurfNames() const;

            //- Return reference to the normal contact models
            PtrList<normalContactModel>& normalModels();

            //- Return const reference to the normal contact models
            const PtrList<normalContactModel>& normalModels() const;

            //- Return reference to the friction contact models
            PtrList<frictionContactModel>& frictionModels();

            //- Return const reference to the friction contact models
            const PtrList<frictionContactModel>& frictionModels() const;

            //- Return const reference to deformed zone
            const globalPolyPatch& zone() const;

            //- Return reference to deformed zone
            globalPolyPatch& zone();

            //- Return const reference to deformed shadow zones
            const PtrList<standAlonePatch>& shadowZones() const;

            // //- Return reference to deformed shadow zones
            PtrList<standAlonePatch>& shadowZones();

            //- Return const reference to deformed shadow zones motion
            const PtrList<triSurfaceMotion>& shadowZonesMotion() const;

            // //- Return reference to deformed shadow zones motion
            PtrList<triSurfaceMotion>& shadowZonesMotion();

            //- Return const reference to patch-to-patch interpolation
            const PtrList<newGgiStandAlonePatchInterpolation>&
                zoneToZones() const;

            //- Return reference to patch-to-patch interpolation
            PtrList<newGgiStandAlonePatchInterpolation>& zoneToZones();

            //- Return reference to contact field
            const scalarField& contact() const
            {
                return contact_;
            }

            //- Return reference to contactPerShadow fields
            PtrList<scalarField>& contactPerShadow();

            //- Return const reference to contactPerShadow fields
            const PtrList<scalarField>& contactPerShadow() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // Member functions

            //- Return the rate of dissipated energy due to friction for the
            //  master patch
            virtual tmp<scalarField> frictionHeatRate() const;

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#ifdef NoRepository
//#   include "solidRigidContactFvPatchVectorFieldTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
