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

Class
    solidContactFvPatchVectorField

Description
    Contact boundary conditions for finite volume solidMechanics solvers.

    There is a master and a slave patch. The master creates the contactModel
    and the slave gets the pointer to it.

    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the normalContactModel and
    frictionContactModel.

    More details in:

    P. Cardiff, A. Karać, A. Ivanković: Development of a Finite Volume contact
    solver based on the penalty method. Computational Materials Science,
    03/2014, 64:283–284. DOI: 10.1016/j.commatsci.2012.03.011.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

SourceFiles
    solidContactFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef solidContactFvPatchVectorField_H
#define solidContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"
#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"
#include "extendedGgiInterpolation.H"
#include "ExtendedGGIInterpolationTemplate.H"
#include "Switch.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class solidContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidContactFvPatchVectorField
:
        public solidTractionFvPatchVectorField
{
    // Private data

        //- Name of the displacement field
        const word fieldName_;

        //- is this patch the master
        const Switch master_;

        //- Shadow patch name
        word shadowPatchName_;

        //- Shadow patch index
        label shadowPatchIndex_;

        //- Zone index
        mutable label* zoneIndexPtr_;

        //- Shadow zone index
        mutable label* shadowZoneIndexPtr_;

        // is the contact boundary active
        Switch rigidMaster_;

        // Copy of contact dictionary
        const dictionary dict_;

        //- Normal contact model pointer
        mutable normalContactModel* normalModelPtr_;

        //- Friction contact model pointer
        mutable frictionContactModel* frictionModelPtr_;

        // Mesh.allPoints() where the contact zone points are kept in the
        // deformed configuration
        // PC: should no longer be needed as our zones will store their own
        // points
        //mutable vectorField* allPointsDeformedPtr_;

        // Zone points at end of previous time step for incremental solver,
        // zone points from initial configuration for total displacement solver
        // NOT needed as the mesh already has these!
        //mutable vectorField* oldZonePointsPtr_;

        // Shadow zone points at end of previous time step for incremental
        // solver,
        // shadow zone points from initial configuration for total displacement
        // solver
        // NOT needed as the mesh already has these!
        //mutable vectorField* oldShadowZonePointsPtr_;

        // Master face zone kept in deformed configuration
        mutable primitiveFacePatch* zonePtr_;

        // Slave face zone kept in deformed configuration
        mutable primitiveFacePatch* shadowZonePtr_;

        //- Zone-to-zone interpolation
        //PrimitivePatch<face, List, const pointField&>
        mutable extendedGgiZoneInterpolation* zoneToZonePtr_;

        //- GGI quick reject algorithm
        extendedGgiInterpolation::quickReject quickReject_;

        //- Region where the contact is expected to be. This aims to save
        // time in the GGI distance calculations
        const boundBox regionOfInterest_;

        //- Current time index
        label curTimeIndex_;

        //- Global face zones in the mesh
        mutable labelList* globalFaceZonesPtr_;

        // Point map from global to local for global face zones
        mutable IOList<labelList>* globalToLocalFaceZonePointMapPtr_;


    // Private Member Functions

            //- Move the master and slave zones to the deformed configuration
            void moveZonesToDeformedConfiguration();

            //- Calculate the zone index
            void calcZoneIndex() const;

            //- Calculate the shadow zone index
            void calcShadowZoneIndex() const;

            //- Create the normal contact model
            void calcNormalModel() const;

            //- Create the friction contact model
            void calcFrictionModel() const;

            //- Calculate zone
            void calcZone() const;

            //- Calculate shadow zone
            void calcShadowZone() const;

            //- Calculate interpolation
            void calcZoneToZone() const;

            //- Calculate slaveFaceNormals
    //void calcSlaveFaceNormals() const;

            //- Calculate global face zones
            void calcGlobalFaceZones() const;

            //- Calculate global to local face zone point map
            void calcGlobalToLocalFaceZonePointMap() const;

            // Contact dictionary
            const dictionary& dict() const
            {
                return dict_;
            };

            //- Get zone field from patch field
            template<class Type>
            tmp<Field<Type> > zoneField
            (
                const label zoneIndex,
                const label patchIndex,
                const Field<Type>& patchField
            ) const;

            //- Get patch field from zone field
            template<class Type>
            tmp<Field<Type> > patchField
            (
                const label patchIndex,
                const label zoneIndex,
                const Field<Type>& zoneField
            ) const;

            //- Interpolate zone face values to points
            template<class Type>
            tmp<Field<Type> > zoneFaceToPointInterpolate
            (
                const label zoneIndex,
                const Field<Type>& zoneField
            ) const;

            //- Interpolate zone point values to faces
            template<class Type>
            tmp<Field<Type> > zonePointToFaceInterpolate
            (
                const label zoneIndex,
                const Field<Type>& zonePointField
            ) const;

public:

    //- Runtime type information
    TypeName("solidContact");


    // Constructors

        //- Construct from patch and internal field
        solidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidContactFvPatchVectorField onto
        //  a new patch
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&,
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
                new solidContactFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~solidContactFvPatchVectorField();


    // Member functions

        // Access

            //- Return if this patch is the master
            Switch master() const
            {
                return master_;
            }

            //- Return the name of the shadow patch
            const word& shadowPatchName() const
            {
                return shadowPatchName_;
            }

            //- Return the name of the shadow patch
            label shadowPatchIndex() const
            {
                return shadowPatchIndex_;
            }

            //- Return index of the zone
            label zoneIndex() const;

            //- Return index of the shadow zone
            label shadowZoneIndex() const;

            //- Return reference to the normal contact model
            normalContactModel& normalModel();

            //- Return const reference to the normal contact model
            const normalContactModel& normalModel() const;

            //- Return reference to the friction contact model
            frictionContactModel& frictionModel();

            //- Return const reference to the friction contact model
            const frictionContactModel& frictionModel() const;

            //- Return const reference to deformed master zone
            const primitiveFacePatch& zone() const;

            //- Return reference to deformed master zone
            primitiveFacePatch& zone();

            //- Return const reference to deformed master zone
            const primitiveFacePatch& shadowZone() const;

            //- Return reference to deformed master zone
            primitiveFacePatch& shadowZone();

            //- Return const reference to patch-to-patch interpolation
            const extendedGgiZoneInterpolation& zoneToZone() const;

            //- Return reference to patch-to-patch interpolation
            extendedGgiZoneInterpolation& zoneToZone();

            //- Return reference to slave face normals
            //vectorField& slaveFaceNormals();

            //- Return const reference to global face zones
            const labelList& globalFaceZones() const;

            //- Return const reference to global face zone point map
            const IOList<labelList>& globalToLocalFaceZonePointMap() const;


        // Mapping functions

            //- Clear out demand driven data
            virtual void clearOut();

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

            //- Dissipated energy due to friction
            //virtual tmp<scalarField> Qc() const;

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "solidContactFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
