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

    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the normalContactModel and
    frictionContactModel.

    The buondary condition makes a copy of the master and slave face zones and
    keeps them in the deformed configuration; this allows the calculation of
    contact distances and interpolation between the zones.
    Note: we use const_cast to move the standAlonePatch zone points.

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
#include "extendedGgiInterpolation.H"
#include "Switch.H"
#include "boundBox.H"
#include "standAlonePatch.H"

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

        //- Is this patch the master
        const Switch master_;

        //- Write the deformed contact zones to VTK surface files
        const Switch writeZoneVTK_;

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

        // Master face zone kept in deformed configuration
        mutable standAlonePatch* zonePtr_;

        // Slave face zone kept in deformed configuration
        mutable standAlonePatch* shadowZonePtr_;

        //- Zone-to-zone interpolation
        mutable extendedGgiStandAlonePatchInterpolation* zoneToZonePtr_;

        //- GGI quick reject algorithm
        extendedGgiInterpolation::quickReject quickReject_;

        //- Region where the contact is expected to be. This aims to save
        // time in the GGI distance calculations
        //const boundBox regionOfInterest_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach being used
        bool movingMesh() const;

        //- Is a nonlinear geometry (i.e. finite strain) approach being used
        bool nonLinearGeometry() const;

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

        //- Create the zone
        void calcZone() const;

        //- Create the shadow zone
        void calcShadowZone() const;

        //- Create the zone-to-zone interpolator
        void calcZoneToZone() const;

        // Return the Ccntact dictionary
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
            const standAlonePatch& zone() const;

            //- Return reference to deformed master zone
            standAlonePatch& zone();

            //- Return const reference to deformed master zone
            const standAlonePatch& shadowZone() const;

            //- Return reference to deformed master zone
            standAlonePatch& shadowZone();

            //- Return const reference to patch-to-patch interpolation
            const extendedGgiStandAlonePatchInterpolation& zoneToZone() const;

            //- Return reference to patch-to-patch interpolation
            extendedGgiStandAlonePatchInterpolation& zoneToZone();

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
