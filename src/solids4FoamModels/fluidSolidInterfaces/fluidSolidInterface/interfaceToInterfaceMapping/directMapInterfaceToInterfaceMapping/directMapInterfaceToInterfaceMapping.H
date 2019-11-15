/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    directMapInterfaceToInterfaceMapping

Description
    interfaceToInterfaceMapping wrapper for direct mapping i.e. conformal
    meshes

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    directMapInterfaceToInterfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef directMapInterfaceToInterfaceMapping_H
#define directMapInterfaceToInterfaceMapping_H

#include "interfaceToInterfaceMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

/*---------------------------------------------------------------------------*\
                        Class directMapInterfaceToInterfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class directMapInterfaceToInterfaceMapping
:
    public interfaceToInterfaceMapping
{
    // Private data

        //- Face map from zoneB to zoneA
        //  For each zoneB face, give the index of the corresponding zoneA face
        mutable autoPtr<labelIOList> zoneAToZoneBFaceMapPtr_;

        //- Face map from zoneA to zoneB
        //  For each zoneA face, give the index of the corresponding zoneB face
        mutable autoPtr<labelIOList> zoneBToZoneAFaceMapPtr_;

        //- Point map from zoneB to zoneA
        //  For each zoneB point, give the index of the corresponding zoneA
        //  point
        mutable autoPtr<labelIOList> zoneAToZoneBPointMapPtr_;

        //- Point map from zoneA to zoneB
        //  For each zoneA point, give the index of the corresponding zoneB
        //  point
        mutable autoPtr<labelIOList> zoneBToZoneAPointMapPtr_;

        //- List of min edge lengths for interfaces points at fluid side
        mutable autoPtr<scalarField> minEdgeLengthsPtr_;

    // Private static data

        //  Relative tolerance (fraction of minimum edge length used when
        //  checking for matched faces/points
        static const scalar relTol_;

    // Private Member Functions

        //- Check zone face and points sizes are conformal
        void checkZoneSizes() const;

        //- Calculate zoneA to zoneB face map
        void calcZoneAToZoneBFaceMap() const;

        //- Return a reference to the zoneA to zoneB face map
        const labelIOList& zoneAToZoneBFaceMap() const;

        //- Calculate zoneB to zoneA face map
        void calcZoneBToZoneAFaceMap() const;

        //- Return a reference to the zoneB to zoneA face map
        const labelIOList& zoneBToZoneAFaceMap() const;

        //- Calculate zoneA to zoneB point map
        void calcZoneAToZoneBPointMap() const;

        //- Return a reference to the zoneA to zoneB point map
        const labelIOList& zoneAToZoneBPointMap() const;

        //- Calculate zoneB to zoneA point map
        void calcZoneBToZoneAPointMap() const;

        //- Return a reference to the zoneB to zoneA point map
        const labelIOList& zoneBToZoneAPointMap() const;

        //- Calculate minimal edge lengths
        void calcMinEdgeLengths() const;

        //- Return minimal edge lengths
        const scalarField& minEdgeLengths() const;

        //- Transfer/map/interpolate from one zone faces to another zone
        //  faces for Type
        template<class Type>
        void transferFacesZoneToZone
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField               // to field
        ) const;

        //- Transfer/map/interpolate from one zone points to another zone
        //  points for Type
        template<class Type>
        void transferPointsZoneToZone
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField               // to field
        ) const;

        //- Disallow default bitwise copy construct
        directMapInterfaceToInterfaceMapping
        (
            const directMapInterfaceToInterfaceMapping&
        );

        //- Disallow default bitwise assignment
        void operator=(const directMapInterfaceToInterfaceMapping&);


public:

    //- Runtime type information
    TypeName("directMap");


    // Constructors

        //- Construct from components
        directMapInterfaceToInterfaceMapping
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~directMapInterfaceToInterfaceMapping()
        {}


    // Member Functions

        // Edit

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for scalars
            virtual void transferFacesZoneToZone
            (
                const standAlonePatch& fromZone,    // from zone
                const standAlonePatch& toZone,      // to zone
                const Field<scalar>& fromField,     // from field
                Field<scalar>& toField              // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for scalars
            virtual void transferPointsZoneToZone
            (
                const standAlonePatch& fromZone,    // from zone
                const standAlonePatch& toZone,      // to zone
                const Field<scalar>& fromField,     // from field
                Field<scalar>& toField              // to field
            ) const;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for vectors
            virtual void transferFacesZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for vectors
            virtual void transferPointsZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceToInterfaceMappings


} // End namespace Foam

#ifdef NoRepository
#   include "directMapInterfaceToInterfaceMappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
