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

Class
    ggiInterfaceToInterfaceMapping

Description
    interfaceToInterfaceMapping wrapper for GGI interpolation

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    ggiInterfaceToInterfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef ggiInterfaceToInterfaceMapping_H
#define ggiInterfaceToInterfaceMapping_H

#include "interfaceToInterfaceMapping.H"
#ifdef FOAMEXTEND
    #include "ggiInterpolation.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

/*---------------------------------------------------------------------------*\
                        Class ggiInterfaceToInterfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class ggiInterfaceToInterfaceMapping
:
    public interfaceToInterfaceMapping
{
    // Private data

#ifdef FOAMEXTEND
        //- GGI interpolator between two standAlone patches
        mutable autoPtr< GGIInterpolation<standAlonePatch, standAlonePatch> >
            interpolatorPtr_;

    // Private Member Functions

        //- Make the GGI interpolator
        void makeInterpolator() const;

        //- Return a const access to the GGI interpolator
        const GGIInterpolation<standAlonePatch, standAlonePatch>&
            interpolator() const;

        //- Check zoneA to zoneB face interpolation error
        void checkZoneAToZoneBError() const;

        //- Check zoneB to zoneA point interpolation error
        void checkZoneBToZoneAError() const;

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
#endif

        //- Disallow default bitwise copy construct
        ggiInterfaceToInterfaceMapping(const ggiInterfaceToInterfaceMapping&);

        //- Disallow default bitwise assignment
        void operator=(const ggiInterfaceToInterfaceMapping&);


public:

    //- Runtime type information
    TypeName("GGI");


    // Constructors

        //- Construct from components
        ggiInterfaceToInterfaceMapping
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~ggiInterfaceToInterfaceMapping()
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

#ifdef FOAMEXTEND
#ifdef NoRepository
#   include "ggiInterfaceToInterfaceMappingTemplates.C"
#endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
