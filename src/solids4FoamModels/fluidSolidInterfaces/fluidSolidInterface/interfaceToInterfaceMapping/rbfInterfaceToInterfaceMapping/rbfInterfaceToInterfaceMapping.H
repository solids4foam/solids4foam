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
    rbfInterfaceToInterfaceMapping

Description
    interfaceToInterfaceMapping wrapper using radial basis functions

Author
    Philip Cardiff, UCD. All rights reserved.
    This class is a wrapper for the code from David Blom

SourceFiles
    rbfInterfaceToInterfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef rbfInterfaceToInterfaceMapping_H
#define rbfInterfaceToInterfaceMapping_H

#include "interfaceToInterfaceMapping.H"
#include "RBFInterpolation.H"

using namespace rbf;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace interfaceToInterfaceMappings
{

/*---------------------------------------------------------------------------*\
                        Class rbfInterfaceToInterfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class rbfInterfaceToInterfaceMapping
:
    public interfaceToInterfaceMapping
{
    // Private data

        //- zoneA to zoneB interpolator
        mutable std::shared_ptr<RBFInterpolation> zoneAToZoneBInterpolatorPtr_;

        //- List of solid zone to fluid zone interpolators
        mutable std::shared_ptr<RBFInterpolation> zoneBToZoneAInterpolatorPtr_;

    // Private Member Functions

        //- Make zoneA to zoneB interpolator
        void makeZoneAToZoneBInterpolator() const;

        //- Return reference to zoneA to zoneB interpolator
        const std::shared_ptr<RBFInterpolation>&
            zoneAToZoneBInterpolator() const;

        //- Make zoneB to zoneA interpolator
        void makeZoneBToZoneAInterpolator() const;

        //- Return reference to zoneB to zoneA interpolator
        const std::shared_ptr<RBFInterpolation>&
            zoneBToZoneAInterpolator() const;

        //- Transfer/map/interpolate from one zone faces to another zone
        //  faces for Type
        template<class Type>
        void transferFacesZoneToZone
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField                // to field
        ) const;

        //- Transfer/map/interpolate from one zone points to another zone
        //  points for Type
        template<class Type>
        void transferPointsZoneToZone
        (
            const standAlonePatch& fromZone,    // from zone
            const standAlonePatch& toZone,      // to zone
            const Field<Type>& fromField,       // from field
            Field<Type>& toField                // to field
        ) const;

        //- Disallow default bitwise copy construct
        rbfInterfaceToInterfaceMapping
        (
            const rbfInterfaceToInterfaceMapping&
        );

        //- Disallow default bitwise assignment
        void operator=(const rbfInterfaceToInterfaceMapping&);


public:

    //- Runtime type information
    TypeName("RBF");


    // Constructors

        //- Construct from components
        rbfInterfaceToInterfaceMapping
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~rbfInterfaceToInterfaceMapping()
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
#   include "rbfInterfaceToInterfaceMappingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
