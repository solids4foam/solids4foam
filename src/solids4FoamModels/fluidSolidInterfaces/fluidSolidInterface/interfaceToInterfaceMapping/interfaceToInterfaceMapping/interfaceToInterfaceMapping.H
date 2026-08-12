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
    interfaceToInterfaceMapping

Description
    Virtual base class for mapping/interpolating between two interfaces

    The interface-to-interface correspondence is topological in body-fitted
    FSI: the same solid face lies beside the same fluid face for the duration
    of the run.  It is therefore constructed once from the undeformed geometry
    rather than from the current interface, which follows the fluid mesh motion
    and the solid deformation.  Building it from the current interface would
    make the correspondence depend on when it happens to be first used.

    The undeformed geometry is read from the mesh points in the "constant"
    instance, so every mesh whose interface is mapped must have its points
    available there and must not change topology during the run.  Both
    conditions are checked in makeReferenceZone(), which aborts if they are not
    met: silently falling back to the current geometry would reintroduce the
    ordering dependence described above, and a point-count mismatch may signal
    a topology change for which no correspondence is well defined.

Author
    Philip Cardiff, UCD.

SourceFiles
    interfaceToInterfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceToInterfaceMapping_H
#define interfaceToInterfaceMapping_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "primitivePatch.H"
#include "globalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceToInterfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class interfaceToInterfaceMapping
{
    // Private data

        //- Copy of settings dictionary
        const dictionary dict_;

        //- Const reference to the patchA
        const primitivePatch& patchA_;

        //- Const reference to the patchB
        const primitivePatch& patchB_;

        //- Const reference to the globalPatchA
        const globalPolyPatch& globalPatchA_;

        //- Const reference to the globalPatchB
        const globalPolyPatch& globalPatchB_;

        // Demand-driven data

            //- Reference-configuration copies used to construct the mapping
            mutable autoPtr<standAlonePatch> zoneARefPtr_;
            mutable autoPtr<standAlonePatch> zoneBRefPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceToInterfaceMapping(const interfaceToInterfaceMapping&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceToInterfaceMapping&);

        //- Build a zone from the constant mesh points
        autoPtr<standAlonePatch> makeReferenceZone
        (
            const globalPolyPatch& globalPatch
        ) const;


protected:

    // Protected Member Functions

        //- Return zoneA in the reference configuration
        const standAlonePatch& zoneARef() const;

        //- Return zoneB in the reference configuration
        const standAlonePatch& zoneBRef() const;


public:

    //- Runtime type information
    TypeName("interfaceToInterfaceMapping");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfaceToInterfaceMapping,
            dictionary,
            (
                const word& type,
                const dictionary& dict,
                const primitivePatch& patchA,
                const primitivePatch& patchB,
                const globalPolyPatch& globalPatchA,
                const globalPolyPatch& globalPatchB
            ),
            (type, dict, patchA, patchB, globalPatchA, globalPatchB)
        );


    // Constructors

        //- Construct from components
        interfaceToInterfaceMapping
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Selectors

        //- Select constructed from fluid and solid meshes
        static autoPtr<interfaceToInterfaceMapping> New
        (
            const word& type,
            const dictionary& dict,
            const primitivePatch& patchA,
            const primitivePatch& patchB,
            const globalPolyPatch& globalPatchA,
            const globalPolyPatch& globalPatchB
        );


    // Destructor

        virtual ~interfaceToInterfaceMapping();


public:


    // Member Functions

        // Access

            //- Return const access to the settings dict
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Const access to the patchA
            const primitivePatch& patchA() const
            {
                return patchA_;
            }

            //- Const access to the patchB
            const primitivePatch& patchB() const
            {
                return patchB_;
            };

            //- Const access to the globalPatchA
            const globalPolyPatch& globalPatchA() const
            {
                return globalPatchA_;
            }

            //- Const access to the globalPatchB
            const globalPolyPatch& globalPatchB() const
            {
                return globalPatchB_;
            };

            //- Const access to the zoneA
            const standAlonePatch& zoneA() const
            {
                return globalPatchA_.globalPatch();
            }

            //- Const access to the zoneB
            const standAlonePatch& zoneB() const
            {
                return globalPatchB_.globalPatch();
            };

            //- Check field sizes are correct before interpolation/mapping
            void checkFieldSizes
            (
                const label& fromZoneSize,
                const label& toZoneSize,
                const label fromFieldSize,
                const label toFieldSize
            ) const;

            //  Note: unfortunatley, a template virtual function is not allowed
            //  in C++, so instead we will make an explicit implementation for
            //  scalars and vectors; derived classes can use templates to make
            //  their lives easier if desired

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for scalars
            virtual void transferFacesZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<scalar>& fromField,       // from field
                Field<scalar>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for scalars
            virtual void transferPointsZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<scalar>& fromField,       // from field
                Field<scalar>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone faces to another zone
            //  faces for vectors
            virtual void transferFacesZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const = 0;

            //- Transfer/map/interpolate from one zone points to another zone
            //  points for vectors
            virtual void transferPointsZoneToZone
            (
                const standAlonePatch& fromZone,      // from zone
                const standAlonePatch& toZone,        // to zone
                const Field<vector>& fromField,       // from field
                Field<vector>& toField                // to field
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
