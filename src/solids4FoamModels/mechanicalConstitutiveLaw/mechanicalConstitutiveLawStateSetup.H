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
    Foam::mechanicalConstitutiveLawStateSetup

Description
    Sets up the state of each mechanical constitutive law of a
    mechanicalConstitutiveLawManager, on a law's integration points, on its
    faces of a boundary patch, or in a scratch state that no field addresses.

    What a law declared in declareState is allocated and set to its declared
    default; the law then initialises what it wants to itself; and a
    prescribed field is read over that from a field the case supplied. A
    composite's sub-laws are prepared the same way, in child states. For a
    restart, the persistent state of every law is registered for writing and
    read back, together with where each point sits on the mesh.

    Split out of the manager, which holds one and forwards to it. It holds
    references to the manager's mesh, laws, law names, law cells, law
    boundary faces, the registered objects that write state, and whether the
    kinematic history needed for a restart is available, and holds no state
    of its own.

SourceFiles
    mechanicalConstitutiveLawStateSetup.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawStateSetup_H
#define mechanicalConstitutiveLawStateSetup_H

#include "fvMesh.H"
#include "volFields.H"
#include "regIOobject.H"
#include "HashSet.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawState.H"
#include "mechanicalConstitutiveLawStateSpec.H"
#include "integrationPointTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class mechanicalConstitutiveLawStateSetup Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawStateSetup
{
    // Private Data

        //- The mesh the laws are defined on
        const fvMesh& mesh_;

        //- The laws, as the manager holds them
        const PtrList<mechanicalConstitutiveLaw>& laws_;

        //- Name of each law, which qualifies the state written for restart
        const wordList& lawNames_;

        //- The cells of each law
        const List<labelList>& lawCells_;

        //- Each law's faces on each patch, indexing into the patch
        const List<List<labelList>>& lawBoundaryFaces_;

        //- Registered objects that write the persistent state; the manager
        //  owns them, since a registered object must outlive its writes
        PtrList<regIOobject>& stateWriteProxies_;

        //- Whether a continued run has the kinematic history a constitutive
        //  history is measured against, as the solid model says
        const bool& restartKinematicsAvailable_;


    // Private Member Functions

        //- Where a prescribed field is read from.
        //
        //  The internal points of a law and its faces on one patch are read
        //  by different functions but through the same per-type dispatch, so
        //  the dispatch is written once and told which to use.
        //
        //  Note that writing state out for restart is not a third source.
        //  This reads into a field that already exists, at
        //  addresses it is given; writing state out has to register an
        //  IOobject and scatter per-integration-point state into a geometric
        //  field, which is the other direction and a different shape
        enum class prescribedSource
        {
            //- The law's own integration points, from a volField
            internalPoints,

            //- The law's faces on one boundary patch, from its patch values
            patchFaces
        };

        //- Read every prescribed field a law declared, from one source
        void readPrescribedFields
        (
            const mechanicalConstitutiveLawStateSpec& spec,
            mechanicalConstitutiveLawState& state,
            const prescribedSource source,
            const label lawI,
            const integrationPointTopology* topoPtr,
            const labelList* ipIDsPtr,
            const label patchI
        ) const;

        //- Refuse a restart that would restore constitutive history without
        //  the kinematic history it is measured against
        void checkRestartKinematics() const;

        //- One law's persistent state, recursing into any child states
        void setupStateRestartLaw
        (
            const mechanicalConstitutiveLaw& law,
            const word& lawName,
            const word& topologyName,
            const wordList& childPath,
            const word& entityName,
            const labelUList& entities,
            const List<mechanicalConstitutiveLawState*>& parts,
            const bool isRestart
        ) const;

        //- Allocate each declared field and set it to its declared default
        void applyStateDefaults
        (
            const mechanicalConstitutiveLawStateSpec& spec,
            mechanicalConstitutiveLawState& state
        ) const;

        //- As applyStateSpec, for a named law, which is how a composite's
        //  sub-laws are reached
        void applyStateSpec
        (
            const mechanicalConstitutiveLaw& law,
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            mechanicalConstitutiveLawState& state
        ) const;

        //- As applyStateSpecPatch, for a named law
        void applyStateSpecPatch
        (
            const mechanicalConstitutiveLaw& law,
            const label lawI,
            const label patchI,
            mechanicalConstitutiveLawState& state
        ) const;

        //- As applyStateSpecScratch, for a named law
        void applyStateSpecScratch
        (
            const mechanicalConstitutiveLaw& law,
            mechanicalConstitutiveLawState& state
        ) const;

        //- The volField a prescribed state is to be read from, if there is
        //  one. It may be a field another model has already registered, or it
        //  may exist only as a file in the current time directory, which is
        //  how a user supplies one. If it is neither then nothing was
        //  prescribed and the caller leaves the declared default in place
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        prescribedField(const word& name) const;

        //- Fill a law's internal integration points from a prescribed field.
        //  The value is given per cell; gatherToIntegrationPoints says how it
        //  reaches each point
        template<class Type>
        void readPrescribed
        (
            const word& name,
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            Field<Type>& fld
        ) const;

        //- Read a prescribed field into a law's state on one boundary patch.
        //  A boundary state is evaluated exactly as the internal points are,
        //  so it needs the prescribed value too; without this the law would
        //  see a silently zero field on every patch
        template<class Type>
        void readPrescribedPatch
        (
            const word& name,
            const label lawI,
            const label patchI,
            Field<Type>& fld
        ) const;


public:

    // Constructors

        //- Construct from references to the manager's members, which must
        //  outlive this
        mechanicalConstitutiveLawStateSetup
        (
            const fvMesh& mesh,
            const PtrList<mechanicalConstitutiveLaw>& laws,
            const wordList& lawNames,
            const List<labelList>& lawCells,
            const List<List<labelList>>& lawBoundaryFaces,
            PtrList<regIOobject>& stateWriteProxies,
            const bool& restartKinematicsAvailable
        );


    // Member Functions

        //- True if this law, or any law below it, keeps history
        bool declaresPersistentState
        (
            const mechanicalConstitutiveLaw& law
        ) const;

        //- Register the persistent state of every law on this topology for
        //  writing, and on a restart read it back
        //  Named by the key the topology was registered under, not by its
        //  class name: two topologies of one class under different keys keep
        //  their own history, and naming the files by class would have them
        //  overwrite and then restore each other's
        void setupStateRestart
        (
            const List<labelList>& lawIntegrationPointIDs,
            PtrList<mechanicalConstitutiveLawState>& states,
            PtrList<PtrList<mechanicalConstitutiveLawState>>& boundaryStates,
            const bool boundaryAware,
            const integrationPointTopology& topo,
            const word& topologyKey
        ) const;

        //- Apply what a law declared in declareState: allocate each field,
        //  set its default, and read any prescribed field from the mesh
        void applyStateSpec
        (
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            mechanicalConstitutiveLawState& state
        ) const;

        //- As applyStateSpec, for a scratch state that no field addresses:
        //  declared defaults and the law's own initialisation, no
        //  prescribed read
        void applyStateSpecScratch
        (
            const label lawI,
            mechanicalConstitutiveLawState& state
        ) const;

        //- As applyStateSpec, for a law's state on one boundary patch
        void applyStateSpecPatch
        (
            const label lawI,
            const label patchI,
            mechanicalConstitutiveLawState& state
        ) const;

        //- Write the history held on a cell-centred topology as volFields,
        //  for viewing: one field per persistent state variable, covering
        //  every material that declares it. A field another model has
        //  registered is not overwritten, and is reported once through
        //  reportedStateFields
        void writeStateFields
        (
            const List<labelList>& lawIntegrationPointIDs,
            const PtrList<mechanicalConstitutiveLawState>& states,
            const PtrList<PtrList<mechanicalConstitutiveLawState>>&
                boundaryStates,
            const bool boundaryAware,
            HashSet<word>& reportedStateFields
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
