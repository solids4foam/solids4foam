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
    Foam::mechanicalConstitutiveLawManager

Description
    Manager class responsible for applying mechanical constitutive laws
    to the computational fields.

    The mechanicalConstitutiveLawManager owns one or more mechanical
    constitutive laws together with their associated constitutive state
    objects. It is responsible for extracting kinematic fields from the
    solver, constructing the appropriate kinematics containers, and
    invoking constitutive evaluation.

    This class forms the interface between OpenFOAM field types
    (e.g. volTensorField, pointTensorField) and the discretisation-
    agnostic mechanicalConstitutiveLaw interface. All concerns related
    to field location, mesh topology, and region selection are handled
    exclusively at this level.

    In its initial form, the manager supports a single constitutive law
    applied to the entire domain. This functionality may be extended
    incrementally to support multiple laws, multiple regions, boundary
    contributions, and tangent-related quantities.

    Notes:
    - Constitutive laws are stateless; all history variables are stored
      in mechanicalConstitutiveLawState objects owned by the manager.
    - The manager is the only component that depends on OpenFOAM
      GeometricField types.

SourceFiles
    mechanicalConstitutiveLawManager.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalConstitutiveLawManager_H
#define mechanicalConstitutiveLawManager_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawState.H"
#include "mechanicalConstitutiveLawStateSpec.H"
#include "mechanicalConstitutiveLawStressCollapseRule.H"

#ifdef OPENFOAM_ORG
#include "CompactListList.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class integrationPointTopology;
class CompactListListBase;
struct mat66;

//- OpenFOAM.org's CompactListList has two template parameters, so its header
//  is included above. The other forks have one template parameter; forward
//  declare it to avoid instantiating foam-extend's broken const operator[].
#ifndef OPENFOAM_ORG
template<class T> class CompactListList;
#endif

/*---------------------------------------------------------------------------* \
              Class mechanicalConstitutiveLawManager Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawManager
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Pointer to the ooint mesh, if needed
        mutable const pointMesh* pMeshPtr_;

        //- Current time index so we can track when there is a new time step
        label curTimeIndex_;

        //- Constitutive laws
        PtrList<mechanicalConstitutiveLaw> laws_;

        //- Name of each law, as the mechanicalProperties dictionary names it.
        //  This is the material name, and it qualifies the state written for
        //  restart: two materials may each declare epsilonP and mean something
        //  different by it
        wordList lawNames_;

        //- Registered objects that write the persistent state.
        //  Held here because a registered object must outlive the writes it
        //  serves, and the registry does not own these
        mutable PtrList<regIOobject> stateWriteProxies_;

        //- Cell indices defining the material support for each mechanical
        //  constitutive law.
        //  Materials in solids4foam are always defined per cell (via cellZones),
        //  independent of the discretisation order or the locations at which
        //  stresses or strains are evaluated.
        //  These lists are used to construct integration-point mappings and to
        //  assemble cell-based material properties (e.g. density, bulk modulus).
        List<labelList> lawCells_;

        //- Boundary patch face indices per mechanical constitutive law.
        //  For each mechanicalConstitutiveLaw, this stores a list over boundary
        //  patches, where each entry contains the local face indices (within that
        //  patch) whose owner cell belongs to the corresponding law.
        //  In other words:
        //    lawBoundaryFaces_[lawI][patchI] gives the list of boundary faces on
        //    patch patchI whose adjacent cell is assigned to law lawI.
        List<List<labelList>> lawBoundaryFaces_;

        //- Scratch fields for stress accumulation on surface fields
        mutable autoPtr<surfaceSymmTensorField> surfaceStressSumPtr_;
        mutable autoPtr<surfaceScalarField> surfaceStressWeightPtr_;
        mutable autoPtr<surfaceScalarField> surfaceTangentWeightPtr_;

        //- Scratch fields for stress accumulation on point fields
        mutable autoPtr<pointSymmTensorField> pointStressSumPtr_;
        mutable autoPtr<pointScalarField>     pointStressWeightPtr_;
        mutable autoPtr<pointScalarField>     pointTangentWeightPtr_;

        //- Scratch stress storage used by the tangent-only updates.
        //  A constitutive law always produces a stress alongside a tangent, so
        //  a tangent-only update needs somewhere to put it that is not the
        //  caller's storage. The contents are never read. The list is sized
        //  on demand, which in practice means once per topology size
        mutable List<symmTensor> scratchStress_;

        //- Density field
        //- Mechanical constitutive law index of each cell.
        //  Needed by the flat-list boundary update for patches that carry no
        //  fvPatch faces, where lawBoundaryFaces_ cannot supply the law
        labelList cellToLaw_;

        mutable autoPtr<volScalarField> rhoPtr_;

        //- Bulk modulus field
        mutable autoPtr<volScalarField> kappaPtr_;

        //- Cache of integration-point topologies, owned by the manager
        mutable HashTable<autoPtr<integrationPointTopology>> topologyCache_;

        //- Topology-specific constitutive data owned by the manager.
        //  A topologyEntry groups together all data that depends on the *location*
        //  of integration points (e.g. cell-centred, face-centred, vertex-centred,
        //  or higher-order multi-point layouts).
        //  For a given integrationPointTopology, this includes:
        //    - the mapping from mechanical laws to integration-point indices
        //      (lawIntegrationPointIDs),
        //    - the constitutive state objects storing history variables at those
        //      integration points, and
        //    - optional boundary constitutive states where applicable.
        //  This indirection allows a single mechanicalConstitutiveLawManager to
        //  evaluate stresses at multiple locations (cells, faces, points, etc.)
        //  using the same material definitions, while keeping state variables
        //  strictly separated per topology.
        //  Topology entries are constructed lazily on first use and cached for the
        //  lifetime of the manager.
        struct topologyEntry
        {
            //- The topology that defines this entry (non-owning)
            const integrationPointTopology& topology_;

            //- Integration-point indices per mechanical law
            List<labelList> lawIntegrationPointIDs_;

            //- Constitutive state per mechanical law
            PtrList<mechanicalConstitutiveLawState> states_;

            //- Boundary constitutive state per law per patch (optional)
            PtrList<PtrList<mechanicalConstitutiveLawState>> boundaryStates_;

            //- True if this topology supports boundary integration points
            bool boundaryAware_;

            //- Scratch storage for live coupling inputs, per law, keyed by
            //  field name. A law that asks for one is handed a view of its own
            //  integration points, and the storage is kept here rather than
            //  allocated per evaluation
            List<HashTable<autoPtr<scalarField>>> lawScalarInputs_;

            //- Constructor
            topologyEntry(const integrationPointTopology& topo)
            :
                topology_(topo),
                lawIntegrationPointIDs_(),
                states_(),
                boundaryStates_(),
                boundaryAware_(false),
                lawScalarInputs_()
            {}
        };

        //- Cache of topology-specific integration-point data and state
        //  Keyed by the address of the integrationPointTopology instance
        mutable HashTable<autoPtr<topologyEntry>> topologyEntries_;


    // Private Member Functions

        //- Check the provided mesh is the same as the mesh used during
        //  construction
        inline void checkMeshConsistency
        (
            const polyMesh& expectedMesh,
            const polyMesh& fieldMesh,
            const word& fieldName
        )
        {
            if (&fieldMesh != &expectedMesh)
            {
                FatalError
                    << "Field '" << fieldName << "' is associated with a "
                    << "different mesh." << nl
                    << "Expected mesh address: " << expectedMesh.name() << nl
                    << "Field mesh address:    " << fieldMesh.name()
                    << exit(FatalError);
            }
        }

        //- Check the provided compact list sizes are consistent.
        //  Defined in the .C so that this header does not need the complete
        //  CompactListList type
        void checkCompactLayoutConsistency
        (
            const CompactListList<tensor>& a,
            const CompactListList<tensor>& b,
            const CompactListList<symmTensor>& out,
            const List<scalar>* tangentPtr,
            const word& context
        );

        //- Check that a flat integration-point list has the size the topology
        //  requires
        inline void checkIntegrationPointListSize
        (
            const label expected,
            const label actual,
            const word& listName,
            const word& context
        ) const
        {
            if (actual != expected)
            {
                FatalErrorInFunction
                    << "Inconsistent list size in " << context << nl
                    << "List '" << listName << "' has " << actual
                    << " entries but the integration-point topology has "
                    << expected << " integration points"
                    << exit(FatalError);
            }
        }

        //- Check that the tangent storage required by the request has been
        //  supplied
        void checkTangentStorage
        (
            const bool haveScalarTangent,
            const bool haveFourthOrderTangent,
            const tangentRequest req,
            const word& context
        );

        //- Return scratch stress storage sized for the given number of
        //  integration points
        List<symmTensor>& scratchStress(const label nIntegrationPoints) const;

        //- Return (and lazily construct) a compact cell-based integration-point
        //  topology inferred from a CompactListList layout.
        //  The topology is owned and cached by the manager and remains valid
        //  for the lifetime of the manager.
        //  The layout is used only for addressing; the value type is irrelevant
        const integrationPointTopology& compactCellTopologyFor
        (
            const CompactListList<tensor>& layout
        ) const;

        //- Return (and lazily construct) a compact face-based integration-point
        //  topology inferred from a CompactListList layout.
        //  Intended for discretisations with face-based integration points.
        //  The layout is used only for addressing; the value type is irrelevant
        const integrationPointTopology& compactFaceTopologyFor
        (
            const CompactListList<tensor>& layout
        ) const;

        //- Ensure a topology entry exists and return it
        topologyEntry& topology
        (
            const integrationPointTopology& topo
        ) const;

        //- Commit old-time state for all topology entries if needed
        //- Where a prescribed field is read from.
        //
        //  The internal points of a law and its faces on one patch are read
        //  by different functions but through the same per-type dispatch, so
        //  the dispatch is written once and told which to use.
        //
        //  Note that the restart writing of section 5 of DESIGN-state-io.md is
        //  not a third source. This reads into a field that already exists, at
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

        //- Work out what to evaluate for one law on one boundary patch:
        //  which integration points, and which state.
        //
        //  Returns false when there is nothing to evaluate there. The state
        //  may be the law's real boundary state, a shadow of it when the
        //  caller is preserving state, or a scratch state for an empty patch
        //  addressed through its polyPatch; the two autoPtrs own the latter
        //  two and must outlive the evaluation
        bool resolveBoundaryEvaluation
        (
            const integrationPointTopology& topo,
            topologyEntry& tp,
            const label lawI,
            const label patchI,
            const bool preserveState,
            labelList& ipIDs,
            autoPtr<mechanicalConstitutiveLawState>& bScratchPtr,
            autoPtr<mechanicalConstitutiveLawState>& bShadowPtr,
            mechanicalConstitutiveLawState*& bState
        ) const;

        //- Gather the live coupling inputs a law asked for, as views over its
        //  own integration points, and hand back an inputs object carrying
        //  them. A law that asks for nothing gets an inputs object holding
        //  only the time step, which is what every law had before
        mechanicalConstitutiveLawInputs lawInputs
        (
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            const scalar dt,
            topologyEntry& tp
        ) const;

        //- As lawInputs, for a law's faces on one boundary patch. The values
        //  come from the field's patch values, not from the cells behind it
        mechanicalConstitutiveLawInputs lawInputsPatch
        (
            const label lawI,
            const label patchI,
            const labelList& faces,
            const scalar dt,
            topologyEntry& tp
        ) const;

        //- An inputs object for an evaluation path that does not gather
        //  coupling inputs yet. Fatal if any law asks for one, because the
        //  alternative is handing that law an empty inputs object and letting
        //  it read a zero it cannot distinguish from a real value
        mechanicalConstitutiveLawInputs inputsWithoutCoupling
        (
            const scalar dt,
            const word& path
        ) const;

        //- True if this law, or any law below it, keeps history
        bool declaresPersistentState
        (
            const mechanicalConstitutiveLaw& law
        ) const;

        //- Refuse a restart that would restore constitutive history without
        //  the kinematic history it is measured against
        void checkRestartKinematics() const;

        //- Register the persistent state of every law on this topology for
        //  writing, and on a restart read it back
        void setupStateRestart
        (
            topologyEntry& entry,
            const integrationPointTopology& topo
        ) const;

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

        //- Apply what a law declared in declareState: allocate each field,
        //  set its default, and read any prescribed field from the mesh
        void applyStateSpec
        (
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
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

        //- Whether a file for this IOobject exists and holds this type.
        //  foam-extend has no typed header check, so the class name is
        //  compared after reading the header, which is what the typed check
        //  does anyway
        template<class FieldType>
        static bool headerIsA(IOobject& io)
        {
#ifdef FOAMEXTEND
            return io.headerOk() && io.headerClassName() == FieldType::typeName;
#else
            return io.typeHeaderOk<FieldType>(true);
#endif
        }

        //- The volField a prescribed state is to be read from, if there is
        //  one. It may be a field another model has already registered, or it
        //  may exist only as a file in the current time directory, which is
        //  how a user supplies one. If it is neither then nothing was
        //  prescribed and the caller leaves the declared default in place
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        prescribedField(const word& name) const
        {
            typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

            // The file comes first, and the registry only after it.
            //
            // The other way round looks more natural and is wrong: the legacy
            // mechanical law registers its own sigma0, read from the current
            // time directory alone, so on a restart there is a registered
            // field of zeros sitting in front of the file the user actually
            // wrote. Reading the file first means what the case says on disk
            // is what the law gets, and the registry is left to serve the
            // fields that only ever exist in memory because another model
            // computes them
            IOobject io
            (
                name,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (!headerIsA<VolFieldType>(io))
            {
                // A prescribed field describes the case, not the state at a
                // particular time, so it is written once into 0 and a restart
                // will not find it beside the fields it restarts from. Look
                // there before giving up, which is what the legacy Guccione
                // law does for its fibre directions
                io.instance() = "0";
            }

            if (headerIsA<VolFieldType>(io))
            {
                return tmp<VolFieldType>(new VolFieldType(io, mesh_));
            }

            if (mesh_.foundObject<VolFieldType>(name))
            {
                return tmp<VolFieldType>
                (
                    mesh_.lookupObject<VolFieldType>(name)
                );
            }

            return tmp<VolFieldType>();
        }

        //- Fill a law's internal integration points from a prescribed field
        //
        //  The value is specified per cell, which is how a user prepares one
        //  and the only form setFields can write. How it reaches a point
        //  depends on where the point is. A point inside a cell takes that
        //  cell's value. A point that is a mesh face takes the interpolated
        //  value, because that is what legacy does and because a face lies
        //  between two cells rather than in one.
        //
        //  A topology that shares a point between cells without being able to
        //  say it is a face still broadcasts, and there the last cell to
        //  reach the point wins: exact where the field is continuous across
        //  the interface, ambiguous where it is not. That is the mapping
        //  policy question recorded in DESIGN-state-io.md section 4
        template<class Type>
        void readPrescribed
        (
            const word& name,
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            Field<Type>& fld
        ) const
        {
            typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

            const tmp<VolFieldType> tsrc(prescribedField<Type>(name));

            if (!tsrc.valid())
            {
                return;
            }

            const VolFieldType& src = tsrc();

            Info<< "    Prescribed state '" << name
                << "' read from the field of the same name" << endl;

            gatherToIntegrationPoints(src, lawI, topo, ipIDs, fld);
        }

        //- Fill a law's integration points from a volField on the mesh.
        //
        //  The value is specified per cell, which is how a user prepares one
        //  and how another model produces one. How it reaches a point depends
        //  on where the point is. A point inside a cell takes that cell's
        //  value. A point that is a mesh face takes the interpolated value,
        //  because a face lies between two cells rather than in one, and
        //  because that is what the legacy code does.
        //
        //  A topology that shares a point between cells without being able to
        //  say it is a face still broadcasts, and there the last cell to reach
        //  the point wins: exact where the field is continuous across the
        //  interface, ambiguous where it is not. That is the mapping policy
        //  question recorded in DESIGN-state-io.md section 4
        template<class Type>
        void gatherToIntegrationPoints
        (
            const GeometricField<Type, fvPatchField, volMesh>& src,
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            Field<Type>& fld
        ) const
        {
            if (topo.integrationPointsAreFaces())
            {
                // The integration point index is the face index here, which is
                // what lets the interpolated field be read straight off
                const GeometricField<Type, fvsPatchField, surfaceMesh> srcf
                (
                    linearInterpolate(src)
                );

                const label nInternal = mesh_.nInternalFaces();

                forAll(ipIDs, localI)
                {
                    const label faceI = ipIDs[localI];

                    if (faceI < nInternal)
                    {
                        fld[localI] = srcf[faceI];
                    }
                    else
                    {
                        const label patchI =
                            mesh_.boundaryMesh().whichPatch(faceI);

                        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

                        fld[localI] =
                            srcf.boundaryField()[patchI][faceI - pp.start()];
                    }
                }

                return;
            }

            // A state field is indexed by the law's own integration points,
            // 0 to n-1, while the topology hands out global indices, so the
            // two have to be related before anything can be written
            labelList globalToLocal(topo.nIntegrationPoints(), -1);
            forAll(ipIDs, localI)
            {
                globalToLocal[ipIDs[localI]] = localI;
            }

            const labelList& cells = lawCells_[lawI];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                const labelUList ips = topo.cellIntegrationPointIDs(cellI);

                forAll(ips, j)
                {
                    const label localI = globalToLocal[ips[j]];

                    // A point of this cell that is not one of this law's own
                    // points cannot happen while materials are assigned per
                    // cell, but the guard costs nothing and the alternative is
                    // writing outside the field
                    if (localI >= 0)
                    {
                        fld[localI] = src[cellI];
                    }
                }
            }
        }

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
        ) const
        {
            typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

            const tmp<VolFieldType> tsrc(prescribedField<Type>(name));

            if (!tsrc.valid())
            {
                return;
            }

            const VolFieldType& src = tsrc();

            // The patch values of the supplied field, which for the uniform
            // and zero-gradient cases a user writes are the owner cell values
            const fvPatchField<Type>& psrc = src.boundaryField()[patchI];

            // This law's faces on this patch, indexing into the patch
            const labelList& faces = lawBoundaryFaces_[lawI][patchI];

            forAll(faces, i)
            {
                fld[i] = psrc[faces[i]];
            }
        }

        void updateOldTimeIfNeeded();


        //- No copy construct
        mechanicalConstitutiveLawManager
        (
            const mechanicalConstitutiveLawManager&
        ) = delete;

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawManager&) = delete;


public:

    // Public static predicates on a tangent request

        //- Check if the provided tangent request is scalar or scalarDeviatoric
        inline static bool needsScalarTangent(const tangentRequest req)
        {
            return
                req == tangentRequest::scalar
             || req == tangentRequest::scalarDeviatoric;
        }

        //- Check if the provided tangent request is fourthOrder or
        //  fourthOrderFiniteDifference
        inline static bool needsFourthOrderTangent(const tangentRequest req)
        {
            return
                req == tangentRequest::fourthOrder
             || req == tangentRequest::fourthOrderFiniteDifference;
        }

        //- Flat-list small-strain evaluation.
        //  Shared implementation of the public stress and tangent-only
        //  updates. With preserveState, each law is evaluated against a
        //  shadow of its state, so history is read but never written
        void evaluateSmallStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& gradD,
            const UList<tensor>& gradD0,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr,
            UList<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq,
            const bool preserveState,
            const bool coldState = false
        );

        //- Finite-strain counterpart of evaluateSmallStrain
        void evaluateFiniteStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& F,
            const UList<tensor>& F0,
            const UList<tensor>& Finv,
            const UList<tensor>& Finv0,
            const UList<scalar>& J,
            const UList<scalar>& J0,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr,
            UList<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq,
            const bool preserveState,
            UList<scalar>* volumetricPtr = nullptr
        );

        //- Check that the requested tangent can be provided at the integration
        //  points of the given topology for the current set of laws
        void checkTangentRequest
        (
            const integrationPointTopology& topo,
            const tangentRequest req
        );

        //- Create (if required) and return scratch containers for surface
        //  fields
        surfaceSymmTensorField& surfaceStressSum() const;
        surfaceScalarField& surfaceStressWeight() const;
        surfaceScalarField& surfaceTangentWeight() const;

        //- Create (if required) and return scratch containers for point
        //  fields
        pointSymmTensorField& pointStressSum() const;
        pointScalarField& pointStressWeight() const;
        pointScalarField& pointTangentWeight() const;

        //- Clear surface scratch fields
        inline void clearScratchFields() const
        {
            surfaceStressSumPtr_.clear();
            surfaceStressWeightPtr_.clear();
            surfaceTangentWeightPtr_.clear();

            pointStressSumPtr_.clear();
            pointStressWeightPtr_.clear();
            pointTangentWeightPtr_.clear();

            scratchStress_.clear();
        }

    //- Runtime type information
    TypeName("mechanicalConstitutiveLawManager");


    // Constructors

        //- Construct from a mesh and dict
        mechanicalConstitutiveLawManager
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~mechanicalConstitutiveLawManager();


    // Member Functions

    // Access

        //- Return (and lazily construct) a run-time selectable
        //  integration-point topology by type name.
        //  The topology is owned and cached by the manager and remains valid
        //  for the lifetime of the manager.
        //  Topologies that need more than the mesh to be built go through
        //  registerTopology() instead
        const integrationPointTopology& topologyFor
        (
            const word& topologyTypeName
        ) const;

        //- Const-access to the density field
        const volScalarField& rho() const;

        //- Const-access to the bulk modulus field
        const volScalarField& kappa() const;

    // Edit

        //- Register an externally constructed integration-point topology and
        //  return a reference to it.
        //  This is required for topologies that cannot be built from the mesh
        //  alone, and so cannot go through the run-time selection table, e.g.
        //  dualFaceIntegrationPointTopology, which also needs the
        //  dual-face-to-cell map.
        //  The manager takes ownership and the topology remains valid for the
        //  lifetime of the manager.
        //  The key shares a namespace with the run-time selectable topology
        //  type names, so it must not collide with one.
        //  The call is idempotent: registering a key that is already present
        //  returns the topology already held and discards the supplied one.
        //  This matters because constitutive state is keyed on the topology
        //  object, so replacing it would silently start a second, empty set of
        //  history variables
        const integrationPointTopology& registerTopology
        (
            const word& key,
            autoPtr<integrationPointTopology> topoPtr
        ) const;

        //- Reset cached material property fields (e.g. rho, bulk modulus).
        //  This does not modify constitutive laws or their state; it only
        //  clears cached volScalarField representations so they will be
        //  rebuilt on next access.
        void resetMaterialPropertyFields();

        //- Update the engineering (small-strain) stress at the integration
        //  points of the supplied topology.
        //  This is the primitive form of the small-strain update: the
        //  GeometricField and CompactListList overloads below are implemented
        //  in terms of it.
        //  The caller owns all storage and is responsible for its indexing,
        //  which must match the integration-point indices returned by the
        //  topology. The lists need not be associated with the manager's mesh,
        //  which is what allows evaluation on a dual or otherwise derived mesh.
        //  Only internal integration points are visited: boundary states are
        //  the business of the GeometricField overloads.
        //  No stress or tangent collapse is performed. A topology whose
        //  integration points are shared between cells is therefore allowed
        //  only with a single mechanical constitutive law, where no
        //  integration point can belong to two materials and there is nothing
        //  to collapse.
        //  A fourth-order tangent request additionally requires
        //  topo.supportsFourthOrderTangent().
        //  Unlike the GeometricField overloads, a tangent request that is not
        //  matched by the corresponding storage is an error rather than a
        //  silently skipped tangent
        void updateStressSmallStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& gradD,
            const UList<tensor>& gradD0,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr = nullptr,
            UList<mat66>* fourthOrderTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Finite-strain counterpart of the flat-list small-strain update
        //  above, with the same storage, indexing and tangent requirements
        void updateStressFiniteStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& F,
            const UList<tensor>& F0,
            const UList<tensor>& Finv,
            const UList<tensor>& Finv0,
            const UList<scalar>& J,
            const UList<scalar>& J0,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr = nullptr,
            UList<mat66>* fourthOrderTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Evaluate only the requested scalar tangent at cell centres,
        //  modifying nothing else.
        //  Each law is evaluated against a shadow of its constitutive state,
        //  so history is read but never written, and the stress goes to
        //  manager-owned scratch storage. Neither the caller's fields nor any
        //  history variable is disturbed.
        //  Intended for solid models that need an implicit stiffness
        //  coefficient at construction time, before any stress update has been
        //  performed, or at a chosen refresh cadence between updates
        void updateScalarTangent
        (
            const volTensorField& gradD,
            const volTensorField& gradD0,
            const scalar dt,
            volScalarField& scalarTangent,
            const tangentRequest tangentReq = tangentRequest::scalar,
            const bool coldState = false
        );

        //- Finite-strain counterpart of updateScalarTangent.
        //  Returns a scalar tangent at cell centres, evaluated from the
        //  deformation gradient rather than the displacement gradient, for a
        //  finite-strain solid model that needs an implicit stiffness
        //  coefficient. Like its small-strain counterpart it evaluates each
        //  law against a shadow of its state, so it neither writes a stress
        //  nor disturbs history
        void updateScalarTangentFiniteStrain
        (
            const volTensorField& F,
            const volTensorField& F0,
            const volTensorField& Finv,
            const volTensorField& Finv0,
            const volScalarField& J,
            const volScalarField& J0,
            const scalar dt,
            volScalarField& scalarTangent,
            const tangentRequest tangentReq = tangentRequest::scalar
        );

        //- Evaluate only the requested small-strain tangent at the integration
        //  points of the supplied topology, leaving the caller's stress
        //  storage untouched.
        //  This is for solid models that take their Jacobian from the manager
        //  while their residual stress still comes from elsewhere, which is the
        //  intended migration path.
        //  Each law is evaluated against a shadow of its constitutive state,
        //  so a tangent query never disturbs history, whichever topology it is
        //  made on.
        //  One qualification. Like every other update, this participates in
        //  the once-per-time-step old-time rollover, so if it is the first
        //  call of a new time step it is what commits the previous step's
        //  converged state. That is bookkeeping the framework owes whoever
        //  evaluates first, not something the query computed: shadows never
        //  write current-time fields, so the values committed are the same
        //  whichever call triggers the rollover
        void updateTangentSmallStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& gradD,
            const UList<tensor>& gradD0,
            const scalar dt,
            UList<scalar>* scalarTangentPtr,
            UList<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq,
            const bool coldState = false
        );

        //- Finite-strain counterpart of updateTangentSmallStrain
        void updateTangentFiniteStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& F,
            const UList<tensor>& F0,
            const UList<tensor>& Finv,
            const UList<tensor>& Finv0,
            const UList<scalar>& J,
            const UList<scalar>& J0,
            const scalar dt,
            UList<scalar>* scalarTangentPtr,
            UList<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq
        );

        //- Update the engineering (small-strain) stress field using the
        //  configured mechanical constitutive laws.
        //  The stress is evaluated from the current and previous displacement
        //  gradients and the supplied time increment, and is updated in-place.
        //  Constitutive history variables are stored and updated internally via
        //  mechanicalConstitutiveLawState objects.
        //  Optionally computes an approximate scalar material tangent when
        //  requested, e.g. for segregated or mixed displacement–pressure
        //  formulations.
        void updateStressSmallStrain
        (
            const volTensorField& gradD,
            const volTensorField& gradD0,
            const scalar dt,
            volSymmTensorField& stress,
            volScalarField* scalarTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Update the engineering (small-strain) stress field stored on mesh
        //  faces using the configured mechanical constitutive laws.
        //  Unlike cell-based approach, faces can appear in more than one
        //  material law (e.g., material interface), so a stress "collapse" rule
        //  must be defined to specify how the multiple stresses are combined
        void updateStressSmallStrain
        (
            const surfaceTensorField& gradD,
            const surfaceTensorField& gradD0,
            const scalar dt,
            surfaceSymmTensorField& stress,
            const stressCollapseRule collapseRule = stressCollapseRule::none,
            surfaceScalarField* scalarTangentPtr = nullptr,
            List<mat66>* fourthOrderTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Update the engineering (small-strain) stress field stored on mesh
        //  faces using the configured mechanical constitutive laws.
        //  Unlike cell-based approach, points can appear in more than one
        //  material law (e.g., material interface), so a stress "collapse" rule
        //  must be defined to specify how the multiple stresses are combined
        void updateStressSmallStrain
        (
            const pointTensorField& gradD,
            const pointTensorField& gradD0,
            const scalar dt,
            pointSymmTensorField& stress,
            const stressCollapseRule collapseRule = stressCollapseRule::none,
            pointScalarField* scalarTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Update the engineering (small-strain) stress values stored in a
        //  compact integration-point layout using the configured mechanical
        //  constitutive laws.
        //  Stress is evaluated from the current and previous displacement
        //  gradients and the supplied time increment, and written in-place to
        //  the packed integration-point storage.
        //  Constitutive history variables are updated internally via
        //  mechanicalConstitutiveLawState objects associated with each material.
        //  This interface is intended for higher-order or multi-point
        //  formulations, e.g. Gauss-point-based discretisations.
        //  Optionally computes an approximate scalar material tangent defined
        //  per integration point when requested.
        //  Requirement:
        //    - All CompactListList inputs must share identical offsets and total
        //      size.
        void updateStressSmallStrain
        (
            const CompactListList<tensor>& gradD,
            const CompactListList<tensor>& gradD0,
            const scalar dt,
            CompactListList<symmTensor>& stress,
            List<scalar>* scalarTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Update the true (Cauchy) stress field using the configured
        //  mechanical constitutive laws for a finite-strain formulation.
        //  The stress is evaluated from the current and previous deformation
        //  gradients, their inverses, determinants, and the supplied time
        //  increment, and is updated in-place.
        //  Constitutive history variables are stored and updated internally via
        //  mechanicalConstitutiveLawState objects.
        //  Optionally computes an approximate scalar material tangent when
        //  requested, e.g. for segregated or mixed displacement–pressure
        //  formulations
        void updateStressFiniteStrain
        (
            const volTensorField& F,
            const volTensorField& F0,
            const volScalarField& J,
            const volScalarField& J0,
            const volTensorField& Finv,
            const volTensorField& Finv0,
            const scalar dt,
            volSymmTensorField& stress,
            volScalarField* scalarTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Update the true (Cauchy) stress field stored on mesh faces using the
        //  configured mechanical constitutive laws for a finite-strain formulation.
        //  This interface is intended for  discretisations where stresses are
        //  evaluated at face-based integration points.
        //  The stress is evaluated from the current and previous deformation
        //  gradients, their inverses, determinants, and the supplied time
        //  increment, and is updated in-place.
        //  Constitutive history variables are stored and updated internally via
        //  mechanicalConstitutiveLawState objects.
        //  Optionally computes an approximate scalar material tangent when
        //  requested.
        // void updateStressFiniteStrain
        // (
        //     const surfaceTensorField& F,
        //     const surfaceTensorField& F0,
        //     const surfaceScalarField& J,
        //     const surfaceScalarField& J0,
        //     const surfaceTensorField& Finv,
        //     const surfaceTensorField& Finv0,
        //     const scalar dt,
        //     surfaceSymmTensorField& stress,
        //     surfaceScalarField* scalarTangentPtr = nullptr,
        //     const tangentRequest tangentReq = tangentRequest::none
        // );

        //- Update the true (Cauchy) stress values stored in a compact
        //  integration-point layout using the configured mechanical
        //  constitutive laws for a finite-strain formulation.
        //  Stress is evaluated from the current and previous deformation
        //  gradients, their inverses, determinants, and the supplied time
        //  increment, and written in-place to the packed integration-point
        //  storage.
        //  Constitutive history variables are updated internally via
        //  mechanicalConstitutiveLawState objects associated with each material.
        //  This interface is intended for higher-order or multi-point
        //  formulations, e.g. Gauss-point-based discretisations.
        //  Optionally computes an approximate scalar material tangent defined
        //  per integration point when requested.
        //  Requirement:
        //    - All CompactListList inputs must share identical offsets and
        //      total size.
        void updateStressFiniteStrain
        (
            const CompactListList<tensor>& F,
            const CompactListList<tensor>& F0,
            const CompactListList<tensor>& Finv,
            const CompactListList<tensor>& Finv0,
            const CompactListList<scalar>& J,
            const CompactListList<scalar>& J0,
            const scalar dt,
            CompactListList<symmTensor>& stress,
            List<scalar>* scalarTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Evaluate finite-strain stress, asking each law for its isochoric
        //  stress and volumetric response separately.
        //
        //  For a mixed displacement-pressure formulation, which replaces the
        //  volumetric part with a solved pressure and so needs the two apart.
        //  A law that cannot separate them is refused by name: taking the
        //  trace of a total stress cannot tell a volumetric response the
        //  formulation may replace from one it must keep
        void updateStressFiniteStrainSplit
        (
            const volTensorField& F,
            const volTensorField& F0,
            const volTensorField& Finv,
            const volTensorField& Finv0,
            const volScalarField& J,
            const volScalarField& J0,
            const scalar dt,
            volSymmTensorField& isochoricStress,
            volScalarField& volumetricResponse
        );

        //- Refuse if any law cannot separate the two, naming it
        void checkVolumetricSplitSupported(const word& context) const;

        //- Let each material law know the end of each time step has been
        //  reached.
        //  Can be used for diagnostics, logging, or state bookkeeping.
        void endTimeStep();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
