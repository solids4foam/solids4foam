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
#include "integrationPointTopology.H"
#include "mechanicalConstitutiveLawCaseInputs.H"
#include "mechanicalConstitutiveLawStateSetup.H"
#include "mechanicalConstitutiveLawInputGatherer.H"

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
struct smallStrainKinematicsFields;
struct finiteStrainKinematicsFields;

//- OpenFOAM.org's CompactListList has two template parameters, so its header
//  is included above. The other forks have one template parameter; forward
//  declare it to avoid instantiating foam-extend's broken const operator[].
#ifndef OPENFOAM_ORG
template<class T> class CompactListList;
#endif

/*---------------------------------------------------------------------------*\
              Class mechanicalConstitutiveLawManager Declaration
\*---------------------------------------------------------------------------*/

class mechanicalConstitutiveLawManager
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Pointer to the point mesh, if needed
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

        //- The optional coupling inputs already reported as registered but
        //  unread, so that each is reported once rather than per evaluation
        mutable HashSet<word> reportedUnreadInputs_;

        //- Cell indices defining the material support for each mechanical
        //  constitutive law.
        //  Materials in solids4foam are always defined per cell (via
        //  cellZones), independent of the discretisation order or the
        //  locations at which stresses or strains are evaluated.
        //  These lists are used to construct integration-point mappings and
        //  to assemble cell-based material properties (e.g. density, bulk
        //  modulus).
        List<labelList> lawCells_;

        //- Boundary patch face indices per mechanical constitutive law.
        //  For each mechanicalConstitutiveLaw, this stores a list over
        //  boundary patches, where each entry contains the local face indices
        //  (within that patch) whose owner cell belongs to the corresponding
        //  law.
        //  In other words:
        //    lawBoundaryFaces_[lawI][patchI] gives the list of boundary faces
        //    on patch patchI whose adjacent cell is assigned to law lawI.
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

        //- Mechanical constitutive law index of each cell.
        //  Needed by the flat-list boundary update for patches that carry no
        //  fvPatch faces, where lawBoundaryFaces_ cannot supply the law
        labelList cellToLaw_;

        //- Unregistered material density used to initialise solver fields
        mutable autoPtr<volScalarField> rhoPtr_;

        //- Bulk modulus field
        mutable autoPtr<volScalarField> kappaPtr_;

        //- Cache of integration-point topologies, owned by the manager
        mutable HashTable<autoPtr<integrationPointTopology>> topologyCache_;

        //- Topology-specific constitutive data owned by the manager.
        //  A topologyEntry groups together all data that depends on the
        //  *location* of integration points (e.g. cell-centred, face-centred,
        //  vertex-centred, or higher-order multi-point layouts).
        //  For a given integrationPointTopology, this includes:
        //    - the mapping from mechanical laws to integration-point indices
        //      (lawIntegrationPointIDs),
        //    - the constitutive state objects storing history variables at
        //      those integration points, and
        //    - optional boundary constitutive states where applicable.
        //  This indirection allows a single mechanicalConstitutiveLawManager to
        //  evaluate stresses at multiple locations (cells, faces, points, etc.)
        //  using the same material definitions, while keeping state variables
        //  strictly separated per topology.
        //  Topology entries are constructed lazily on first use and cached
        //  for the lifetime of the manager.
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

            //- Each law's convergence scale for the current evaluation,
            //  reduced over every rank. Held here so that the boundary
            //  evaluations use the same scale as the internal ones rather
            //  than one taken over whichever faces this rank happens to hold
            scalarList lawConvergenceScales_;

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
                lawConvergenceScales_(),
                lawScalarInputs_()
            {}
        };

        //- Cache of topology-specific integration-point data and state
        //  Keyed as topologyCache_ is, by the name the topology was
        //  registered under
        mutable HashTable<autoPtr<topologyEntry>> topologyEntries_;

        //- Shape of the compact layout each compact topology was built
        //  from, keyed as topologyCache_ is.
        //  The key names the role and is identical on every rank; this
        //  records the addressing behind it, so a second layout arriving
        //  under the same role is refused rather than silently given the
        //  first one's history
        mutable HashTable<word> compactFingerprints_;

        //- Signature of the mesh the per-patch addressing above was built
        //  for: the number of cells, the number of faces on each patch, and a
        //  checksum of each patch's face-to-cell addressing.
        //  A topology-changing mesh (e.g. crackerFvMesh, which moves internal
        //  faces onto a crack patch) changes these under the manager;
        //  comparing against them is how that is detected. The checksums
        //  catch a change that leaves every size as it was
        labelList addressingMeshSizes_;

        //- Scalar inputs read from another case directory, for the laws
        //  whose dictionaries ask for that; empty if none does
        autoPtr<mechanicalConstitutiveLawCaseInputs> caseInputsPtr_;

        //- Gathers the coupling inputs each law reads, from the registry, a
        //  case directory or a file. Holds references to the members above,
        //  so it is declared after them
        mechanicalConstitutiveLawInputGatherer inputGatherer_;

        //- Whether a continued run has the kinematic history a constitutive
        //  history is measured against. Said by the solid model, which knows
        //  what it wrote; true until it says otherwise
        bool restartKinematicsAvailable_;

        //- Sets up each law's state: declared defaults, prescribed fields
        //  and restart. Holds references to the members above, so it is
        //  declared after them
        mechanicalConstitutiveLawStateSetup stateSetup_;


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
                    << "Expected mesh: " << expectedMesh.name() << nl
                    << "Field mesh:    " << fieldMesh.name()
                    << exit(FatalError);
            }
        }

        //- Check that two compact layouts have the same rows, not merely the
        //  same number of values.
        //
        //  Equal flattened lengths are not enough: row sizes of [2, 1] and
        //  [1, 2] both hold three values, and flat index 1 is then the second
        //  point of the first cell in one and the first point of the second
        //  cell in the other. Nothing downstream can notice, because every
        //  list is the length it should be and the addressing is implied by
        //  the offsets rather than carried with the values
        static void checkCompactRowSizes
        (
            const labelList& ref,
            const word& refName,
            const labelList& other,
            const word& otherName,
            const word& context
        );

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

        //- The flat-list evaluation, written once for both strain measures.
        //  Fields is smallStrainKinematicsFields or
        //  finiteStrainKinematicsFields; evaluateSmallStrain and
        //  evaluateFiniteStrain forward to it
        template<class Fields>
        void evaluateFlat
        (
            const integrationPointTopology& topo,
            const Fields& fields,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr,
            UList<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq,
            const bool preserveState,
            const bool coldState,
            UList<scalar>* volumetricPtr
        );

        //- The boundary faces of the volField stress updates, written once
        //  for both strain measures. patchFields(patchI) gives the patch's
        //  kinematic lists as Fields; lead is the field whose patches say
        //  which are coupled
        template<class Fields, class PatchFieldsFn>
        void updateStressVolBoundary
        (
            topologyEntry& tp,
            const volTensorField& lead,
            const PatchFieldsFn& patchFields,
            const scalar dt,
            volSymmTensorField& stress,
            volScalarField* scalarTangentPtr,
            const tangentRequest tangentReq,
            volScalarField* volumetricResponsePtr
        );

        //- The surfaceField stress updates, written once for both strain
        //  measures: each law's internal faces and its boundary faces are
        //  evaluated, and a face reached by more than one law collapsed by
        //  collapseRule. patchFields(patchI) gives the patch's kinematic
        //  lists as Fields. A null tangent or volumetric pointer is not
        //  computed
        template<class Fields, class PatchFieldsFn>
        void updateStressSurface
        (
            const integrationPointTopology& topo,
            topologyEntry& tp,
            const Fields& fields,
            const PatchFieldsFn& patchFields,
            const scalar dt,
            surfaceSymmTensorField& stress,
            const stressCollapseRule collapseRule,
            surfaceScalarField* scalarTangentPtr,
            List<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq,
            surfaceScalarField* volumetricResponsePtr
        );

        //- Check each kinematic list has the size the topology requires
        void checkKinematicsListSizes
        (
            const label nIP,
            const smallStrainKinematicsFields& fields,
            const word& context
        ) const;

        void checkKinematicsListSizes
        (
            const label nIP,
            const finiteStrainKinematicsFields& fields,
            const word& context
        ) const;

        //- Check that an integration point reached by more than one
        //  constitutive law has been given a rule for combining them.
        //  stressCollapseRule::none says there is no such rule, so a second
        //  contribution is an error rather than something to average away
        inline void checkCollapsePermitted
        (
            const stressCollapseRule collapseRule,
            const scalar weight,
            const label ipI,
            const word& what
        ) const
        {
            if (collapseRule != stressCollapseRule::none || weight < 1.5)
            {
                return;
            }

            FatalErrorInFunction
                << what << ' ' << ipI << " was reached by "
                << label(weight) << " mechanical constitutive laws, but the "
                << "stress collapse rule is '"
                << stressCollapseRuleName(collapseRule)
                << "', which admits only one." << nl
                << "An integration point shared between materials has as many "
                << "stresses as it has materials, and no rule has been given "
                << "for combining them." << nl
                << "Pass stressCollapseRule::average or "
                << "stressCollapseRule::harmonic to say how they combine."
                << exit(FatalError);
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

        //- Each law's convergence scale, reduced over every rank.
        //  Asked of every law whether or not this rank holds any of its
        //  points, so that the number of reductions is the same everywhere.
        //  Fields is smallStrainKinematicsFields or
        //  finiteStrainKinematicsFields
        template<class Fields>
        scalarList convergenceScales
        (
            topologyEntry& tp,
            const Fields& fields
        ) const;

        //- The key a topology is registered under in topologyCache_.
        //
        //  Identity is the registration key, not the object's address. The
        //  address is not identity: it is reused once an object dies, and it
        //  differs between ranks, so anything ordered by it - the collective
        //  reduction in endTimeStep(), the restart file names - would order
        //  differently on each rank and could pair one topology's history with
        //  another's after a free and reallocate.
        //
        //  The address is still the right way to *find* the entry, because a
        //  registered topology lives as long as the manager: comparing
        //  addresses here is a membership test, which is sound, rather than a
        //  persistent name, which is not.
        //
        //  Fatal if the topology is not one this manager owns. Every topology
        //  reaches the manager through topologyFor(), registerTopology() or
        //  compactCellTopologyFor(), all of which register it, so a topology
        //  that is not in the cache is a caller holding an object the manager
        //  never took ownership of, and its state would be silently unrelated
        //  to any other
        const word& topologyKeyFor(const integrationPointTopology& topo) const;

        //- Ensure a topology entry exists and return it
        topologyEntry& topology
        (
            const integrationPointTopology& topo
        ) const;

        //- As above, by registration key
        topologyEntry& topology(const word& key) const;

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

        //- The coupling-input gathering below is done by inputGatherer_; see
        //  mechanicalConstitutiveLawInputGatherer for what each one does

        mechanicalConstitutiveLawInputs lawInputs
        (
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            const scalar dt,
            topologyEntry& tp
        ) const
        {
            return inputGatherer_.lawInputs
            (
                lawI, topo, ipIDs, dt, tp.lawScalarInputs_
            );
        }

        mechanicalConstitutiveLawInputs lawInputsPatch
        (
            const label lawI,
            const label patchI,
            const labelList& faces,
            const scalar dt,
            topologyEntry& tp
        ) const
        {
            return inputGatherer_.lawInputsPatch
            (
                lawI, patchI, faces, dt, tp.lawScalarInputs_
            );
        }

        void refreshScalarInputs() const
        {
            inputGatherer_.refreshScalarInputs();
        }

        static wordList requiredScalarInputsRecursive
        (
            const mechanicalConstitutiveLaw& law
        )
        {
            return
                mechanicalConstitutiveLawInputGatherer::
                requiredScalarInputsRecursive(law);
        }

        //- Tell a law and every law below it that the step has ended
        void endTimeStepLaw
        (
            const mechanicalConstitutiveLaw& law,
            mechanicalConstitutiveLawState& state,
            const scalar time,
            const label timeIndex,
            DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics,
            DynamicList<const mechanicalConstitutiveLaw*>* reporters = nullptr,
            DynamicList<label>* firsts = nullptr,
            DynamicList<label>* counts = nullptr
        ) const;

        //- Setting up a law's state - declared defaults, prescribed fields
        //  and restart - is done by stateSetup_; see
        //  mechanicalConstitutiveLawStateSetup for what each one does

        bool declaresPersistentState
        (
            const mechanicalConstitutiveLaw& law
        ) const
        {
            return stateSetup_.declaresPersistentState(law);
        }

        void setupStateRestart
        (
            topologyEntry& entry,
            const integrationPointTopology& topo,
            const word& topologyKey
        ) const
        {
            stateSetup_.setupStateRestart
            (
                entry.lawIntegrationPointIDs_,
                entry.states_,
                entry.boundaryStates_,
                entry.boundaryAware_,
                topo,
                topologyKey
            );
        }

        void applyStateSpec
        (
            const label lawI,
            const integrationPointTopology& topo,
            const labelList& ipIDs,
            mechanicalConstitutiveLawState& state
        ) const
        {
            stateSetup_.applyStateSpec(lawI, topo, ipIDs, state);
        }

        void applyStateSpecScratch
        (
            const label lawI,
            mechanicalConstitutiveLawState& state
        ) const
        {
            stateSetup_.applyStateSpecScratch(lawI, state);
        }

        void applyStateSpecPatch
        (
            const label lawI,
            const label patchI,
            mechanicalConstitutiveLawState& state
        ) const
        {
            stateSetup_.applyStateSpecPatch(lawI, patchI, state);
        }

        //- Roll every law's state to old time, once, on the first
        //  evaluation of a new time step. Called from each entry point
        //  rather than by the solid model, so that a law's history advances
        //  exactly once per step however many times its stress is evaluated
        void updateOldTimeIfNeeded();

        //- The current mesh signature, in the layout of addressingMeshSizes_
        labelList currentMeshSizes() const;

        //- Build lawBoundaryFaces_ from cellToLaw_ for the current mesh
        void calcLawBoundaryFaces();

        //- Rebuild the per-patch addressing and boundary states if the mesh
        //  topology has changed since they were built. Checked on every
        //  evaluation, not once per time step, because a cracking mesh
        //  changes topology within a time step. Cell addressing and cell
        //  states are kept, so a change in the number of cells is refused,
        //  as is any change when a law carries persistent state, which would
        //  need mapping onto the new boundary faces
        void updateAddressingIfTopologyChanged();


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
            const bool coldState = false,
            UList<scalar>* volumetricPtr = nullptr
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
        //- Scratch accumulator for collapsing scalar tangents onto faces,
        //  allocated on first use and reused thereafter
        surfaceScalarField& surfaceTangentWeight() const;

        //- Create (if required) and return scratch containers for point
        //  fields
        pointSymmTensorField& pointStressSum() const;
        pointScalarField& pointStressWeight() const;
        //- As above, for point-centred topologies
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

        //- Const-access to the unregistered material-density fallback field
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
        //  It must also be the same word on every rank. The key is what
        //  names this topology's state and orders the end-of-step
        //  reduction, so a key built from something rank-local - a cell
        //  count, a local digest - would have the ranks reduce in
        //  different orders. endTimeStep() checks this and stops, but it
        //  cannot repair it.
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
            const tangentRequest tangentReq = tangentRequest::none,
            volScalarField* volumetricResponsePtr = nullptr
        );

        //- Update the engineering (small-strain) stress field stored on mesh
        //  faces using the configured mechanical constitutive laws.
        //  Unlike the cell-based approach, a face can belong to more than one
        //  material law (e.g. at a material interface), so a stress "collapse"
        //  rule says how the multiple contributions are combined.
        //  The default, stressCollapseRule::none, is a refusal rather than a
        //  rule: a face reached by two laws is then an error
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
        //  points using the configured mechanical constitutive laws.
        //  Unlike the cell-based approach, a point can belong to more than one
        //  material law (e.g. at a material interface), so a stress "collapse"
        //  rule says how the multiple contributions are combined.
        //  The default, stressCollapseRule::none, is a refusal rather than a
        //  rule: a point reached by two laws is then an error
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
        //  mechanicalConstitutiveLawState objects associated with each
        //  material.
        //  This interface is intended for higher-order or multi-point
        //  formulations, e.g. Gauss-point-based discretisations.
        //  Optionally computes an approximate scalar material tangent defined
        //  per integration point when requested.
        //  Requirement:
        //    - All CompactListList inputs must share identical offsets and
        //      total size.
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
            const tangentRequest tangentReq = tangentRequest::none,
            volScalarField* volumetricResponsePtr = nullptr
        );

        //- Update the true (Cauchy) stress values stored in a compact
        //  integration-point layout using the configured mechanical
        //  constitutive laws for a finite-strain formulation.
        //  Stress is evaluated from the current and previous deformation
        //  gradients, their inverses, determinants, and the supplied time
        //  increment, and written in-place to the packed integration-point
        //  storage.
        //  Constitutive history variables are updated internally via
        //  mechanicalConstitutiveLawState objects associated with each
        //  material.
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

        //- Evaluate small-strain stress, asking each law for its deviatoric
        //  stress and volumetric response separately.
        //
        //  The small-strain volumetric response is K*tr(epsilon), which is
        //  what the mixed formulation's pressure equation solves for, exactly
        //  as dU/dJ is at finite strain. A law that cannot separate the two is
        //  refused by name: the trace of a total stress is K*tr(epsilon) only
        //  for an isotropic law, and for an anisotropic one the difference is
        //  a different material rather than a visible error
        void updateStressSmallStrainSplit
        (
            const volTensorField& gradD,
            const volTensorField& gradD0,
            const scalar dt,
            volSymmTensorField& deviatoricStress,
            volScalarField& volumetricResponse
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

        //- Update the true (Cauchy) stress field stored on mesh faces using
        //  the configured mechanical constitutive laws for a finite-strain
        //  formulation: the face twin of the cell-centred overload, with the
        //  same collapse rule as the small-strain face overload for a face
        //  reached by more than one law. When volumetricResponsePtr is given,
        //  stress receives the isochoric stress and the volumetric response
        //  is written separately, as updateStressFiniteStrainSplit does
        void updateStressFiniteStrain
        (
            const surfaceTensorField& F,
            const surfaceTensorField& F0,
            const surfaceScalarField& J,
            const surfaceScalarField& J0,
            const surfaceTensorField& Finv,
            const surfaceTensorField& Finv0,
            const scalar dt,
            surfaceSymmTensorField& stress,
            const stressCollapseRule collapseRule = stressCollapseRule::none,
            surfaceScalarField* volumetricResponsePtr = nullptr
        );

        //- The face twin of updateStressFiniteStrainSplit
        void updateStressFiniteStrainSplit
        (
            const surfaceTensorField& F,
            const surfaceTensorField& F0,
            const surfaceTensorField& Finv,
            const surfaceTensorField& Finv0,
            const surfaceScalarField& J,
            const surfaceScalarField& J0,
            const scalar dt,
            surfaceSymmTensorField& isochoricStress,
            surfaceScalarField& volumetricResponse,
            const stressCollapseRule collapseRule = stressCollapseRule::none
        );

        //- Refuse if any law cannot separate the two, naming it
        void checkVolumetricSplitSupported(const word& context) const;

        //- True if every law can separate its isochoric stress from its
        //  volumetric response. Asked rather than assumed by anything that
        //  wants to exercise the split without demanding it
        //- How many materials this manager holds
        label nLaws() const
        {
            return laws_.size();
        }

        //- The one law, for a solid model that is single-material by
        //  construction.
        //
        //  A Kirchhoff plate has one thickness and one bending stiffness, and
        //  a block-coupled linear elastic solver assembles one set of
        //  constants into its matrix; both are isotropic single-material
        //  models by their own formulation rather than by limitation. Such a
        //  model needs the law itself, to read material constants from, not a
        //  field. Refuses rather than picks the first when there is more than
        //  one, because a model asking this has no meaning for a second
        //  material
        const mechanicalConstitutiveLaw& singleLaw() const;

        //- True if every law can separate its isochoric stress from its
        //  volumetric response. Asked rather than assumed by anything that
        //  wants to exercise the split without demanding it - the unit test
        //  does; a mixed formulation calls
        //  checkVolumetricSplitSupported() instead and refuses by name
        bool allLawsProvideVolumetricSplit() const;

        //- True if any law, or any law one of them wraps, is fully
        //  incompressible, so that it has no total stress
        bool anyLawIncompressible() const;

        //- True if every law's split stress is a function of the isochoric
        //  deformation alone, so that J*sigma_iso is invariant under a
        //  superposed dilation. Weaker than the above, and separate from it
        bool allLawsHaveDilationInvariantIsochoricStress() const;

        //- Say whether a continued run has the kinematic history that the
        //  laws' constitutive history is measured against. A client that
        //  continues a run without it must say so, and a law that keeps
        //  history then refuses to restart
        void setRestartKinematicsAvailable(const bool available)
        {
            restartKinematicsAvailable_ = available;
        }

        //- Let each material law know the end of each time step has been
        //  reached.
        //  Can be used for diagnostics, logging, or state bookkeeping.
        void endTimeStep();

        //- Write the cell-centred history of every law as volFields, for
        //  viewing.
        //
        //  The restart files hold the state as it is kept, which ParaView
        //  cannot show. These are a projection of it, written and never read:
        //  one field per persistent state variable, named after the variable
        //  and covering every material that declares it, zero in the cells of
        //  the others. A variable a composite and its sub-law both declare is
        //  qualified with the sub-law's name, and one whose name another
        //  model has registered is not written, rather than overwrite that
        //  model's field. Only the cell-centred topology is written: its
        //  points are the cells, so nothing is interpolated. The fields carry
        //  no dimensions, since a law does not declare them
        void writeStateFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
