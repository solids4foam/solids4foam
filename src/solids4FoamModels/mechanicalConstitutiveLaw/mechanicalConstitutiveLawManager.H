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
#include "mechanicalConstitutiveLawStressCollapseRule.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class integrationPointTopology;
class CompactListListBase;

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

        //- Density field
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
            List<PtrList<mechanicalConstitutiveLawState>> boundaryStates_;

            //- True if this topology supports boundary integration points
            bool boundaryAware_;

            //- Constructor
            topologyEntry(const integrationPointTopology& topo)
            :
                topology_(topo),
                lawIntegrationPointIDs_(),
                states_(),
                boundaryStates_(),
                boundaryAware_(false)
            {}
        };

        //- Cache of topology-specific integration-point data and state
        //  Keyed by the address of the integrationPointTopology instance
        mutable HashTable<topologyEntry> topologyEntries_;


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

        //- Check the provided compact list sizes are consistent
        inline void checkCompactLayoutConsistency
        (
            const CompactListList<tensor>& a,
            const CompactListList<tensor>& b,
            const CompactListList<symmTensor>& out,
            const List<scalar>* tangentPtr,
            const word& context
        )
        {
            const label nIP = out.totalSize();

            if (a.totalSize() != nIP || b.totalSize() != nIP)
            {
                FatalError
                    << "Inconsistent CompactListList sizes in " << context << nl
                    << "Expected totalSize = " << nIP << nl
                    << "Got: grad = " << a.totalSize()
                    << ", grad0 = " << b.totalSize()
                    << exit(FatalError);
            }

            if (tangentPtr && tangentPtr->size() != nIP)
            {
                FatalError
                    << "Scalar tangent list has incorrect size in " << context
                    << nl
                    << "Expected: " << nIP
                    << ", got: " << tangentPtr->size()
                    << exit(FatalError);
            }
        }

        //- Return (and lazily construct) an integration-point topology by type
        //  name
        const integrationPointTopology& topologyFor
        (
            const word& topologyTypeName
        ) const;

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
        void updateOldTimeIfNeeded();

        //- Check if the provided tangent request is scalar or scalarDeviatoric
        inline bool needsScalarTangent(const tangentRequest req)
        {
            return
                req == tangentRequest::scalar
             || req == tangentRequest::scalarDeviatoric;
        }

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
        }

        //- No copy construct
        mechanicalConstitutiveLawManager
        (
            const mechanicalConstitutiveLawManager&
        ) = delete;

        //- No copy assignment
        void operator=(const mechanicalConstitutiveLawManager&) = delete;


public:

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

        //- Const-access to the density field
        const volScalarField& rho() const;

        //- Const-access to the bulk modulus field
        const volScalarField& kappa() const;

    // Edit

        //- Reset cached material property fields (e.g. rho, bulk modulus).
        //  This does not modify constitutive laws or their state; it only
        //  clears cached volScalarField representations so they will be
        //  rebuilt on next access.
        void resetMaterialPropertyFields();

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
