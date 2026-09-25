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

\*---------------------------------------------------------------------------*/

#include "mechanicalConstitutiveLawManager.H"
#include "mechanicalConstitutiveLawStateIO.H"
#include "mechanicalConstitutiveLawEvaluateResponse.H"
#include "mechanicalConstitutiveLawKinematicsFields.H"
#include "IFstream.H"
#include "labelIOList.H"
#include "compatibilityFunctions.H"
#include "integrationPointTopologies.H"
#include "emptyFvPatch.H"
#include "calculatedFvPatchFields.H"
#include "mat66.H"
#include "Switch.H"
#include "CompactListList.H"
#include <cstdint>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalConstitutiveLawManager, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// A 64-bit FNV-1a digest of a list of row sizes, used to key a compact
// integration-point topology on the shape of the layout that produced it.
// Spelled out rather than taken from a fork's hasher, which the three
// supported forks do not agree on
std::uint64_t rowSizesHash(const labelUList& rowSizes)
{
    std::uint64_t h = 14695981039346656037ULL;

    forAll(rowSizes, i)
    {
        const std::uint64_t v = static_cast<std::uint64_t>(rowSizes[i]);

        for (int byteI = 0; byteI < 8; ++byteI)
        {
            h ^= (v >> (8*byteI)) & 0xFFULL;
            h *= 1099511628211ULL;
        }
    }

    return h;
}


// A checksum of an ordered list of keys, folded into a label so it can be
// compared across ranks with the reductions every fork provides
label keyListChecksum(const wordList& keys)
{
    std::uint64_t h = 14695981039346656037ULL;

    forAll(keys, i)
    {
        const word& k = keys[i];

        for (label c = 0; c < k.size(); ++c)
        {
            h ^= static_cast<std::uint64_t>(static_cast<unsigned char>(k[c]));
            h *= 1099511628211ULL;
        }

        // Separator, so {"ab","c"} and {"a","bc"} do not agree
        h ^= 0xFFULL;
        h *= 1099511628211ULL;
    }

    return static_cast<label>(h & 0x7FFFFFFFULL);
}


// Combine one diagnostic into another, by the operation it carries
void combineDiagnostic
(
    mechanicalConstitutiveLawDiagnostic& into,
    const mechanicalConstitutiveLawDiagnostic& from
)
{
    typedef mechanicalConstitutiveLawDiagnostic diagnostic;

    if (into.op() == diagnostic::combineOperation::sum)
    {
        into.value() += from.value();
    }
    else
    {
        into.value() = max(into.value(), from.value());
    }
}


//- Warn about the entries of a law's dictionary, or of one below it, that only
//  the removed legacy laws read
static void reportRemovedLegacyEntries
(
    const dictionary& dict,
    const word& lawName
)
{
    static const char* removed[] =
    {
        "pressureDisplacement",
        "pressureDisplacementCoeff",
        "alternatePressureDefinition",
        "impKcoeff",
        "calculateStressInLocalCoordinateSystem",
        "writeS0N0R",
        "tangentEps",
        "regionName",
        "pressureFieldRegion",
        "solvePressureEquation",
        "maxDeltaErr",
        "writeSubMeshes"
    };

    DynamicList<word> found;

    const label nRemoved = sizeof(removed)/sizeof(removed[0]);

    for (label i = 0; i < nRemoved; ++i)
    {
        const word name(removed[i]);

        if (dict.found(name))
        {
            found.append(name);
        }
    }

    if (found.size())
    {
        WarningInFunction
            << "Mechanical law '" << lawName << "' sets "
            << wordList(found) << ", which only the removed legacy "
            << "mechanical model read. They are ignored and can be removed."
            << endl;
    }

    const wordList keys(dict.toc());

    forAll(keys, i)
    {
        if (dict.isDict(keys[i]))
        {
            reportRemovedLegacyEntries
            (
                dict.subDict(keys[i]), lawName + '.' + keys[i]
            );
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void Foam::mechanicalConstitutiveLawManager::checkCompactRowSizes
(
    const labelList& ref,
    const word& refName,
    const labelList& other,
    const word& otherName,
    const word& context
)
{
    if (other.size() != ref.size())
    {
        FatalErrorInFunction
            << "Inconsistent compact layouts in " << context << nl << nl
            << "    " << refName << " has " << ref.size() << " rows and "
            << otherName << " has " << other.size() << '.' << nl << nl
            << "    The rows are the integration points of one mesh entity, "
            << "so two layouts with different row counts describe different "
            << "meshes."
            << exit(FatalError);
    }

    forAll(ref, rowI)
    {
        if (other[rowI] != ref[rowI])
        {
            FatalErrorInFunction
                << "Inconsistent compact layouts in " << context << nl << nl
                << "    Row " << rowI << " holds " << ref[rowI]
                << " integration points in " << refName << " and "
                << other[rowI] << " in " << otherName << '.' << nl << nl
                << "    The two hold the same number of values in total, so "
                << "every list is the length it should be, but the flat index "
                << "of a point differs between them: values from one entity "
                << "would be read as another's."
                << exit(FatalError);
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::checkCompactLayoutConsistency
(
    const CompactListList<tensor>& a,
    const CompactListList<tensor>& b,
    const CompactListList<symmTensor>& out,
    const List<scalar>* tangentPtr,
    const word& context
)
{
    const label nIP = out.m().size();

    if (a.m().size() != nIP || b.m().size() != nIP)
    {
        FatalError
            << "Inconsistent CompactListList sizes in " << context << nl
            << "Expected size = " << nIP << nl
            << "Got: grad = " << a.m().size()
            << ", grad0 = " << b.m().size()
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

    // The lengths agreeing does not make the layouts the same
    const labelList refRows(out.sizes());

    checkCompactRowSizes(refRows, "the stress", a.sizes(), "grad", context);
    checkCompactRowSizes(refRows, "the stress", b.sizes(), "grad0", context);
}


void Foam::mechanicalConstitutiveLawManager::checkTangentStorage
(
    const bool haveScalarTangent,
    const bool haveFourthOrderTangent,
    const tangentRequest req,
    const word& context
)
{
    if (needsScalarTangent(req) && !haveScalarTangent)
    {
        FatalErrorInFunction
            << "A " << tangentRequestName(req) << " tangent was requested in "
            << context << " but no scalar tangent storage was supplied."
            << exit(FatalError);
    }

    if (needsFourthOrderTangent(req) && !haveFourthOrderTangent)
    {
        FatalErrorInFunction
            << "A " << tangentRequestName(req) << " tangent was requested in "
            << context << " but no fourth-order tangent storage was supplied."
            << exit(FatalError);
    }
}


Foam::List<Foam::symmTensor>&
Foam::mechanicalConstitutiveLawManager::scratchStress
(
    const label nIntegrationPoints
) const
{
    if (scratchStress_.size() != nIntegrationPoints)
    {
        scratchStress_.setSize(nIntegrationPoints);
    }

    return scratchStress_;
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::topologyFor
(
    const word& topologyTypeName
) const
{
    // Already constructed?
    if (topologyCache_.found(topologyTypeName))
    {
        return topology(autoPtrRef(topologyCache_[topologyTypeName])).topology_;
    }

    // Lazily construct via OpenFOAM runtime selection
    autoPtr<integrationPointTopology> topoPtr
    (
        integrationPointTopology::New(topologyTypeName, mesh_)
    );

    if (!topoPtr.valid())
    {
        FatalErrorInFunction
            << "Failed to construct integrationPointTopology of type "
            << topologyTypeName
            << exit(FatalError);
    }

    // Cache and return
    topologyCache_.insert(topologyTypeName, topoPtr);

    return topology(autoPtrRef(topologyCache_[topologyTypeName])).topology_;
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::compactCellTopologyFor
(
    const CompactListList<tensor>& layout
) const
{
    // OpenFOAM.org keeps size() in a member of its own, which a layout filled
    // in through offsets() and m() - as the quadrature layouts are - leaves at
    // zero, so there the rows are counted from the offset table, which holds
    // one more entry than there are rows. The other forks derive size() from
    // their offsets, and foam-extend's hold row ends rather than row starts,
    // so they keep size(), sizes() and index()
#ifdef OPENFOAM_ORG
    const labelUList& offsets = layout.offsets();
    const label nRows = max(offsets.size() - 1, label(0));

    labelList rowSizes(nRows);
    forAll(rowSizes, rowI)
    {
        rowSizes[rowI] = offsets[rowI + 1] - offsets[rowI];
    }
#else
    const label nRows = layout.size();
    const labelList rowSizes(layout.sizes());
#endif

    // Which entity indexes the rows is decided by the row count, and checked.
    // A mesh never has as many cells as faces, so the two cases cannot be
    // confused, and anything else is an error rather than a guess
    const bool cellBased = (nRows == mesh_.nCells());
    const bool faceBased = (nRows == mesh_.nFaces());

    if (!cellBased && !faceBased)
    {
        FatalErrorInFunction
            << "A compact integration-point layout must have one row per cell "
            << "or one row per face." << nl
            << "This one has " << nRows << " rows, while the mesh has "
            << mesh_.nCells() << " cells and " << mesh_.nFaces() << " faces."
            << exit(FatalError);
    }

    // The key names the role, and is a literal, so it is the same word on
    // every rank. That matters because endTimeStep() sorts these keys and
    // reduces per entry: a key built from this rank's cell count or from a
    // digest of this rank's row sizes would sort differently on each rank and
    // pair one topology's quantity against another's.
    //
    // The shape is not identity, then, but it is still worth checking. The
    // digest below is a local fingerprint: it says whether the layout handed
    // in now has the same addressing as the one this topology was built from.
    // One compact cell layout and one compact face layout are supported per
    // manager, which is what the flat-list API uses; a second layout of the
    // same kind is a different topology and must be registered by name
    const word key(cellBased ? "compactCell" : "compactFace");

    const word fingerprint =
        Foam::name(rowSizes.size()) + ":"
      + Foam::name(layout.m().size()) + ":"
      + Foam::name(rowSizesHash(rowSizes));

    // Already constructed?
    if (topologyCache_.found(key))
    {
        if (!compactFingerprints_.found(key))
        {
            FatalErrorInFunction
                << "The name '" << key << "' is already registered, but not "
                << "by a compact layout." << nl << nl
                << "    That name is reserved for the topology the flat-list "
                << "CompactListList interface builds. Register other "
                << "topologies under a name of their own."
                << exit(FatalError);
        }

        const word& seen = compactFingerprints_[key];

        if (seen != fingerprint)
        {
            FatalErrorInFunction
                << "A second compact integration-point layout was passed to "
                << "the manager under the role '" << key << "'." << nl << nl
                << "    The layout first seen had shape " << seen
                << " and this one has " << fingerprint
                << " (rows:values:digest)." << nl << nl
                << "    The constitutive state is held per topology, so "
                << "returning the first topology for the second layout would "
                << "read one layout's history through the other's addressing."
                << nl << nl
                << "    Register the second layout under its own name with "
                << "registerTopology(), choosing a name that every rank "
                << "supplies identically."
                << exit(FatalError);
        }

        return topology(key).topology_;
    }

    // Construct topology lazily

    // Build cell → IP addressing
    // rowSizes rather than layout[cellI].size(): the const operator[] does
    // not compile on foam-extend
    CompactListList<label> cellToIP(rowSizes);

    for (label cellI = 0; cellI < nRows; ++cellI)
    {
        const label n = rowSizes[cellI];
        for (label j = 0; j < n; ++j)
        {
#ifdef OPENFOAM_ORG
            cellToIP(cellI, j) = offsets[cellI] + j;
#else
            cellToIP(cellI, j) = layout.index(cellI, j);
#endif
        }
    }

    autoPtr<integrationPointTopology> topoPtr;

    if (cellBased)
    {
        topoPtr.set
        (
            new compactCellIntegrationPointTopology(mesh_, std::move(cellToIP))
        );
    }
    else
    {
        topoPtr.set
        (
            new compactFaceIntegrationPointTopology(mesh_, std::move(cellToIP))
        );
    }

    // Cache and return
    topologyCache_.insert(key, topoPtr);
    compactFingerprints_.set(key, fingerprint);

    return topology(key).topology_;
}


template<class Fields>
Foam::scalarList Foam::mechanicalConstitutiveLawManager::convergenceScales
(
    topologyEntry& tp,
    const Fields& fields
) const
{
    scalarList scales(laws_.size(), 0.0);

    // Every law, on every rank, whether or not this rank holds any of its
    // points. That is the point of doing it here: the evaluation loops skip a
    // law with no points, so a law reducing for itself would reduce a
    // different number of times on each rank
    forAll(laws_, lawI)
    {
        const typename kinematicsViewsOf<Fields>::type views
        (
            fields, tp.lawIntegrationPointIDs_[lawI]
        );

        scales[lawI] =
            Fields::convergenceScale(laws_[lawI], views.kin, tp.states_[lawI]);
    }

    // One reduction per law, in law order, which every rank shares
    forAll(scales, lawI)
    {
        reduce(scales[lawI], maxOp<scalar>());
    }

    // Kept so that the boundary evaluations use the same scale as the
    // internal ones
    tp.lawConvergenceScales_ = scales;

    return scales;
}


const Foam::word&
Foam::mechanicalConstitutiveLawManager::topologyKeyFor
(
    const integrationPointTopology& topo
) const
{
    forAllIters(topologyCache_, iter)
    {
        if (&autoPtrRef(iter()) == &topo)
        {
            return iter.key();
        }
    }

    FatalErrorInFunction
        << "The integration-point topology of type " << topo.type()
        << " was not registered with this manager." << nl << nl
        << "    Its constitutive state, boundary state and restart entries "
        << "are all held against the key it was registered under, so a "
        << "topology the manager does not own has no state to find and "
        << "would silently be given a fresh one." << nl << nl
        << "    Obtain the topology from topologyFor(), registerTopology() "
        << "or compactCellTopologyFor() rather than constructing one and "
        << "passing it in."
        << exit(FatalError);

    return word::null;
}


Foam::mechanicalConstitutiveLawManager::topologyEntry&
Foam::mechanicalConstitutiveLawManager::topology
(
    const integrationPointTopology& topo
) const
{
    return topology(topologyKeyFor(topo));
}


Foam::mechanicalConstitutiveLawManager::topologyEntry&
Foam::mechanicalConstitutiveLawManager::topology
(
    const word& key
) const
{
    const integrationPointTopology& topo = autoPtrRef(topologyCache_[key]);

    // Return existing entry if already constructed
    if (topologyEntries_.found(key))
    {
        return autoPtrRef(topologyEntries_[key]);
    }

    // ---------------------------------------------------------------------
    // Construct new topology entry (lazy initialisation)
    // ---------------------------------------------------------------------

    DebugInfo
        << "Creating topologyEntry for " << key << endl;

    autoPtr<topologyEntry> entryPtr(new topologyEntry(topo));
    topologyEntries_.insert(key, entryPtr);
    topologyEntry& entry = autoPtrRef(topologyEntries_[key]);

    const label nLaws = laws_.size();

    entry.lawIntegrationPointIDs_.setSize(nLaws);
    entry.states_.setSize(nLaws);
    entry.boundaryStates_.setSize(nLaws);

    // Detect whether this topology supports boundary integration points
    entry.boundaryAware_ = topo.boundaryAware();

    // ---------------------------------------------------------------------
    // Build integration-point addressing per law
    // ---------------------------------------------------------------------

    forAll(laws_, lawI)
    {
        DynamicList<label> ipIDs;
        labelHashSet seen;   // only used when needed

        const labelList& cells = lawCells_[lawI];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const labelUList cellIPs = topo.cellIntegrationPointIDs(cellI);

            forAll(cellIPs, j)
            {
                const label ip = cellIPs[j];

                if (topo.requiresUniqueIntegrationPointsPerMaterial())
                {
                    if (seen.insert(ip))
                    {
                        ipIDs.append(ip);
                    }
                }
                else
                {
                    ipIDs.append(ip);
                }
            }
        }

        entry.lawIntegrationPointIDs_[lawI].transfer(ipIDs);
    }


    // ---------------------------------------------------------------------
    // Allocate and initialise constitutive states (per law)
    // ---------------------------------------------------------------------

    forAll(laws_, lawI)
    {
        entry.states_.set
        (
            lawI,
            new mechanicalConstitutiveLawState
            (
                entry.lawIntegrationPointIDs_[lawI].size()
            )
        );

        // Apply whatever the law declared and let it initialise its own
        // state, in that order. See applyStateSpec
        applyStateSpec
        (
            lawI,
            topo,
            entry.lawIntegrationPointIDs_[lawI],
            entry.states_[lawI]
        );
    }

    // ---------------------------------------------------------------------
    // Allocate and initialise boundary states (if applicable)
    // ---------------------------------------------------------------------

    if (entry.boundaryAware_)
    {
        forAll(laws_, lawI)
        {
            entry.boundaryStates_.set
            (
                lawI,
                new PtrList<mechanicalConstitutiveLawState>
                (
                    mesh_.boundary().size()
                )
            );

            forAll(mesh_.boundary(), patchI)
            {
                const label nFaces =
                    lawBoundaryFaces_[lawI][patchI].size();

                entry.boundaryStates_[lawI].set
                (
                    patchI,
                    new mechanicalConstitutiveLawState(nFaces)
                );

                applyStateSpecPatch
                (
                    lawI, patchI, entry.boundaryStates_[lawI][patchI]
                );
            }
        }
    }

    // Last, so that a restart overwrites the cold start defaults applied above
    // rather than the other way round
    setupStateRestart(entry, topo, key);

    return entry;
}


void Foam::mechanicalConstitutiveLawManager::checkTangentRequest
(
    const integrationPointTopology& topo,
    const tangentRequest req
)
{
    if (!needsFourthOrderTangent(req))
    {
        return;
    }

    if (!topo.supportsFourthOrderTangent())
    {
        FatalErrorInFunction
            << "A " << tangentRequestName(req) << " tangent was requested but "
            << "the integration-point topology " << topo.type()
            << " does not provide one." << nl
            << "Fourth-order tangents are only available where the integration "
            << "points are the locations at which a Jacobian operator "
            << "evaluates fluxes."
            << exit(FatalError);
    }

    if (topo.requiresUniqueIntegrationPointsPerMaterial() && laws_.size() > 1)
    {
        FatalErrorInFunction
            << "A " << tangentRequestName(req) << " tangent was requested on "
            << "topology " << topo.type() << " with " << laws_.size()
            << " mechanical constitutive laws." << nl
            << "Integration points of this topology are shared between cells, "
            << "so an integration point on a material interface belongs to "
            << "more than one law and has no single fourth-order tangent. "
            << "There is no meaningful collapse rule for a fourth-order "
            << "tangent: interface continuity is a normal-direction traction "
            << "and displacement matching problem, not an average."
            << exit(FatalError);
    }
}


Foam::labelList
Foam::mechanicalConstitutiveLawManager::currentMeshSizes() const
{
    const label nPatches = mesh_.boundary().size();

    labelList sizes(2*nPatches + 1);

    sizes[0] = mesh_.nCells();

    forAll(mesh_.boundary(), patchI)
    {
        const labelUList& faceCells = mesh_.boundary()[patchI].faceCells();

        sizes[patchI + 1] = faceCells.size();

        // Order-sensitive, so that faces moved between cells or reordered
        // within a patch change it; kept below the label maximum on every
        // label size
        long long checksum = 0;

        forAll(faceCells, i)
        {
            checksum = (31*checksum + faceCells[i] + 1) % 2147483629LL;
        }

        sizes[nPatches + patchI + 1] = label(checksum);
    }

    return sizes;
}


void Foam::mechanicalConstitutiveLawManager::calcLawBoundaryFaces()
{
    forAll(lawBoundaryFaces_, lawI)
    {
        lawBoundaryFaces_[lawI].clear();
        lawBoundaryFaces_[lawI].setSize(mesh_.boundary().size());

        forAll(lawBoundaryFaces_[lawI], patchI)
        {
            const labelUList& faceCells =
                mesh_.boundary()[patchI].faceCells();

            // Collected in ascending face order: the boundary constitutive
            // state is indexed by position in this list, so the order must be
            // reproducible
            DynamicList<label> curFaces(faceCells.size());

            forAll(faceCells, faceI)
            {
                if (cellToLaw_[faceCells[faceI]] == lawI)
                {
                    curFaces.append(faceI);
                }
            }

            lawBoundaryFaces_[lawI][patchI].transfer(curFaces);
        }
    }

    addressingMeshSizes_ = currentMeshSizes();
}


void Foam::mechanicalConstitutiveLawManager::updateAddressingIfTopologyChanged()
{
    const labelList sizes(currentMeshSizes());

    if (sizes == addressingMeshSizes_)
    {
        return;
    }

    if (sizes[0] != addressingMeshSizes_[0])
    {
        FatalErrorInFunction
            << "The number of cells changed from " << addressingMeshSizes_[0]
            << " to " << sizes[0] << "." << nl
            << "    The mechanicalConstitutiveLaw framework keeps its cell "
            << "addressing and cell states through a topology change, so it "
            << "supports changes that only move faces between patches, as "
            << "crackerFvMesh does, and not ones that add or remove cells."
            << exit(FatalError);
    }

    if (caseInputsPtr_.valid())
    {
        FatalErrorInFunction
            << "The mesh topology changed while one or more mechanical "
            << "constitutive law inputs are read from another case "
            << "directory." << nl
            << "    The input is copied by cell and face index from a static "
            << "source mesh, so case-directory inputs cannot be combined "
            << "with a topology-changing mesh."
            << exit(FatalError);
    }

    forAll(laws_, lawI)
    {
        if (declaresPersistentState(laws_[lawI]))
        {
            FatalErrorInFunction
                << "The mesh topology changed, and the mechanical "
                << "constitutive law " << lawNames_[lawI] << " carries "
                << "persistent state." << nl
                << "    The framework does not map a law's history onto new "
                << "boundary faces, so it cannot continue without silently "
                << "restarting that history on them."
                << exit(FatalError);
        }
    }

    DebugInfo
        << "Mesh topology changed: rebuilding the boundary addressing and "
        << "boundary states" << endl;

    calcLawBoundaryFaces();

    forAllIters(topologyEntries_, topoIter)
    {
        topologyEntry& entry = autoPtrRef(topoIter());

        // A topology whose points are the cells keeps its internal addressing
        // through a change that only moves faces between patches, and keeps a
        // state per patch face, sized from lawBoundaryFaces_. The others index
        // faces or points, which the topology itself would have to be rebuilt
        // for
        if (!entry.topology_.integrationPointsAreCells())
        {
            FatalErrorInFunction
                << "The mesh topology changed, and the integration-point "
                << "topology " << entry.topology_.type() << " is in use." << nl
                << "    Only a topology whose integration points are the "
                << "cells, such as "
                << cellCentredIntegrationPointTopology::typeName
                << ", is kept through a topology change."
                << exit(FatalError);
        }

        if (!entry.boundaryAware_)
        {
            continue;
        }

        // No law carries persistent state, so a cold boundary state is the
        // state these faces would have had anyway
        forAll(laws_, lawI)
        {
            PtrList<mechanicalConstitutiveLawState>& bStates =
                entry.boundaryStates_[lawI];

            bStates.clear();
            bStates.setSize(mesh_.boundary().size());

            forAll(mesh_.boundary(), patchI)
            {
                bStates.set
                (
                    patchI,
                    new mechanicalConstitutiveLawState
                    (
                        lawBoundaryFaces_[lawI][patchI].size()
                    )
                );

                applyStateSpecPatch(lawI, patchI, bStates[patchI]);
            }
        }
    }

    // Scratch and cached fields sized to the old mesh
    surfaceStressSumPtr_.clear();
    surfaceStressWeightPtr_.clear();
    surfaceTangentWeightPtr_.clear();
    pointStressSumPtr_.clear();
    pointStressWeightPtr_.clear();
    pointTangentWeightPtr_.clear();
    resetMaterialPropertyFields();
}


void Foam::mechanicalConstitutiveLawManager::updateOldTimeIfNeeded()
{
    // First, so that the states rolled over below are the current ones
    updateAddressingIfTopologyChanged();

    const label timeIndex = mesh_.time().timeIndex();

    if (timeIndex != curTimeIndex_)
    {
        if (debug)
        {
            InfoInFunction
                << "Updating old-time states for all topology entries"
                << endl;
        }

        // Loop over all topology entries
        forAllIters(topologyEntries_, topoIter)
        {
            topologyEntry& entry = autoPtrRef(topoIter());

            // Internal states
            forAll(entry.states_, lawI)
            {
                entry.states_[lawI].storeOldTime();
            }

            // Boundary states (if present)
            if (entry.boundaryAware_)
            {
                forAll(entry.boundaryStates_, lawI)
                {
                    forAll(entry.boundaryStates_[lawI], patchI)
                    {
                        entry.boundaryStates_[lawI][patchI].storeOldTime();
                    }
                }
            }
        }

        curTimeIndex_ = timeIndex;
    }
}


Foam::surfaceSymmTensorField&
Foam::mechanicalConstitutiveLawManager::surfaceStressSum() const
{
    if (!surfaceStressSumPtr_.valid())
    {
        surfaceStressSumPtr_.reset
        (
            new surfaceSymmTensorField
            (
                IOobject
                (
                    "surfaceStressSum",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );
    }

    return autoPtrRef(surfaceStressSumPtr_);
}


Foam::surfaceScalarField&
Foam::mechanicalConstitutiveLawManager::surfaceStressWeight() const
{
    if (!surfaceStressWeightPtr_.valid())
    {
        surfaceStressWeightPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "surfaceStressWeight",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    return autoPtrRef(surfaceStressWeightPtr_);
}


Foam::surfaceScalarField&
Foam::mechanicalConstitutiveLawManager::surfaceTangentWeight() const
{
    if (!surfaceTangentWeightPtr_.valid())
    {
        surfaceTangentWeightPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "surfaceTangentWeight",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimPressure, 0.0)
            )
        );
    }

    return autoPtrRef(surfaceTangentWeightPtr_);
}


Foam::pointSymmTensorField&
Foam::mechanicalConstitutiveLawManager::pointStressSum() const
{
    if (!pointStressSumPtr_.valid())
    {
        // Allocate the point mesh pointer, if needed
        if (!pMeshPtr_)
        {
            pMeshPtr_ = &pointMesh::New(mesh_);
        }

        pointStressSumPtr_.reset
        (
            new pointSymmTensorField
            (
                IOobject
                (
                    "pointStressSum",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *pMeshPtr_,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
            )
        );
    }

    return autoPtrRef(pointStressSumPtr_);
}


Foam::pointScalarField&
Foam::mechanicalConstitutiveLawManager::pointStressWeight() const
{
    if (!pointStressWeightPtr_.valid())
    {
        // Allocate the point mesh pointer, if needed
        if (!pMeshPtr_)
        {
            pMeshPtr_ = &pointMesh::New(mesh_);
        }

        pointStressWeightPtr_.reset
        (
            new pointScalarField
            (
                IOobject
                (
                    "pointStressWeight",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *pMeshPtr_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    return autoPtrRef(pointStressWeightPtr_);
}


Foam::pointScalarField&
Foam::mechanicalConstitutiveLawManager::pointTangentWeight() const
{
    if (!pointTangentWeightPtr_.valid())
    {
        // Allocate the point mesh pointer, if needed
        if (!pMeshPtr_)
        {
            pMeshPtr_ = &pointMesh::New(mesh_);
        }

        pointTangentWeightPtr_.reset
        (
            new pointScalarField
            (
                IOobject
                (
                    "pointTangentWeight",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *pMeshPtr_,
                dimensionedScalar("zero", dimPressure, 0.0)
            )
        );
    }

    return autoPtrRef(pointTangentWeightPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::mechanicalConstitutiveLawManager
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    pMeshPtr_(nullptr),
    curTimeIndex_(-1),
    laws_(),
    lawCells_(),
    lawBoundaryFaces_(),
    surfaceStressSumPtr_(),
    surfaceStressWeightPtr_(),
    surfaceTangentWeightPtr_(),
    pointStressSumPtr_(),
    pointStressWeightPtr_(),
    pointTangentWeightPtr_(),
    rhoPtr_(),
    kappaPtr_(),
    topologyCache_(),
    topologyEntries_(),
    compactFingerprints_(),
    addressingMeshSizes_(),
    caseInputsPtr_(),
    inputGatherer_
    (
        mesh_, laws_, lawCells_, caseInputsPtr_, reportedUnreadInputs_
    ),
    restartKinematicsAvailable_(true),
    stateSetup_
    (
        mesh_,
        laws_,
        lawNames_,
        lawCells_,
        lawBoundaryFaces_,
        stateWriteProxies_,
        restartKinematicsAvailable_
    )
{
    // Read the mechanical laws
    const PtrList<entry> lawEntries(dict.lookup("mechanical"));

    laws_.setSize(lawEntries.size());
    lawCells_.setSize(lawEntries.size());
    lawBoundaryFaces_.setSize(lawEntries.size());

    // Build the list of law names
    wordList lawNames(lawEntries.size());
    forAll(lawNames, lawI)
    {
        lawNames[lawI] = lawEntries[lawI].keyword();
    }

    // Entries only the removed legacy laws read. A case migrated from it
    // keeps them, and nothing reads them now, so say so rather than let a
    // setting that does nothing look as though it does
    forAll(lawEntries, lawI)
    {
        if (lawEntries[lawI].isDict())
        {
            reportRemovedLegacyEntries
            (
                lawEntries[lawI].dict(), lawEntries[lawI].keyword()
            );
        }
    }

    // Kept, because the state written for restart is named after the material
    lawNames_ = lawNames;

    // Create a map for each cell to its mechanical law
    cellToLaw_.setSize(mesh_.nCells(), -1);
    labelList& cellToLaw = cellToLaw_;

    // Plane stress is defined once for the whole mechanicalProperties
    // dictionary. It is injected into each law's sub-dictionary below rather
    // than looked up from the object registry, so that a mechanical
    // constitutive law depends on nothing but the dictionary it is given
    const Switch planeStress
    (
        dict.lookupOrDefault<Switch>("planeStress", false)
    );

    forAll(lawNames, lawI)
    {
        const word& lawName = lawNames[lawI];

        // Take a copy so that the shared settings can be injected
        dictionary lawDict(lawEntries[lawI].dict());

        if (lawDict.found("planeStress"))
        {
            FatalIOErrorInFunction(lawDict)
                << "'planeStress' is set once for all materials in the "
                << "mechanicalProperties dictionary and must not be given "
                << "inside the '" << lawName << "' sub-dictionary."
                << exit(FatalIOError);
        }

        lawDict.add("planeStress", planeStress);

        // Which directions the mesh actually solves in. A law that behaves
        // differently in two dimensions needs this and cannot ask the mesh
        // itself: a mechanical constitutive law is constructed from a
        // dictionary and nothing else. Injected for the same reason as
        // planeStress, and named as the mesh names it
        if (lawDict.found("solutionD"))
        {
            FatalIOErrorInFunction(lawDict)
                << "'solutionD' is supplied by the mesh and must not be given "
                << "inside the '" << lawName << "' sub-dictionary."
                << exit(FatalIOError);
        }

        lawDict.add("solutionD", mesh_.solutionD());

        // Construct law
        laws_.set
        (
            lawI,
            mechanicalConstitutiveLaw::New(lawDict)
        );

        // Any of the law's scalar inputs that the case reads from another
        // case directory rather than from this run. Read from the dictionary
        // as the user gave it, and for the inputs the whole law tree reads,
        // so that a sub-law's input can be sourced like its parent's
        {
            const HashTable<fileName> sources
            (
                mechanicalConstitutiveLawCaseInputs::readSources
                (
                    lawEntries[lawI].dict(),
                    lawName,
                    requiredScalarInputsRecursive(laws_[lawI])
                )
            );

            if (sources.size() && !caseInputsPtr_.valid())
            {
                caseInputsPtr_.reset
                (
                    new mechanicalConstitutiveLawCaseInputs
                    (
                        mesh_,
                        lawEntries.size()
                    )
                );
            }

            const wordList names(sources.sortedToc());

            forAll(names, i)
            {
                caseInputsPtr_->addSource(lawI, names[i], sources[names[i]]);
            }
        }

        if (lawNames.size() == 1)
        {
            // A single law covers the whole domain, so no cellZone is needed
            lawCells_[lawI] = Foam::identity(mesh_.nCells());
        }
        else // more than one material law
        {
            // Look up cell zone of the same name as the law
            const label zoneID = mesh_.cellZones().findZoneID(lawName);

            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "CellZone " << lawName << " not found"
                    << "When defining more than one mechanical constitutive "
                    << "law, each cell must belong to exactly one cellZone "
                    << "with the same name as the law."
                    << exit(FatalError);
            }

            lawCells_[lawI] = mesh_.cellZones()[zoneID];
        }

        // Check that each cell appears in only one cell zone
        forAll(lawCells_[lawI], i)
        {
            const label cellI = lawCells_[lawI][i];

            if (cellToLaw[cellI] != -1)
            {
                FatalErrorInFunction
                    << "Cell " << cellI
                    << " appears in more than one mechanical constitutive law: "
                    << lawNames[cellToLaw[cellI]] << " and " << lawName
                    << exit(FatalError);
            }

            cellToLaw[cellI] = lawI;
        }
    }

    // Check for any cells without a material law
    forAll(cellToLaw, cellI)
    {
        if (cellToLaw[cellI] == -1)
        {
            FatalErrorInFunction
                << "Cell " << cellI
                << " is not assigned to any mechanical constitutive law. "
                << "When defining more than one material, every cell must "
                << "belong to exactly one cellZone."
                << exit(FatalError);
        }
    }

    // An input read from a case directory by one material must be read that
    // way by every material that reads it. The sourced copy is registered and
    // written under the input's own name, so a material left to find the
    // input in this run would find the other material's copy instead, or
    // read it back from disk at the next write
    if (caseInputsPtr_.valid())
    {
        forAll(laws_, lawI)
        {
            const wordList names(requiredScalarInputsRecursive(laws_[lawI]));

            forAll(names, i)
            {
                if (caseInputsPtr_->found(lawI, names[i]))
                {
                    continue;
                }

                // A field already supplied by another model shadows every
                // case-directory source of the same name at evaluation time
                if (mesh_.foundObject<volScalarField>(names[i]))
                {
                    continue;
                }

                forAll(laws_, otherI)
                {
                    if (caseInputsPtr_->found(otherI, names[i]))
                    {
                        FatalErrorInFunction
                            << "Material " << lawEntries[otherI].keyword()
                            << " reads '" << names[i] << "' from a case "
                            << "directory, but material "
                            << lawEntries[lawI].keyword() << " reads it from "
                            << "this run." << nl
                            << "    Every material that reads an input must "
                            << "read it the same way: give "
                            << lawEntries[lawI].keyword() << " a case "
                            << "directory for '" << names[i] << "' too."
                            << exit(FatalError);
                    }
                }
            }
        }
    }

    // Set lawBoundaryFaces
    calcLawBoundaryFaces();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::~mechanicalConstitutiveLawManager()
{}


// * * * * * * * * * * * * State field output helpers  * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- One law's contribution to a state field written for viewing: its cells
//  and, patch by patch, its boundary faces, with the state holding each
struct stateFieldPart
{
    const mechanicalConstitutiveLawState* internal;
    const labelList* cells;
    List<const mechanicalConstitutiveLawState*> patchStates;
    List<const labelList*> patchFaces;
};

//- A state variable to be written, gathered over every law declaring it
struct stateFieldPlan
{
    word variable;
    word type;
    wordList path;
    List<stateFieldPart> parts;
};

//- Collect the persistent state variables of a law and of the laws below it
void collectStateFieldPlans
(
    const mechanicalConstitutiveLaw& law,
    const wordList& path,
    const stateFieldPart& part,
    HashTable<stateFieldPlan>& plans
)
{
    mechanicalConstitutiveLawStateSpec spec;
    law.declareState(spec);

    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        const mechanicalConstitutiveLawStateSpec::entry& e = es[i];

        // History only, as the restart writes: a prescribed field is already
        // on disk, as the user gave it, and a fixed one is a constant
        if (e.role != mechanicalConstitutiveLawStateSpec::stateRole::persistent)
        {
            continue;
        }

        word key(e.typeName);
        forAll(path, j)
        {
            key += '/' + path[j];
        }
        key += '/' + e.name;

        if (!plans.found(key))
        {
            stateFieldPlan plan;
            plan.variable = e.name;
            plan.type = e.typeName;
            plan.path = path;
            plans.insert(key, plan);
        }

        List<stateFieldPart>& parts = plans[key].parts;
        parts.setSize(parts.size() + 1);
        parts[parts.size() - 1] = part;
    }

    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        wordList childPath(path);
        childPath.setSize(path.size() + 1);
        childPath[path.size()] = childNames[i];

        stateFieldPart childPart(part);
        childPart.internal = &part.internal->child(childNames[i]);

        forAll(childPart.patchStates, patchI)
        {
            if (part.patchStates[patchI])
            {
                childPart.patchStates[patchI] =
                    &part.patchStates[patchI]->child(childNames[i]);
            }
        }

        collectStateFieldPlans
        (
            law.childLaw(childNames[i]), childPath, childPart, plans
        );
    }
}


//- Write one state variable, gathered from every law declaring it, as a
//  volField
template<class Type>
void writeStateVolField
(
    const fvMesh& mesh,
    const word& fieldName,
    const stateFieldPlan& plan
)
{
    GeometricField<Type, fvPatchField, volMesh> vf
    (
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero),
        calculatedFvPatchField<Type>::typeName
    );

    Field<Type>& vfI = primitiveFieldRef(vf);

    forAll(plan.parts, partI)
    {
        const stateFieldPart& part = plan.parts[partI];

        const Field<Type>& f =
            stateFieldAccess<Type>::get(*part.internal, plan.variable);
        const labelList& cells = *part.cells;

        forAll(cells, i)
        {
            vfI[cells[i]] = f[i];
        }

        forAll(part.patchStates, patchI)
        {
            if (!part.patchStates[patchI])
            {
                continue;
            }

            const Field<Type>& pf =
                stateFieldAccess<Type>::get
                (
                    *part.patchStates[patchI], plan.variable
                );
            const labelList& faces = *part.patchFaces[patchI];

            if (pf.size() != faces.size())
            {
                continue;
            }

            fvPatchField<Type>& pvf = boundaryFieldRef(vf)[patchI];

            forAll(faces, i)
            {
                pvf[faces[i]] = pf[i];
            }
        }
    }

    vf.write();
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const Foam::volScalarField& Foam::mechanicalConstitutiveLawManager::rho() const
{
    if (!rhoPtr_.valid())
    {
        // Not registered. This is the manager's private cache; the solid
        // model's own copy (solidModel::makeRho) is the field registered as
        // "rho". Were this one to hold the name, the
        // solid model's copy could not register, and a topology-changing mesh
        // (crackerFvMesh) maps only registered fields, so the solid model's
        // rho would keep its old boundary sizes after the mesh changed
        rhoPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rhoFromLaws",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false  // Do not register
                ),
                mesh_,
                dimensionedScalar("rho", dimDensity, 0.0),
                "zeroGradient"
            )
        );

        volScalarField& rhoField = autoPtrRef(rhoPtr_);

        forAll(laws_, lawI)
        {
            const dimensionedScalar rhoLaw = laws_[lawI].rho();
            const labelList& cells = lawCells_[lawI];

            forAll(cells, i)
            {
                rhoField[cells[i]] = rhoLaw.value();
            }
        }

        rhoField.correctBoundaryConditions();
    }

    return rhoPtr_();
}


const Foam::mechanicalConstitutiveLaw&
Foam::mechanicalConstitutiveLawManager::singleLaw() const
{
    if (laws_.size() != 1)
    {
        FatalErrorInFunction
            << "singleLaw() was asked for, but this manager holds "
            << laws_.size() << " materials." << nl
            << "    The caller is a solid model whose formulation has one "
            << "set of material constants, so a second material has no "
            << "meaning for it."
            << exit(FatalError);
    }

    return laws_[0];
}


const Foam::volScalarField&
Foam::mechanicalConstitutiveLawManager::kappa() const
{
    if (!kappaPtr_.valid())
    {
        kappaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "kappa",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                // A bulk modulus, not a density. This was dimDensity,
                // copied from rho() directly above, and every consumer then
                // used it as a pressure
                dimensionedScalar("kappa", dimPressure, 0.0),
                "zeroGradient"
            )
        );

        volScalarField& kappaField = autoPtrRef(kappaPtr_);

        forAll(laws_, lawI)
        {
            const dimensionedScalar kappaLaw = laws_[lawI].kappa();
            const labelList& cells = lawCells_[lawI];

            forAll(cells, i)
            {
                kappaField[cells[i]] = kappaLaw.value();
            }
        }

        kappaField.correctBoundaryConditions();
    }

    return kappaPtr_();
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::registerTopology
(
    const word& key,
    autoPtr<integrationPointTopology> topoPtr
) const
{
    if (!topoPtr.valid())
    {
        FatalErrorInFunction
            << "Null integrationPointTopology registered under key " << key
            << exit(FatalError);
    }

    // Already registered?
    if (topologyCache_.found(key))
    {
        const integrationPointTopology& existing =
            autoPtrRef(topologyCache_[key]);

        if (existing.type() != topoPtr->type())
        {
            FatalErrorInFunction
                << "An integrationPointTopology of type " << existing.type()
                << " is already registered under the key " << key << ", but a "
                << "topology of type " << topoPtr->type() << " was supplied."
                << nl
                << "Keys must identify a topology uniquely: the constitutive "
                << "state is keyed on the topology object, so two topologies "
                << "sharing a key would share history variables."
                << exit(FatalError);
        }

        // Discard the supplied topology and keep the one already in use, so
        // that its constitutive state is preserved
        return topology(existing).topology_;
    }

    topologyCache_.insert(key, topoPtr);

    return topology(autoPtrRef(topologyCache_[key])).topology_;
}


void Foam::mechanicalConstitutiveLawManager::resetMaterialPropertyFields()
{
    rhoPtr_.clear();
    kappaPtr_.clear();
}


bool Foam::mechanicalConstitutiveLawManager::resolveBoundaryEvaluation
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
) const
{
    const labelList patchIPs
    (
        topo.boundaryIntegrationPointIDs(patchI)
    );

    if (patchIPs.empty())
    {
        return false;
    }

    // An empty patch has no fvPatch faces, so lawBoundaryFaces_
    // is empty for it, yet its polyPatch faces still occupy slots
    // in this topology's index space. Skipping it would leave
    // those slots unwritten, which is what a caller sizing its
    // list to nIntegrationPoints() then reads. So address them
    // through the polyPatch instead, and take the law from the
    // owner cell.
    //
    // These faces take no part in the finite volume
    // discretisation - that is what empty means - so they are
    // evaluated against a shadow of the law's state. They get a
    // well-defined value without committing history for a face
    // the discretisation does not have
    const bool viaPolyPatch =
        isA<emptyFvPatch>(mesh_.boundary()[patchI]);


    if (viaPolyPatch)
    {
        const labelUList& ownCells =
            mesh_.boundaryMesh()[patchI].faceCells();

        DynamicList<label> ids(ownCells.size());
        forAll(ownCells, faceI)
        {
            if (cellToLaw_[ownCells[faceI]] == lawI)
            {
                ids.append(patchIPs[faceI]);
            }
        }

        ipIDs.transfer(ids);
    }
    else
    {
        const labelList& faces = lawBoundaryFaces_[lawI][patchI];

        // This law's faces on this patch, as integration-point
        // indices
        ipIDs.setSize(faces.size());
        forAll(faces, i)
        {
            ipIDs[i] = patchIPs[faces[i]];
        }
    }

    if (ipIDs.empty())
    {
        return false;
    }

    // sized to the fvPatch, which has no faces there. Its faces
    // take no part in the discretisation, so give them a scratch
    // state of the right size rather than a shadow of a state
    // that is the wrong length. Nothing reads it afterwards

    if (viaPolyPatch)
    {
        bScratchPtr.set
        (
            new mechanicalConstitutiveLawState(ipIDs.size())
        );

        // The law must be given the chance to register the state
        // fields it reads, exactly as for a real state. Without
        // this a history-dependent law fails looking up its own
        // history, e.g. epsilonP for the plastic law, and a
        // declared field would be created empty on first access.
        //
        // A prescribed field keeps its declared default here. The
        // faces are addressed through the polyPatch because the
        // fvPatch has none, so there is no patch field to read
        // them from; they take no part in the discretisation, so
        // nothing downstream depends on the difference
        applyStateSpecScratch(lawI, bScratchPtr());
    }
    else if (preserveState)
    {
        bShadowPtr.set
        (
            new mechanicalConstitutiveLawState
            (
                tp.boundaryStates_[lawI][patchI],
                mechanicalConstitutiveLawState::SHADOW
            )
        );
    }

    bState =
        viaPolyPatch
      ? &bScratchPtr()
      : preserveState
      ? &bShadowPtr()
      : &tp.boundaryStates_[lawI][patchI];

    return true;
}


template<class Fields>
void Foam::mechanicalConstitutiveLawManager::evaluateFlat
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
)
{
    // Every processor refreshes every coupling input source here, before any
    // material is skipped for having no points on it, since reading is
    // collective
    refreshScalarInputs();

    const word context = Fields::flatContext();
    const label nIP = topo.nIntegrationPoints();

    checkKinematicsListSizes(nIP, fields, context);
    checkIntegrationPointListSize(nIP, stress.size(), "stress", context);

    checkTangentRequest(topo, tangentReq);

    checkTangentStorage
    (
        scalarTangentPtr != nullptr,
        fourthOrderTangentPtr != nullptr,
        tangentReq,
        context
    );

    // A volumetric split asked for here is checked here. The GeometricField
    // overloads validate it where the request is made; these take the storage
    // directly and validated neither that a law can fill it nor that it is the
    // right length, so an unsupported law returned a total stress the caller
    // would read as isochoric, and a short list was indexed through the
    // topology's addressing
    if (volumetricPtr)
    {
        checkIntegrationPointListSize
        (
            nIP, volumetricPtr->size(), "volumetricResponse", context
        );

        checkVolumetricSplitSupported(context);
    }

    if (scalarTangentPtr)
    {
        checkIntegrationPointListSize
        (
            nIP, scalarTangentPtr->size(), "scalarTangent", context
        );
    }

    if (fourthOrderTangentPtr)
    {
        checkIntegrationPointListSize
        (
            nIP, fourthOrderTangentPtr->size(), "fourthOrderTangent", context
        );
    }

    if (topo.requiresUniqueIntegrationPointsPerMaterial() && laws_.size() > 1)
    {
        FatalErrorInFunction
            << "The flat-list update does not perform stress collapse, but "
            << "topology " << topo.type() << " shares integration points "
            << "between cells, and there are " << laws_.size()
            << " mechanical constitutive laws, so an integration point on a "
            << "material interface would be written more than once."
            << Fields::multiLawHint()
            << exit(FatalError);
    }

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    topologyEntry& tp = topology(topo);

    // One collective per law, before any of them is evaluated.
    //
    // A law that normalises its convergence test by a scale over its points
    // needs that scale to be the same everywhere, and cannot reduce for
    // itself: the loop below skips a law where this rank holds none of its
    // points, so the reductions would not pair up. Asked of every law on
    // every rank, including where it has no points, the count matches
    const scalarList lawScales(convergenceScales(tp, fields));

    // A caller that wants a tangent independent of history gets a state
    // prepared exactly as a fresh run prepares one: declared defaults, the
    // law's own initialisation, and any prescribed field. That is what the
    // material looked like before it was loaded, so the tangent is the same
    // whether the run started here or was continued.
    //
    // Prepared for every law on every rank, before the loop below skips the
    // laws this rank holds no points of: a prescribed field is read from a
    // file, and reading one can be collective
    PtrList<mechanicalConstitutiveLawState> coldStates;
    if (coldState)
    {
        coldStates.setSize(laws_.size());

        forAll(laws_, lawI)
        {
            coldStates.set
            (
                lawI,
                new mechanicalConstitutiveLawState(tp.states_[lawI].size())
            );

            applyStateSpec
            (
                lawI, topo, tp.lawIntegrationPointIDs_[lawI], coldStates[lawI]
            );
        }
    }

    // Loop over mechanical constitutive laws
    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        if (ipIDs.empty())
        {
            continue;
        }

        // Live inputs for this law's evaluation. Built per law, because a
        // coupling input is handed over as a view of that law's own
        // integration points, and built here so that it is passed through
        // every evaluation below, including each finite-difference
        // perturbation, with no per-call forwarding to get wrong
        const mechanicalConstitutiveLawInputs inputs
        (
            lawInputs(lawI, topo, ipIDs, dt, tp)
        );

        inputs.setConvergenceScale(lawScales[lawI]);

        // A tangent query evaluates against a shadow of the law's state: the
        // shadow aliases the old-time fields, so history is read but never
        // written, and the law's outputs land where they are discarded
        autoPtr<mechanicalConstitutiveLawState> shadowPtr;
        if (preserveState && !coldState)
        {
            shadowPtr.set
            (
                new mechanicalConstitutiveLawState
                (
                    tp.states_[lawI],
                    mechanicalConstitutiveLawState::SHADOW
                )
            );
        }

        mechanicalConstitutiveLawState& lawState =
            coldState
          ? coldStates[lawI]
          : (preserveState ? shadowPtr() : tp.states_[lawI]);

        // Views into integration-point data (no copies), and the kinematics
        // wrapper built on them
        const typename kinematicsViewsOf<Fields>::type views(fields, ipIDs);
        UIndirectList<symmTensor> stressView(stress, ipIDs);

        // Constitutive response
        evaluateResponse
        (
            laws_[lawI],
            views.kin,
            inputs,
            lawState,
            stressView,
            ipIDs,
            scalarTangentPtr,
            fourthOrderTangentPtr,
            tangentReq,
            volumetricPtr,
            !preserveState
        );
    }

    // Boundary integration points.
    // The topology's cell-to-integration-point map covers internal points
    // only, so without this every boundary entry of the caller's storage is
    // left exactly as it was found - which, for a caller that sized its list
    // to nIntegrationPoints(), means unwritten memory.
    // A topology with no boundary slots in its flat index space returns an
    // empty list per patch below and nothing happens, which is the right
    // outcome for a cell-centred topology: it is boundaryAware because it
    // keeps a state per patch, not because its index space extends past the
    // cells
    if (tp.boundaryAware_)
    {
        forAll(laws_, lawI)
        {
            forAll(mesh_.boundary(), patchI)
            {
                labelList ipIDs;
                autoPtr<mechanicalConstitutiveLawState> bScratchPtr;
                autoPtr<mechanicalConstitutiveLawState> bShadowPtr;
                mechanicalConstitutiveLawState* bStatePtr = nullptr;

                if
                (
                   !resolveBoundaryEvaluation
                    (
                        topo,
                        tp,
                        lawI,
                        patchI,
                        preserveState,
                        ipIDs,
                        bScratchPtr,
                        bShadowPtr,
                        bStatePtr
                    )
                )
                {
                    continue;
                }

                mechanicalConstitutiveLawState& bState = *bStatePtr;

                // Live inputs for this law on this patch. The boundary points
                // are a different set from the internal ones, so the coupling
                // input has to be gathered for them rather than reused, and
                // they are judged by the same scale as the law's internal
                // points
                const mechanicalConstitutiveLawInputs inputs
                (
                    lawInputs(lawI, topo, ipIDs, dt, tp)
                );

                inputs.setConvergenceScale(lawScales[lawI]);

                const typename kinematicsViewsOf<Fields>::type views
                (
                    fields, ipIDs
                );
                UIndirectList<symmTensor> stressView(stress, ipIDs);

                evaluateResponse
                (
                    laws_[lawI],
                    views.kin,
                    inputs,
                    bState,
                    stressView,
                    ipIDs,
                    scalarTangentPtr,
                    fourthOrderTangentPtr,
                    tangentReq,
                    volumetricPtr,
                    !preserveState
                );
            }
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::evaluateSmallStrain
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
    const bool coldState,
    UList<scalar>* volumetricPtr
)
{
    evaluateFlat
    (
        topo,
        smallStrainKinematicsFields(gradD, gradD0),
        dt,
        stress,
        scalarTangentPtr,
        fourthOrderTangentPtr,
        tangentReq,
        preserveState,
        coldState,
        volumetricPtr
    );
}


void Foam::mechanicalConstitutiveLawManager::evaluateFiniteStrain
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
    UList<scalar>* volumetricPtr
)
{
    // No cold state on the finite-strain path: its only caller, the
    // small-strain implicit stiffness, does not come this way
    evaluateFlat
    (
        topo,
        finiteStrainKinematicsFields(F, F0, Finv, Finv0, J, J0),
        dt,
        stress,
        scalarTangentPtr,
        fourthOrderTangentPtr,
        tangentReq,
        preserveState,
        false,
        volumetricPtr
    );
}


void Foam::mechanicalConstitutiveLawManager::checkKinematicsListSizes
(
    const label nIP,
    const smallStrainKinematicsFields& fields,
    const word& context
) const
{
    checkIntegrationPointListSize(nIP, fields.gradD.size(), "gradD", context);
    checkIntegrationPointListSize(nIP, fields.gradD0.size(), "gradD0", context);
}


void Foam::mechanicalConstitutiveLawManager::checkKinematicsListSizes
(
    const label nIP,
    const finiteStrainKinematicsFields& fields,
    const word& context
) const
{
    checkIntegrationPointListSize(nIP, fields.F.size(), "F", context);
    checkIntegrationPointListSize(nIP, fields.F0.size(), "F0", context);
    checkIntegrationPointListSize(nIP, fields.Finv.size(), "Finv", context);
    checkIntegrationPointListSize(nIP, fields.Finv0.size(), "Finv0", context);
    checkIntegrationPointListSize(nIP, fields.J.size(), "J", context);
    checkIntegrationPointListSize(nIP, fields.J0.size(), "J0", context);
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const integrationPointTopology& topo,
    const UList<tensor>& gradD,
    const UList<tensor>& gradD0,
    const scalar dt,
    UList<symmTensor>& stress,
    UList<scalar>* scalarTangentPtr,
    UList<mat66>* fourthOrderTangentPtr,
    const tangentRequest tangentReq
)
{
    evaluateSmallStrain
    (
        topo,
        gradD,
        gradD0,
        dt,
        stress,
        scalarTangentPtr,
        fourthOrderTangentPtr,
        tangentReq,
        false           // commit the constitutive state
    );
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
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
    const tangentRequest tangentReq
)
{
    evaluateFiniteStrain
    (
        topo,
        F,
        F0,
        Finv,
        Finv0,
        J,
        J0,
        dt,
        stress,
        scalarTangentPtr,
        fourthOrderTangentPtr,
        tangentReq,
        false           // commit the constitutive state
    );
}


void Foam::mechanicalConstitutiveLawManager::updateTangentSmallStrain
(
    const integrationPointTopology& topo,
    const UList<tensor>& gradD,
    const UList<tensor>& gradD0,
    const scalar dt,
    UList<scalar>* scalarTangentPtr,
    UList<mat66>* fourthOrderTangentPtr,
    const tangentRequest tangentReq,
    const bool coldState
)
{
    if (tangentReq == tangentRequest::none)
    {
        FatalErrorInFunction
            << "updateTangentSmallStrain was called with tangentRequest::none, "
            << "so there is nothing to compute."
            << exit(FatalError);
    }

    // A constitutive law produces a stress alongside its tangent, so give it
    // somewhere to put one that is not the caller's storage
    evaluateSmallStrain
    (
        topo,
        gradD,
        gradD0,
        dt,
        scratchStress(topo.nIntegrationPoints()),
        scalarTangentPtr,
        fourthOrderTangentPtr,
        tangentReq,
        true,           // preserve the constitutive state
        coldState
    );
}


void Foam::mechanicalConstitutiveLawManager::updateTangentFiniteStrain
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
)
{
    if (tangentReq == tangentRequest::none)
    {
        FatalErrorInFunction
            << "updateTangentFiniteStrain was called with "
            << "tangentRequest::none, so there is nothing to compute."
            << exit(FatalError);
    }

    // A constitutive law produces a stress alongside its tangent, so give it
    // somewhere to put one that is not the caller's storage
    evaluateFiniteStrain
    (
        topo,
        F,
        F0,
        Finv,
        Finv0,
        J,
        J0,
        dt,
        scratchStress(topo.nIntegrationPoints()),
        scalarTangentPtr,
        fourthOrderTangentPtr,
        tangentReq,
        true            // preserve the constitutive state
    );
}


void Foam::mechanicalConstitutiveLawManager::updateScalarTangent
(
    const volTensorField& gradD,
    const volTensorField& gradD0,
    const scalar dt,
    volScalarField& scalarTangent,
    const tangentRequest tangentReq,
    const bool coldState
)
{
    checkMeshConsistency(mesh_, gradD.mesh(), gradD.name());
    checkMeshConsistency(mesh_, gradD0.mesh(), gradD0.name());
    checkMeshConsistency(mesh_, scalarTangent.mesh(), scalarTangent.name());

    if (!needsScalarTangent(tangentReq))
    {
        FatalErrorInFunction
            << "updateScalarTangent was asked for a "
            << tangentRequestName(tangentReq) << " tangent." << nl
            << "This interface returns a scalar tangent at cell centres, so "
            << "the request must be scalar or scalarDeviatoric."
            << exit(FatalError);
    }

    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    scalarField& tangent = Foam::primitiveFieldRef(scalarTangent);

    updateTangentSmallStrain
    (
        topo,
        Foam::primitiveField(gradD),
        Foam::primitiveField(gradD0),
        dt,
        &tangent,
        nullptr,
        tangentReq,
        coldState
    );

    // The flat-list primitive fills internal integration points only, so give
    // the boundary usable values. A scalar tangent is a per-cell material
    // property, and a boundary face belongs to the material of its owner cell,
    // so taking the patch-internal value is exact rather than an
    // approximation. Without this the boundary stays at whatever the field was
    // constructed with, which a caller forming 1/tangent then divides by
    forAll(scalarTangent.boundaryField(), patchI)
    {
        if (!scalarTangent.boundaryField()[patchI].coupled())
        {
            Foam::boundaryFieldRef(scalarTangent)[patchI] =
                scalarTangent.boundaryField()[patchI].patchInternalField();
        }
    }

    // Sync the coupled patches
    scalarTangent.correctBoundaryConditions();
}


void Foam::mechanicalConstitutiveLawManager::updateScalarTangentFiniteStrain
(
    const volTensorField& F,
    const volTensorField& F0,
    const volTensorField& Finv,
    const volTensorField& Finv0,
    const volScalarField& J,
    const volScalarField& J0,
    const scalar dt,
    volScalarField& scalarTangent,
    const tangentRequest tangentReq
)
{
    checkMeshConsistency(mesh_, F.mesh(), F.name());
    checkMeshConsistency(mesh_, F0.mesh(), F0.name());
    checkMeshConsistency(mesh_, Finv.mesh(), Finv.name());
    checkMeshConsistency(mesh_, Finv0.mesh(), Finv0.name());
    checkMeshConsistency(mesh_, J.mesh(), J.name());
    checkMeshConsistency(mesh_, J0.mesh(), J0.name());
    checkMeshConsistency(mesh_, scalarTangent.mesh(), scalarTangent.name());

    if (!needsScalarTangent(tangentReq))
    {
        FatalErrorInFunction
            << "updateScalarTangentFiniteStrain was asked for a "
            << tangentRequestName(tangentReq) << " tangent." << nl
            << "This interface returns a scalar tangent at cell centres, so "
            << "the request must be scalar or scalarDeviatoric."
            << exit(FatalError);
    }

    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    scalarField& tangent = Foam::primitiveFieldRef(scalarTangent);

    updateTangentFiniteStrain
    (
        topo,
        Foam::primitiveField(F),
        Foam::primitiveField(F0),
        Foam::primitiveField(Finv),
        Foam::primitiveField(Finv0),
        Foam::primitiveField(J),
        Foam::primitiveField(J0),
        dt,
        &tangent,
        nullptr,
        tangentReq
    );

    // As in updateScalarTangent: a boundary face belongs to the material of
    // its owner cell, so the patch-internal value is exact, and a caller
    // forming 1/tangent must not be handed whatever the field was constructed
    // with
    forAll(scalarTangent.boundaryField(), patchI)
    {
        if (!scalarTangent.boundaryField()[patchI].coupled())
        {
            Foam::boundaryFieldRef(scalarTangent)[patchI] =
                scalarTangent.boundaryField()[patchI].patchInternalField();
        }
    }

    // Sync the coupled patches
    scalarTangent.correctBoundaryConditions();
}


template<class Fields, class PatchFieldsFn>
void Foam::mechanicalConstitutiveLawManager::updateStressVolBoundary
(
    topologyEntry& tp,
    const volTensorField& lead,
    const PatchFieldsFn& patchFields,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq,
    volScalarField* volumetricResponsePtr
)
{
    // Boundary constitutive response uses independent state objects, allowing
    // history-dependent laws to operate correctly on boundary faces
    if (!tp.boundaryAware_)
    {
        return;
    }

    forAll(laws_, lawI)
    {
        forAll(lead.boundaryField(), patchI)
        {
            // A coupled patch is not evaluated: it takes its values from the
            // neighbouring processor when the boundary conditions are
            // corrected
            if (lead.boundaryField()[patchI].coupled())
            {
                continue;
            }

            // Select all faces on the patch whose adjacent cell is in this
            // material
            const labelList& faces = lawBoundaryFaces_[lawI][patchI];

            if (faces.empty() || isA<emptyFvPatch>(mesh_.boundary()[patchI]))
            {
                continue;
            }

            // Live inputs for this law on this patch. The patch values are what
            // a boundary face sees, not the values in the cells behind it
            const mechanicalConstitutiveLawInputs inputs
            (
                lawInputsPatch(lawI, patchI, faces, dt, tp)
            );

            // The same scale the internal points were evaluated with. Taking
            // it over this rank's faces instead would make the convergence
            // tolerance depend on where the mesh was cut
            if (lawI < tp.lawConvergenceScales_.size())
            {
                inputs.setConvergenceScale(tp.lawConvergenceScales_[lawI]);
            }

            // Views into the kinematic and stress fields for this material,
            // which do not copy data, and the kinematics built on them
            const typename kinematicsViewsOf<Fields>::type views
            (
                patchFields(patchI), faces
            );
            UIndirectList<symmTensor> stressView
            (
                Foam::boundaryFieldRef(stress)[patchI], faces
            );

            // This path computes no fourth-order tangent, so none is offered
            evaluateResponse
            (
                laws_[lawI],
                views.kin,
                inputs,
                tp.boundaryStates_[lawI][patchI],
                stressView,
                faces,
                scalarTangentPtr
              ? &scalarTangentPtr->boundaryField()[patchI]
              : nullptr,
                static_cast<const UList<mat66>*>(nullptr),
                tangentReq,
                volumetricResponsePtr
              ? &volumetricResponsePtr->boundaryField()[patchI]
              : nullptr
            );
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const volTensorField& gradD,
    const volTensorField& gradD0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq,
    volScalarField* volumetricResponsePtr
)
{
    // Check gradD is defined on the correct mesh
    checkMeshConsistency(mesh_, gradD.mesh(), gradD.name());
    checkMeshConsistency(mesh_, gradD0.mesh(), gradD0.name());
    checkMeshConsistency(mesh_, stress.mesh(), stress.name());
    if (scalarTangentPtr)
    {
        checkMeshConsistency
        (
            mesh_, scalarTangentPtr->mesh(), scalarTangentPtr->name()
        );
    }

    if (volumetricResponsePtr)
    {
        checkMeshConsistency
        (
            mesh_,
            volumetricResponsePtr->mesh(),
            volumetricResponsePtr->name()
        );

        // Checked where the request is made, as on the finite-strain path
        checkVolumetricSplitSupported("updateStressSmallStrain");
    }

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Look up the map and state for cell-based topologies
    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    // Update the internal field via the flat-list primitive: a cell-centred
    // topology has one integration point per cell, so the internal fields are
    // already in the flat form it expects
    evaluateSmallStrain
    (
        topo,
        Foam::primitiveField(gradD),
        Foam::primitiveField(gradD0),
        dt,
        Foam::primitiveFieldRef(stress),
        scalarTangentPtr
      ? &Foam::primitiveFieldRef(*scalarTangentPtr)
      : nullptr,
        nullptr,
        tangentReq,
        false,          // commit the constitutive state
        false,          // not a cold-state query
        volumetricResponsePtr
      ? &Foam::primitiveFieldRef(*volumetricResponsePtr)
      : nullptr
    );

    topologyEntry& tp = topology(topo);

    // The boundary faces, each with its own state
    updateStressVolBoundary<smallStrainKinematicsFields>
    (
        tp,
        gradD,
        [&](const label patchI)
        {
            return smallStrainKinematicsFields
            (
                gradD.boundaryField()[patchI],
                gradD0.boundaryField()[patchI]
            );
        },
        dt,
        stress,
        scalarTangentPtr,
        tangentReq,
        volumetricResponsePtr
    );

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
    }

    // The volumetric response is corrected with the others. Coupled patches
    // are skipped when it is evaluated, exactly as the stress is, so without
    // this its processor values are whatever the field was constructed with.
    // Only the internal field is read today, which makes that harmless rather
    // than correct - and harmless-for-now is a poor thing to leave for
    // whoever next reads a boundary value
    if (volumetricResponsePtr)
    {
        volumetricResponsePtr->correctBoundaryConditions();
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const surfaceTensorField& gradD,
    const surfaceTensorField& gradD0,
    const scalar dt,
    surfaceSymmTensorField& stress,
    const stressCollapseRule collapseRule,
    surfaceScalarField* scalarTangentPtr,
    List<mat66>* fourthOrderTangentPtr,
    const tangentRequest tangentReq
)
{
    // In parallel a material interface can lie on a processor boundary,
    // where each side holds only its own material's contribution and nothing
    // exchanges them before the collapse, so the two sides would disagree
    // with each other and with serial. Refused rather than answered wrongly,
    // as the finite-strain face overload does
    if (Pstream::parRun() && laws_.size() > 1)
    {
        FatalErrorInFunction
            << "The face small-strain stress update does not support more "
            << "than one material in parallel: contributions at a material "
            << "interface on a processor boundary are not reconciled"
            << exit(FatalError);
    }

    // Check gradD is defined on the correct mesh
    checkMeshConsistency(mesh_, gradD.mesh(), gradD.name());
    checkMeshConsistency(mesh_, gradD0.mesh(), gradD0.name());
    checkMeshConsistency(mesh_, stress.mesh(), stress.name());
    if (scalarTangentPtr)
    {
        checkMeshConsistency
        (
            mesh_, scalarTangentPtr->mesh(), scalarTangentPtr->name()
        );
    }

    // Both storages are checked, whichever of them is given: a scalar
    // tangent does not excuse a fourth-order one of the wrong size
    if (fourthOrderTangentPtr)
    {
        if (fourthOrderTangentPtr->size() != gradD.mesh().nInternalFaces())
        {
            FatalErrorInFunction
                << "Inconsistent fourthOrderTangent size" << nl
                << "Expected size = " << gradD.mesh().nInternalFaces() << nl
                << "Got: " << fourthOrderTangentPtr->size()
                << exit(FatalError);
        }
    }

    checkTangentStorage
    (
        scalarTangentPtr != nullptr,
        fourthOrderTangentPtr != nullptr,
        tangentReq,
        "updateStressSmallStrain (surfaceField)"
    );

    // Look up the map and state for face-based topologies
    const integrationPointTopology& topo =
        topologyFor(faceCentredIntegrationPointTopology::typeName);

    topologyEntry& tp = topology(topo);

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Every processor refreshes every coupling input source here, before any
    // material is skipped for having no points on it, since reading is
    // collective
    refreshScalarInputs();

    // The scale each law's convergence test is normalised by, taken over its
    // internal faces on every rank, and used on its boundary faces too, as
    // the other paths do. Every rank calls this, whatever points it holds
    const scalarList scales
    (
        convergenceScales
        (
            tp,
            smallStrainKinematicsFields
            (
                gradD.internalField(), gradD0.internalField()
            )
        )
    );

    surfaceSymmTensorField& stressSum = surfaceStressSum();
    surfaceScalarField& weightSum = surfaceStressWeight();

    stressSum = dimensionedSymmTensor("0", dimPressure, symmTensor::zero);
    weightSum = 0.0;

    surfaceScalarField* tangentWeightPtr = nullptr;

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        surfaceScalarField& tangentWeight = surfaceTangentWeight();
        tangentWeight = dimensionedScalar("0", dimPressure, 0.0);
        tangentWeightPtr = &tangentWeight;
    }

    checkTangentRequest(topo, tangentReq);

    // Loop over constitutive laws
    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Live inputs for this law, as a view of its own integration points
        const mechanicalConstitutiveLawInputs inputs
        (
            lawInputs(lawI, topo, ipIDs, dt, tp)
        );

        inputs.setConvergenceScale(scales[lawI]);

        const smallStrainKinematicsViews views
        (
            smallStrainKinematicsFields
            (
                gradD.internalField(), gradD0.internalField()
            ),
            ipIDs
        );

        UIndirectList<symmTensor> stressView
        (
            stress.internalField(), ipIDs
        );

        evaluateResponse
        (
            laws_[lawI],
            views.kin,
            inputs,
            tp.states_[lawI],
            stressView,
            ipIDs,
            scalarTangentPtr,
            fourthOrderTangentPtr,
            tangentReq
        );

        // Update stress accumulation fields used for stress collapse
        forAll(ipIDs, i)
        {
            const label faceI = ipIDs[i];

            stressSum[faceI] += stress[faceI];
            weightSum[faceI] += 1.0;

            if (scalarTangentPtr && needsScalarTangent(tangentReq))
            {
                const scalar K = (*scalarTangentPtr)[faceI];

                if (collapseRule == stressCollapseRule::harmonic)
                {
                    (*tangentWeightPtr)[faceI] += 1.0/max(K, SMALL);
                }
                else
                {
                    // 'average', and 'none', which is only reached with a
                    // single contribution and so is the same sum
                    (*tangentWeightPtr)[faceI] += K;
                }
            }
        }

        // Optionally, update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary
        // faces
        if (tp.boundaryAware_)
        {
            forAll(gradD.boundaryField(), patchI)
            {
                // Coupled patches are included. A processor face carries a
                // real stress that only this rank can compute: unlike a
                // volField, a surface field's coupled patch is not filled in
                // by correctBoundaryConditions(), so skipping it would leave
                // the stress there at zero

                // Select all faces on the patch for which the adjacent
                // cell is in this material
                const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                if
                (
                    faces.empty()
                 || isA<emptyFvPatch>(mesh_.boundary()[patchI])
                )
                {
                    continue;
                }

                // Live inputs for this law on this patch: the patch values,
                // which are a different set from the internal faces'
                const mechanicalConstitutiveLawInputs patchInputs
                (
                    lawInputsPatch(lawI, patchI, faces, dt, tp)
                );

                patchInputs.setConvergenceScale(scales[lawI]);

                // "View" into the kinematic and stress fields for this
                // material => does not copy data
                const smallStrainKinematicsViews views
                (
                    smallStrainKinematicsFields
                    (
                        gradD.boundaryField()[patchI],
                        gradD0.boundaryField()[patchI]
                    ),
                    faces
                );

                UIndirectList<symmTensor> stressView
                (
                    Foam::boundaryFieldRef(stress)[patchI], faces
                );

                // No fourth-order tangent is computed on this boundary,
                // so a request for one becomes a request for nothing. The
                // law must not be told a fourth-order tangent is wanted
                // when there is nowhere to put it
                const tangentRequest boundaryReq =
                    scalarTangentPtr && needsScalarTangent(tangentReq)
                  ? tangentReq
                  : tangentRequest::none;

                evaluateResponse
                (
                    laws_[lawI],
                    views.kin,
                    patchInputs,
                    tp.boundaryStates_[lawI][patchI],
                    stressView,
                    faces,
                    scalarTangentPtr
                  ? &scalarTangentPtr->boundaryField()[patchI]
                  : nullptr,
                    static_cast<const UList<mat66>*>(nullptr),
                    boundaryReq
                );
            }
        }
    }

    // Collapse accumulated stress on internal faces

    forAll(stress.internalField(), faceI)
    {
        const scalar w = weightSum[faceI];

        if (w <= SMALL)
        {
            FatalErrorInFunction
                << "Face " << faceI << " received no constitutive contributions"
                << exit(FatalError);
        }

        checkCollapsePermitted(collapseRule, w, faceI, "Face");

        // Stress collapse as arithmetic mean. With 'none' the weight is one,
        // so this is the single contribution unchanged
        stress[faceI] = stressSum[faceI]/w;

        // Tangent collapse, if requested
        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            if (collapseRule == stressCollapseRule::harmonic)
            {
                // Harmonically average the tangent
                (*scalarTangentPtr)[faceI] =
                    w/max((*tangentWeightPtr)[faceI], SMALL);
            }
            else
            {
                // Arithmetic mean of the contributing tangents
                (*scalarTangentPtr)[faceI] = (*tangentWeightPtr)[faceI]/w;
            }
        }
    }

// This one guard is real, and only this one. OpenFOAM.org's fvsPatchField
// has no evaluate(), so correctBoundaryConditions() does not compile for a
// SURFACE field there. It compiles and is needed for volFields
#ifndef OPENFOAM_ORG
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
    }
#endif
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const pointTensorField& gradD,
    const pointTensorField& gradD0,
    const scalar dt,
    pointSymmTensorField& stress,
    const stressCollapseRule collapseRule,
    pointScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Check gradD is defined on the correct mesh
    checkMeshConsistency(mesh_, gradD.mesh()(), gradD.name());
    checkMeshConsistency(mesh_, gradD0.mesh()(), gradD0.name());
    checkMeshConsistency(mesh_, stress.mesh()(), stress.name());
    if (scalarTangentPtr)
    {
        checkMeshConsistency
        (
            mesh_, scalarTangentPtr->mesh()(), scalarTangentPtr->name()
        );
    }

    // Topology + state

    const integrationPointTopology& topo =
        topologyFor(pointCentredIntegrationPointTopology::typeName);

    checkTangentRequest(topo, tangentReq);

    topologyEntry& tp = topology(topo);

    updateOldTimeIfNeeded();

    // Every processor refreshes every coupling input source here, before any
    // material is skipped for having no points on it, since reading is
    // collective
    refreshScalarInputs();

    // Accumulation fields

    pointSymmTensorField& stressSum = pointStressSum();
    pointScalarField& weightSum = pointStressWeight();

    stressSum = dimensionedSymmTensor("0", dimPressure, symmTensor::zero);
    weightSum = 0.0;

    pointScalarField* tangentWeightPtr = nullptr;

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        pointScalarField& tangentWeight = pointTangentWeight();
        tangentWeight = dimensionedScalar("zero", dimPressure, 0.0);
        tangentWeightPtr = &tangentWeight;
    }

    // Constitutive evaluation

    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Live inputs for this law, as a view of its own integration points
        const mechanicalConstitutiveLawInputs inputs
        (
            lawInputs(lawI, topo, ipIDs, dt, tp)
        );

        const smallStrainKinematicsViews views
        (
            smallStrainKinematicsFields
            (
                gradD.internalField(), gradD0.internalField()
            ),
            ipIDs
        );

        UIndirectList<symmTensor> stressView
        (
            stress.internalField(), ipIDs
        );

        // This path computes no fourth-order tangent, so none is offered
        evaluateResponse
        (
            laws_[lawI],
            views.kin,
            inputs,
            tp.states_[lawI],
            stressView,
            ipIDs,
            scalarTangentPtr
          ? &scalarTangentPtr->internalField()
          : nullptr,
            static_cast<const UList<mat66>*>(nullptr),
            tangentReq
        );

        // Accumulate per point

        forAll(ipIDs, i)
        {
            const label pointI = ipIDs[i];

            stressSum[pointI] += stress[pointI];
            weightSum[pointI] += 1.0;

            if (tangentWeightPtr && scalarTangentPtr)
            {
                const scalar K = (*scalarTangentPtr)[pointI];

                if (collapseRule == stressCollapseRule::harmonic)
                {
                    (*tangentWeightPtr)[pointI] += 1.0/max(K, SMALL);
                }
                else
                {
                    // 'average', and 'none', which is only reached with a
                    // single contribution and so is the same sum
                    (*tangentWeightPtr)[pointI] += K;
                }
            }
        }
    }

    // A point on a processor boundary is held by every rank that touches it,
    // and each rank has accumulated only what it can see, so the collapsed
    // value differs between ranks and matches a serial run on none of them.
    //
    // Summing the accumulators across ranks is NOT the fix, and was tried:
    // this topology asks for unique integration points per material, so a
    // point receives one contribution per material present on the rank rather
    // than one per cell. Summing then counts a material once for every rank
    // that holds it. Take a point shared by materials A and B with stresses
    // 10 and 40. In serial the collapse is 25. Decomposed so that rank 0 has
    // both materials and rank 1 has A as well, summing gives (10 + 40 + 10)/3
    // = 20 - a different wrong answer, and one that changes with the
    // decomposition just as the unsynchronised version does.
    //
    // What this needs is per-material contributions and multiplicities
    // reconciled across ranks before the materials are combined, or weights
    // that are genuinely additive. Neither is written, and nothing calls
    // either collapse overload today - no solid model does, and the
    // vertex-centred models evaluate on dual-mesh faces with no collapse - so
    // this refuses rather than guessing
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Point-centred stress collapse is not implemented for a "
            << "decomposed mesh." << nl << nl
            << "    A point on a processor boundary is held by several ranks "
            << "and each has accumulated only its own contributions, so the "
            << "collapsed value would depend on the decomposition. Summing "
            << "the accumulators does not fix it: this topology contributes "
            << "once per material rather than once per cell, so a sum counts "
            << "each material once per rank that holds it." << nl << nl
            << "    Run this case in serial, or use a topology that does not "
            << "collapse."
            << exit(FatalError);
    }

    // Collapse
    //
    // Both rules divide plain sums: the harmonic rule sums reciprocals, so
    // what makes it harmonic is what was accumulated rather than how it is
    // reduced. That is why one sum per accumulator serves both

    forAll(stress.internalField(), pointI)
    {
        const scalar w = weightSum[pointI];

        if (w <= SMALL)
        {
            FatalErrorInFunction
                << "Point " << pointI
                << " received no constitutive contributions"
                << exit(FatalError);
        }

        checkCollapsePermitted(collapseRule, w, pointI, "Point");

        // With 'none' the weight is one, so this is the single contribution
        // unchanged
        stress[pointI] = stressSum[pointI]/w;

        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            if (collapseRule == stressCollapseRule::harmonic)
            {
                (*scalarTangentPtr)[pointI] =
                    w/max((*tangentWeightPtr)[pointI], SMALL);
            }
            else
            {
                (*scalarTangentPtr)[pointI] =
                    (*tangentWeightPtr)[pointI]/w;
            }
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const CompactListList<tensor>& gradD,
    const CompactListList<tensor>& gradD0,
    const scalar dt,
    CompactListList<symmTensor>& stress,
    List<scalar>* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Check field sizes are consistent
    checkCompactLayoutConsistency
    (
        gradD,
        gradD0,
        stress,
        scalarTangentPtr,
        "updateStressSmallStrain (CompactListList)"
    );

    // Look up the map and state for compact list cell-based topologies
    const integrationPointTopology& topo = compactCellTopologyFor(gradD);

    // Evaluate on the packed integration-point storage. A compact cell
    // topology maps each cell to its own integration points, so the packed
    // lists are already in the flat form the primitive expects
    updateStressSmallStrain
    (
        topo,
        gradD.m(),
        gradD0.m(),
        dt,
        stress.m(),
        scalarTangentPtr,
        nullptr,
        tangentReq
    );
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
(
    const volTensorField& F,
    const volTensorField& F0,
    const volScalarField& J,
    const volScalarField& J0,
    const volTensorField& Finv,
    const volTensorField& Finv0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq,
    volScalarField* volumetricResponsePtr
)
{
    // Check F is defined on the correct mesh
    checkMeshConsistency(mesh_, F.mesh(), F.name());
    checkMeshConsistency(mesh_, F0.mesh(), F0.name());
    checkMeshConsistency(mesh_, Finv.mesh(), Finv.name());
    checkMeshConsistency(mesh_, Finv0.mesh(), Finv0.name());
    checkMeshConsistency(mesh_, J.mesh(), J.name());
    checkMeshConsistency(mesh_, J0.mesh(), J0.name());
    checkMeshConsistency(mesh_, stress.mesh(), stress.name());
    if (scalarTangentPtr)
    {
        checkMeshConsistency
        (
            mesh_, scalarTangentPtr->mesh(), scalarTangentPtr->name()
        );
    }

    if (volumetricResponsePtr)
    {
        checkMeshConsistency
        (
            mesh_,
            volumetricResponsePtr->mesh(),
            volumetricResponsePtr->name()
        );

        // The split is reachable through this overload as well as through
        // updateStressFiniteStrainSplit, so the capability is checked where
        // it is asked for rather than at one of the ways of asking. A law
        // that cannot separate the two would otherwise ignore the request and
        // hand back a total stress the caller would treat as isochoric
        checkVolumetricSplitSupported("updateStressFiniteStrain");
    }

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Look up the map and state for cell-based topologies
    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    // Update the internal field via the flat-list primitive: a cell-centred
    // topology has one integration point per cell, so the internal fields are
    // already in the flat form it expects
    evaluateFiniteStrain
    (
        topo,
        Foam::primitiveField(F),
        Foam::primitiveField(F0),
        Foam::primitiveField(Finv),
        Foam::primitiveField(Finv0),
        Foam::primitiveField(J),
        Foam::primitiveField(J0),
        dt,
        Foam::primitiveFieldRef(stress),
        scalarTangentPtr
      ? &Foam::primitiveFieldRef(*scalarTangentPtr)
      : nullptr,
        nullptr,
        tangentReq,
        false,          // commit the constitutive state
        volumetricResponsePtr
      ? &Foam::primitiveFieldRef(*volumetricResponsePtr)
      : nullptr
    );

    topologyEntry& tp = topology(topo);

    // The boundary faces, each with its own state
    updateStressVolBoundary<finiteStrainKinematicsFields>
    (
        tp,
        F,
        [&](const label patchI)
        {
            return finiteStrainKinematicsFields
            (
                F.boundaryField()[patchI],
                F0.boundaryField()[patchI],
                Finv.boundaryField()[patchI],
                Finv0.boundaryField()[patchI],
                J.boundaryField()[patchI],
                J0.boundaryField()[patchI]
            );
        },
        dt,
        stress,
        scalarTangentPtr,
        tangentReq,
        volumetricResponsePtr
    );

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
    }

    // The volumetric response is corrected with the others. Coupled patches
    // are skipped when it is evaluated, exactly as the stress is, so without
    // this its processor values are whatever the field was constructed with.
    // Only the internal field is read today, which makes that harmless rather
    // than correct - and harmless-for-now is a poor thing to leave for
    // whoever next reads a boundary value
    if (volumetricResponsePtr)
    {
        volumetricResponsePtr->correctBoundaryConditions();
    }
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
(
    const CompactListList<tensor>& F,
    const CompactListList<tensor>& F0,
    const CompactListList<tensor>& Finv,
    const CompactListList<tensor>& Finv0,
    const CompactListList<scalar>& J,
    const CompactListList<scalar>& J0,
    const scalar dt,
    CompactListList<symmTensor>& stress,
    List<scalar>* scalarTangentPtr,
    const tangentRequest tangentReq
)
{
    // Check field sizes are consistent
    const word context("updateStressFiniteStrain (CompactListList)");

    checkCompactLayoutConsistency(F, F0, stress, scalarTangentPtr, context);

    // The four the check above does not see. They are addressed by the same
    // flat index as F, so a different row shape puts one entity's inverse or
    // Jacobian against another's deformation gradient
    const labelList refRows(stress.sizes());

    checkCompactRowSizes(refRows, "the stress", Finv.sizes(), "Finv", context);
    checkCompactRowSizes
    (
        refRows, "the stress", Finv0.sizes(), "Finv0", context
    );
    checkCompactRowSizes(refRows, "the stress", J.sizes(), "J", context);
    checkCompactRowSizes(refRows, "the stress", J0.sizes(), "J0", context);

    // Look up the map and state for compact list cell-based topologies
    const integrationPointTopology& topo = compactCellTopologyFor(F);

    // Evaluate on the packed integration-point storage. A compact cell
    // topology maps each cell to its own integration points, so the packed
    // lists are already in the flat form the primitive expects
    updateStressFiniteStrain
    (
        topo,
        F.m(),
        F0.m(),
        Finv.m(),
        Finv0.m(),
        J.m(),
        J0.m(),
        dt,
        stress.m(),
        scalarTangentPtr,
        nullptr,
        tangentReq
    );
}


void Foam::mechanicalConstitutiveLawManager::checkVolumetricSplitSupported
(
    const word& context
) const
{
    forAll(laws_, lawI)
    {
        if (!laws_[lawI].providesVolumetricSplit())
        {
            FatalErrorInFunction
                << context << " asked for the isochoric stress and the "
                << "volumetric response separately, and the law for material '"
                << lawNames_[lawI] << "', of type " << laws_[lawI].type()
                << ", cannot supply them." << nl << nl
                << "    A mixed displacement-pressure formulation replaces the "
                << "volumetric part of the stress with a solved pressure, so "
                << "it needs the law's isochoric stress rather than the "
                << "trace-free part of its total. The two are the same only "
                << "when the law's energy is written on an isochoric measure; "
                << "otherwise they differ, and the difference is a different "
                << "material rather than a visible error." << nl << nl
                << "    Either use a law that separates them, or solve this "
                << "case in the displacement formulation."
                << exit(FatalError);
        }
    }
}


bool Foam::mechanicalConstitutiveLawManager::anyLawIncompressible() const
{
    forAll(laws_, lawI)
    {
        if (incompressibleLawTree(laws_[lawI]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::mechanicalConstitutiveLawManager::allLawsProvideVolumetricSplit
() const
{
    forAll(laws_, lawI)
    {
        if (!laws_[lawI].providesVolumetricSplit())
        {
            return false;
        }
    }

    return true;
}


bool Foam::mechanicalConstitutiveLawManager::
allLawsHaveDilationInvariantIsochoricStress() const
{
    forAll(laws_, lawI)
    {
        if (!laws_[lawI].isochoricStressIsDilationInvariant())
        {
            return false;
        }
    }

    return true;
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrain
(
    const surfaceTensorField& F,
    const surfaceTensorField& F0,
    const surfaceScalarField& J,
    const surfaceScalarField& J0,
    const surfaceTensorField& Finv,
    const surfaceTensorField& Finv0,
    const scalar dt,
    surfaceSymmTensorField& stress,
    const stressCollapseRule collapseRule,
    surfaceScalarField* volumetricResponsePtr
)
{
    // The face twin of the cell-centred finite-strain overload, laid out as
    // the small-strain face overload is: a face can be reached by two laws at
    // a material interface, so contributions are accumulated and collapsed
    // rather than written once
    //
    // In parallel a material interface can lie on a processor boundary, where
    // each side holds only its own material's contribution and nothing
    // exchanges them before the collapse, so the two sides would disagree
    // with each other and with serial. Refused rather than answered wrongly:
    // the only caller is single-material
    if (Pstream::parRun() && laws_.size() > 1)
    {
        FatalErrorInFunction
            << "The face finite-strain stress update does not support more "
            << "than one material in parallel: contributions at a material "
            << "interface on a processor boundary are not reconciled"
            << exit(FatalError);
    }
    checkMeshConsistency(mesh_, F.mesh(), F.name());
    checkMeshConsistency(mesh_, F0.mesh(), F0.name());
    checkMeshConsistency(mesh_, J.mesh(), J.name());
    checkMeshConsistency(mesh_, J0.mesh(), J0.name());
    checkMeshConsistency(mesh_, Finv.mesh(), Finv.name());
    checkMeshConsistency(mesh_, Finv0.mesh(), Finv0.name());
    checkMeshConsistency(mesh_, stress.mesh(), stress.name());

    if (volumetricResponsePtr)
    {
        checkMeshConsistency
        (
            mesh_,
            volumetricResponsePtr->mesh(),
            volumetricResponsePtr->name()
        );

        // Checked where the split is asked for, as the cell-centred overload
        // does: a law that cannot separate the two would otherwise hand back
        // a total stress the caller would treat as isochoric
        checkVolumetricSplitSupported("updateStressFiniteStrain");
    }

    // Look up the map and state for face-based topologies
    const integrationPointTopology& topo =
        topologyFor(faceCentredIntegrationPointTopology::typeName);

    topologyEntry& tp = topology(topo);

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Every processor refreshes every coupling input source here, before any
    // material is skipped for having no points on it, since reading is
    // collective
    refreshScalarInputs();

    // The scale each law's convergence test is normalised by, taken over its
    // internal faces on every rank, and used on its boundary faces too. It
    // was read from whatever a flat-list call on this topology had left,
    // which is usually nothing and otherwise another evaluation's
    const scalarList scales
    (
        convergenceScales
        (
            tp,
            finiteStrainKinematicsFields
            (
                F.internalField(),
                F0.internalField(),
                Finv.internalField(),
                Finv0.internalField(),
                J.internalField(),
                J0.internalField()
            )
        )
    );

    surfaceSymmTensorField& stressSum = surfaceStressSum();
    surfaceScalarField& weightSum = surfaceStressWeight();

    stressSum = dimensionedSymmTensor("0", dimPressure, symmTensor::zero);
    weightSum = 0.0;

    // The volumetric response is collapsed as the stress is. Only allocated
    // when the split is asked for
    autoPtr<scalarField> volumetricSumPtr;

    if (volumetricResponsePtr)
    {
        volumetricSumPtr.reset(new scalarField(mesh_.nInternalFaces(), 0.0));
    }

    checkTangentRequest(topo, tangentRequest::none);

    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        // Live inputs for this law, as a view of its own integration points
        const mechanicalConstitutiveLawInputs inputs
        (
            lawInputs(lawI, topo, ipIDs, dt, tp)
        );

        inputs.setConvergenceScale(scales[lawI]);

        const finiteStrainKinematicsViews views
        (
            finiteStrainKinematicsFields
            (
                F.internalField(),
                F0.internalField(),
                Finv.internalField(),
                Finv0.internalField(),
                J.internalField(),
                J0.internalField()
            ),
            ipIDs
        );

        UIndirectList<symmTensor> stressView(stress.internalField(), ipIDs);

        evaluateResponse
        (
            laws_[lawI],
            views.kin,
            inputs,
            tp.states_[lawI],
            stressView,
            ipIDs,
            static_cast<const UList<scalar>*>(nullptr),
            static_cast<const UList<mat66>*>(nullptr),
            tangentRequest::none,
            volumetricResponsePtr
        );

        forAll(ipIDs, i)
        {
            const label faceI = ipIDs[i];

            stressSum[faceI] += stress[faceI];
            weightSum[faceI] += 1.0;

            if (volumetricResponsePtr)
            {
                volumetricSumPtr()[faceI] += (*volumetricResponsePtr)[faceI];
            }
        }

        // Boundary faces carry their own state, as in the other overloads.
        // Coupled patches are included: unlike a volField, a surface field's
        // coupled patch is not filled in by correctBoundaryConditions(), so
        // skipping it would leave the stress there at zero
        if (tp.boundaryAware_)
        {
            forAll(F.boundaryField(), patchI)
            {
                const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                if
                (
                    faces.empty()
                 || isA<emptyFvPatch>(mesh_.boundary()[patchI])
                )
                {
                    continue;
                }

                // Live inputs for this law on this patch: the patch values,
                // which are a different set from the internal faces', judged
                // by the same scale
                const mechanicalConstitutiveLawInputs patchInputs
                (
                    lawInputsPatch(lawI, patchI, faces, dt, tp)
                );

                patchInputs.setConvergenceScale(scales[lawI]);

                const finiteStrainKinematicsViews views
                (
                    finiteStrainKinematicsFields
                    (
                        F.boundaryField()[patchI],
                        F0.boundaryField()[patchI],
                        Finv.boundaryField()[patchI],
                        Finv0.boundaryField()[patchI],
                        J.boundaryField()[patchI],
                        J0.boundaryField()[patchI]
                    ),
                    faces
                );

                UIndirectList<symmTensor> stressView
                (
                    Foam::boundaryFieldRef(stress)[patchI], faces
                );

                evaluateResponse
                (
                    laws_[lawI],
                    views.kin,
                    patchInputs,
                    tp.boundaryStates_[lawI][patchI],
                    stressView,
                    faces,
                    static_cast<const UList<scalar>*>(nullptr),
                    static_cast<const UList<mat66>*>(nullptr),
                    tangentRequest::none,
                    volumetricResponsePtr
                  ? &volumetricResponsePtr->boundaryField()[patchI]
                  : nullptr
                );
            }
        }
    }

    // Collapse the accumulated contributions on internal faces
    forAll(stress.internalField(), faceI)
    {
        const scalar w = weightSum[faceI];

        if (w <= SMALL)
        {
            FatalErrorInFunction
                << "Face " << faceI << " received no constitutive contributions"
                << exit(FatalError);
        }

        checkCollapsePermitted(collapseRule, w, faceI, "Face");

        stress[faceI] = stressSum[faceI]/w;

        if (volumetricResponsePtr)
        {
            (*volumetricResponsePtr)[faceI] = volumetricSumPtr()[faceI]/w;
        }
    }

// As in the small-strain face overload: OpenFOAM.org's fvsPatchField has no
// evaluate(), so correctBoundaryConditions() does not compile for a surface
// field there
#ifndef OPENFOAM_ORG
    stress.correctBoundaryConditions();

    if (volumetricResponsePtr)
    {
        volumetricResponsePtr->correctBoundaryConditions();
    }
#endif
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrainSplit
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
    const stressCollapseRule collapseRule
)
{
    updateStressFiniteStrain
    (
        F,
        F0,
        J,
        J0,
        Finv,
        Finv0,
        dt,
        isochoricStress,
        collapseRule,
        &volumetricResponse
    );
}


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrainSplit
(
    const volTensorField& gradD,
    const volTensorField& gradD0,
    const scalar dt,
    volSymmTensorField& deviatoricStress,
    volScalarField& volumetricResponse
)
{
    // The small-strain twin of updateStressFiniteStrainSplit, and the same
    // path as the ordinary update with somewhere to put the volumetric
    // response, so boundary faces are covered rather than left holding a
    // total stress while the interior holds a deviatoric one
    updateStressSmallStrain
    (
        gradD,
        gradD0,
        dt,
        deviatoricStress,
        nullptr,
        tangentRequest::none,
        &volumetricResponse
    );
}


void Foam::mechanicalConstitutiveLawManager::updateStressFiniteStrainSplit
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
)
{
    // The same path as the ordinary finite-strain update, with somewhere to
    // put the volumetric response. Sharing it is the point: the boundary
    // faces carry their own constitutive state and their own stress, and a
    // split that filled only the internal field would leave the boundary
    // holding a total stress while the interior held an isochoric one. The
    // traction is built from the boundary values, so that would be wrong
    // where it is least visible
    updateStressFiniteStrain
    (
        F,
        F0,
        J,
        J0,
        Finv,
        Finv0,
        dt,
        isochoricStress,
        nullptr,
        tangentRequest::none,
        &volumetricResponse
    );
}


void Foam::mechanicalConstitutiveLawManager::endTimeStepLaw
(
    const mechanicalConstitutiveLaw& law,
    mechanicalConstitutiveLawState& state,
    const scalar time,
    const label timeIndex,
    DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics,
    DynamicList<const mechanicalConstitutiveLaw*>* reporters,
    DynamicList<label>* firsts,
    DynamicList<label>* counts
) const
{
    const label before = diagnostics.size();

    law.endTimeStep(state, time, timeIndex, diagnostics);

    // Who said what. A composite's sub-law reports its own quantities, and
    // without this the wrapper would be handed them and drop them on the
    // floor - it inherits the empty reporting hook, so a nested plastic law's
    // numbers would be reduced and then silently discarded
    if (reporters && diagnostics.size() > before)
    {
        reporters->append(&law);
        firsts->append(before);
        counts->append(diagnostics.size() - before);
    }

    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        endTimeStepLaw
        (
            law.childLaw(childNames[i]),
            state.child(childNames[i]),
            time,
            timeIndex,
            diagnostics,
            reporters,
            firsts,
            counts
        );
    }
}



void Foam::mechanicalConstitutiveLawManager::writeStateFields() const
{
    const topologyEntry* entryPtr = nullptr;

    forAllIters(topologyEntries_, topoIter)
    {
        const topologyEntry& entry = autoPtrRef(topoIter());

        if (entry.topology_.integrationPointsAreCells())
        {
            entryPtr = &entry;
            break;
        }
    }

    if (!entryPtr)
    {
        return;
    }

    const topologyEntry& entry = *entryPtr;

    HashTable<stateFieldPlan> plans;

    forAll(laws_, lawI)
    {
        if (lawI >= entry.states_.size() || !entry.states_.set(lawI))
        {
            continue;
        }

        stateFieldPart part;
        part.internal = &entry.states_[lawI];
        part.cells = &entry.lawIntegrationPointIDs_[lawI];
        part.patchStates.setSize(mesh_.boundary().size(), nullptr);
        part.patchFaces.setSize(mesh_.boundary().size(), nullptr);

        if
        (
            entry.boundaryAware_
         && lawI < entry.boundaryStates_.size()
         && lawI < lawBoundaryFaces_.size()
        )
        {
            const PtrList<mechanicalConstitutiveLawState>& bStates =
                entry.boundaryStates_[lawI];

            forAll(bStates, patchI)
            {
                if (bStates.set(patchI) && patchI < part.patchStates.size())
                {
                    part.patchStates[patchI] = &bStates[patchI];
                    part.patchFaces[patchI] = &lawBoundaryFaces_[lawI][patchI];
                }
            }
        }

        collectStateFieldPlans(laws_[lawI], wordList(), part, plans);
    }

    if (plans.empty())
    {
        return;
    }

    // A variable declared at more than one level of a composite would write
    // two fields of the same name, so the deeper ones are qualified with the
    // path to them
    HashTable<label> nPaths;
    HashSet<word> seen;

    forAllIters(plans, iter)
    {
        const stateFieldPlan& plan = iter();

        word pathKey(plan.variable);
        forAll(plan.path, j)
        {
            pathKey += '/' + plan.path[j];
        }

        if (!seen.found(pathKey))
        {
            seen.insert(pathKey);

            if (nPaths.found(plan.variable))
            {
                nPaths[plan.variable]++;
            }
            else
            {
                nPaths.insert(plan.variable, 1);
            }
        }
    }

    HashSet<word> written;
    wordList names;

    wordList keys(plans.toc());
    sort(keys);

    forAll(keys, keyI)
    {
        const stateFieldPlan& plan = plans[keys[keyI]];

        word fieldName(plan.variable);

        if (nPaths[plan.variable] > 1 && plan.path.size())
        {
            fieldName.clear();
            forAll(plan.path, j)
            {
                fieldName += plan.path[j] + '_';
            }
            fieldName += plan.variable;
        }

        if (written.found(fieldName))
        {
            WarningInFunction
                << "State variable '" << plan.variable << "' is declared "
                << "with more than one type, so only one is written as '"
                << fieldName << "'" << endl;

            continue;
        }

        if (mesh_.foundObject<regIOobject>(fieldName))
        {
            if (!reportedUnreadInputs_.found("state:" + fieldName))
            {
                reportedUnreadInputs_.insert("state:" + fieldName);

                Info<< "    State variable '" << plan.variable
                    << "' is not written as a field: another model has "
                    << "registered '" << fieldName << "'" << endl;
            }

            continue;
        }

        if (plan.type == "scalar")
        {
            writeStateVolField<scalar>(mesh_, fieldName, plan);
        }
        else if (plan.type == "vector")
        {
            writeStateVolField<vector>(mesh_, fieldName, plan);
        }
        else if (plan.type == "tensor")
        {
            writeStateVolField<tensor>(mesh_, fieldName, plan);
        }
        else if (plan.type == "symmTensor")
        {
            writeStateVolField<symmTensor>(mesh_, fieldName, plan);
        }
        else
        {
            continue;
        }

        written.insert(fieldName);
        names.append(fieldName);
    }

    if (names.size())
    {
        Info<< "Writing constitutive state fields " << names << endl;
    }
}

void Foam::mechanicalConstitutiveLawManager::endTimeStep()
{
    const scalar time = mesh_.time().value();
    const label timeIndex = mesh_.time().timeIndex();

    // Gather, then reduce. The laws report locally and communicate nothing,
    // so this is where the collectives happen and where their count is made
    // the same on every rank.
    //
    // The shape comes from the internal state, which every rank has for every
    // law whether or not any of its integration points live here. A law
    // appends the same quantities in the same order every time it is asked,
    // so that first call fixes how many values there are and what they are
    // called; the boundary states then combine into those same slots. The
    // number of boundary states differs from rank to rank - processor patches
    // exist only where the mesh was cut - and that is exactly why they cannot
    // be allowed to influence how many collectives are called.
    //
    // The topologies are visited by name, in sorted order.
    //
    // Both tables are keyed by the name the topology was registered under,
    // which every rank supplies identically - a type name, a caller's word,
    // or the role a compact layout fills - so sorting them gives one order
    // that all ranks agree on. An address would not: it differs between
    // ranks, so sorting would make each rank's order stable without making
    // the orders agree
    wordList topologyKeys(topologyCache_.toc());
    Foam::sort(topologyKeys);

    // Agreeing on the order of the keys a rank has is not the same as having
    // the same keys. A topology is built on first use, so a rank that has not
    // reached that use holds one fewer. The reductions below are per topology,
    // so that rank calls fewer collectives than the others and the run stops
    // dead with no indication of why. Say what happened instead.
    //
    // The type registered under each key goes into the digest as well as the
    // key. registerTopology() takes the name from the caller, so two ranks can
    // supply the same word for different topologies; the keys would then agree
    // while the point sets behind them do not, and the reductions would pair
    // one rank's quantity with another's
    if (Pstream::parRun())
    {
        wordList keysAndTypes(2*topologyKeys.size());

        forAll(topologyKeys, keyI)
        {
            const word& key = topologyKeys[keyI];

            keysAndTypes[2*keyI] = key;
            keysAndTypes[2*keyI + 1] = autoPtrRef(topologyCache_[key]).type();
        }

        const label check = keyListChecksum(keysAndTypes);

        if
        (
            returnReduce(check, minOp<label>())
         != returnReduce(check, maxOp<label>())
        )
        {
            FatalErrorInFunction
                << "The ranks of this run do not hold the same set of "
                << "integration-point topologies." << nl << nl
                << "    This rank holds " << topologyKeys.size() << ": "
                << topologyKeys << nl << nl
                << "    The type registered under each key is compared too, "
                << "so this also fires where the names agree and the "
                << "topologies behind them do not." << nl << nl
                << "    Topologies are created when they are first used, so "
                << "this means some rank reached an evaluation that the "
                << "others did not. The end-of-step diagnostics reduce once "
                << "per topology, so the ranks would call different numbers "
                << "of collectives and the run would hang here rather than "
                << "report anything." << nl << nl
                << "    Every rank must register the same topologies, in the "
                << "sense of reaching the same evaluation calls, even where "
                << "it owns no cells of the material concerned."
                << exit(FatalError);
        }
    }

    forAll(topologyKeys, keyI)
    {
        topologyEntry& tp = topology(topologyCache_[topologyKeys[keyI]]());

        forAll(laws_, lawI)
        {
            DynamicList<mechanicalConstitutiveLawDiagnostic> diagnostics;
            DynamicList<const mechanicalConstitutiveLaw*> reporters;
            DynamicList<label> firsts;
            DynamicList<label> counts;

            // The law and every law below it. A composite's sub-law keeps its
            // own history and was never told the step had ended
            endTimeStepLaw
            (
                laws_[lawI],
                tp.states_[lawI],
                time,
                timeIndex,
                diagnostics,
                &reporters,
                &firsts,
                &counts
            );

            // Not skipped where the law reported nothing. endTimeStep() is
            // the hook a law commits or rolls over anything storeOldTime()
            // does not, and nothing in its interface ties that to also
            // reporting a diagnostic. Returning early here would call it on
            // the internal points and not on the boundary ones, which is a
            // wrong answer exactly where it is hardest to see. Only the
            // reduction and the reporting below depend on there being
            // anything to report

            // The boundary states
            if (tp.boundaryAware_ && lawI < tp.boundaryStates_.size())
            {
                forAll(tp.boundaryStates_[lawI], patchI)
                {
                    DynamicList<mechanicalConstitutiveLawDiagnostic> patchDiags;

                    // Every patch, coupled or not. This is the law's hook for
                    // committing its own end-of-step work, and a coupled
                    // patch's state is as real as any other where the
                    // topology evaluates one
                    endTimeStepLaw
                    (
                        laws_[lawI],
                        tp.boundaryStates_[lawI][patchI],
                        time,
                        timeIndex,
                        patchDiags,
                        nullptr,
                        nullptr,
                        nullptr
                    );

                    // A coupled patch is not counted, though. Where it is
                    // evaluated at all - the face-centred topology does,
                    // because a processor face carries a stress only this
                    // rank can compute - the face is held by both ranks, so
                    // counting it here would put the same faces in the total
                    // twice, once from each side, and make the total depend
                    // on how the mesh was cut: this case reports 380
                    // integration points in serial and would report 458 on
                    // four ranks
                    if (mesh_.boundary()[patchI].coupled())
                    {
                        continue;
                    }

                    if (patchDiags.size() != diagnostics.size())
                    {
                        FatalErrorInFunction
                            << "The law for material '" << lawNames_[lawI]
                            << "' reported " << patchDiags.size()
                            << " quantities on patch " << patchI << " and "
                            << diagnostics.size() << " on the internal points."
                            << nl << nl
                            << "    A law must report the same quantities in "
                            << "the same order every time it is asked, "
                            << "because the number of them fixes how many "
                            << "collectives are called and that has to match "
                            << "on every rank."
                            << exit(FatalError);
                    }

                    // Nothing to combine where the law reports nothing; the
                    // hook above has run either way
                    forAll(diagnostics, d)
                    {
                        // Names and operations too, not just the count. Two
                        // quantities that arrive in a different order, or one
                        // that is summed here and maximised there, would
                        // combine silently and give a number that means
                        // nothing
                        if
                        (
                            patchDiags[d].name() != diagnostics[d].name()
                         || patchDiags[d].op() != diagnostics[d].op()
                        )
                        {
                            FatalErrorInFunction
                                << "The law for material '" << lawNames_[lawI]
                                << "' reported '" << patchDiags[d].name()
                                << "' where the internal points reported '"
                                << diagnostics[d].name() << "'." << nl << nl
                                << "    A law must report the same quantities "
                                << "in the same order every time it is asked."
                                << exit(FatalError);
                        }

                        combineDiagnostic(diagnostics[d], patchDiags[d]);
                    }
                }
            }

            // One collective per quantity, in an order fixed by the law
            forAll(diagnostics, d)
            {
                mechanicalConstitutiveLawDiagnostic& diag = diagnostics[d];

                if
                (
                    diag.op()
                 == mechanicalConstitutiveLawDiagnostic::combineOperation::sum
                )
                {
                    reduce(diag.value(), sumOp<scalar>());
                }
                else
                {
                    reduce(diag.value(), maxOp<scalar>());
                }
            }

            // Each law that reported says what it did, and is handed back
            // only its own quantities. The manager knows how to combine the
            // numbers but not what they mean, and the switch that decides
            // whether anyone wants to hear it belongs to the law
            forAll(reporters, r)
            {
                const SubList<mechanicalConstitutiveLawDiagnostic> slice
                (
                    diagnostics, counts[r], firsts[r]
                );

                reporters[r]->reportDiagnostics(slice);
            }
        }
    }
}



// ************************************************************************* //
