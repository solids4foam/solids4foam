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
#include "compatibilityFunctions.H"
#include "integrationPointTopologies.H"
#include "emptyFvPatch.H"
#include "mat66.H"
#include "Switch.H"
#include "CompactListList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalConstitutiveLawManager, 0);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


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
        return topology(topologyCache_[topologyTypeName]()).topology_;
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

    return topology(topologyCache_[topologyTypeName]()).topology_;
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::compactCellTopologyFor
(
    const CompactListList<tensor>& layout
) const
{
    // Which entity indexes the rows is decided by the row count, and checked.
    // A mesh never has as many cells as faces, so the two cases cannot be
    // confused, and anything else is an error rather than a guess
    const bool cellBased = (layout.size() == mesh_.nCells());
    const bool faceBased = (layout.size() == mesh_.nFaces());

    if (!cellBased && !faceBased)
    {
        FatalErrorInFunction
            << "A compact integration-point layout must have one row per cell "
            << "or one row per face." << nl
            << "This one has " << layout.size() << " rows, while the mesh has "
            << mesh_.nCells() << " cells and " << mesh_.nFaces() << " faces."
            << exit(FatalError);
    }

    // Unique key per layout instance
    const word key =
        (cellBased ? "compactCell:" : "compactFace:") + Foam::name
        (
            static_cast<std::uint64_t>
            (
                reinterpret_cast<std::uintptr_t>(&layout)
            )
        );

    // Already constructed?
    if (topologyCache_.found(key))
    {
        return topology(topologyCache_[key]()).topology_;
    }

    // Construct topology lazily

    // We know this is cell-based compact storage:
    //  - one sub-list per cell
    //  - integration-point counts encoded in sub-list sizes

    // Build cell → IP addressing
    const labelList rowSizes(layout.sizes());

    CompactListList<label> cellToIP(rowSizes);

    for (label cellI = 0; cellI < layout.size(); ++cellI)
    {
        // sizes() rather than layout[cellI].size(): the const operator[] does
        // not compile on foam-extend
        const label n = rowSizes[cellI];
        for (label j = 0; j < n; ++j)
        {
            cellToIP(cellI, j) = layout.index(cellI, j);
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

    return topology(topologyCache_[key]()).topology_;
}


const Foam::integrationPointTopology&
Foam::mechanicalConstitutiveLawManager::compactFaceTopologyFor
(
    const CompactListList<tensor>& layout
) const
{
    FatalErrorInFunction
        << "compactFaceTopologyFor is not yet implemented"
        << exit(FatalError);

    // Dummy return to silence compiler warnings
    return topologyFor(cellCentredIntegrationPointTopology::typeName);
}


Foam::mechanicalConstitutiveLawManager::topologyEntry&
Foam::mechanicalConstitutiveLawManager::topology
(
    const integrationPointTopology& topo
) const
{
    // Use the address of the topology object as a unique key
    const word key = Foam::name
    (
        static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(&topo))
    );

    // Return existing entry if already constructed
    if (topologyEntries_.found(key))
    {
        return topologyEntries_[key]();
    }

    // ---------------------------------------------------------------------
    // Construct new topology entry (lazy initialisation)
    // ---------------------------------------------------------------------

    DebugInfo
        << "Creating topologyEntry for " << key << endl;

    autoPtr<topologyEntry> entryPtr(new topologyEntry(topo));
    topologyEntries_.insert(key, entryPtr);
    topologyEntry& entry = topologyEntries_[key]();

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

                applyStateSpecPatch
                (
                    lawI, patchI, entry.boundaryStates_[lawI][patchI]
                );
            }
        }
    }

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


Foam::mechanicalConstitutiveLawInputs
Foam::mechanicalConstitutiveLawManager::lawInputsPatch
(
    const label lawI,
    const label patchI,
    const labelList& faces,
    const scalar dt,
    topologyEntry& tp
) const
{
    mechanicalConstitutiveLawInputs inputs(dt);

    const wordList names(laws_[lawI].requiredScalarInputs());

    if (names.empty())
    {
        return inputs;
    }

    if (tp.lawScalarInputs_.empty())
    {
        tp.lawScalarInputs_.setSize(laws_.size());
    }

    HashTable<autoPtr<scalarField>>& store = tp.lawScalarInputs_[lawI];

    forAll(names, i)
    {
        const word& name = names[i];

        // The boundary values are a different set from the internal ones, so
        // they get storage of their own rather than sharing it and being
        // overwritten by whichever was gathered last
        const word key(name + "@patch");

        if (!store.found(key))
        {
            store.insert(key, autoPtr<scalarField>(new scalarField()));
        }

        scalarField& fld = store[key]();
        fld.setSize(faces.size(), 0.0);

        tmp<volScalarField> tsrc;

        if (mesh_.foundObject<volScalarField>(name))
        {
            tsrc = tmp<volScalarField>
            (
                mesh_.lookupObject<volScalarField>(name)
            );
        }
        else
        {
            tsrc = prescribedField<scalar>(name);
        }

        if (!tsrc.valid())
        {
            FatalErrorInFunction
                << "Mechanical constitutive law '" << laws_[lawI].type()
                << "' reads the coupling input '" << name
                << "', which is neither registered nor present as a field in "
                << "the current or initial time directory."
                << exit(FatalError);
        }

        const fvPatchField<scalar>& psrc = tsrc().boundaryField()[patchI];

        forAll(faces, faceI)
        {
            fld[faceI] = psrc[faces[faceI]];
        }

        inputs.setScalar(name, fld);
    }

    return inputs;
}


Foam::mechanicalConstitutiveLawInputs
Foam::mechanicalConstitutiveLawManager::inputsWithoutCoupling
(
    const scalar dt,
    const word& path
) const
{
    forAll(laws_, lawI)
    {
        const wordList names(laws_[lawI].requiredScalarInputs());

        if (!names.empty())
        {
            FatalErrorInFunction
                << "Mechanical constitutive law '" << laws_[lawI].type()
                << "' reads the coupling input(s) " << names
                << ", but the evaluation path '" << path
                << "' does not gather them yet." << nl
                << "Use a solid model that evaluates through a path which "
                << "does, or add the gather to this one."
                << exit(FatalError);
        }
    }

    return mechanicalConstitutiveLawInputs(dt);
}


Foam::mechanicalConstitutiveLawInputs
Foam::mechanicalConstitutiveLawManager::lawInputs
(
    const label lawI,
    const integrationPointTopology& topo,
    const labelList& ipIDs,
    const scalar dt,
    topologyEntry& tp
) const
{
    mechanicalConstitutiveLawInputs inputs(dt);

    const wordList names(laws_[lawI].requiredScalarInputs());

    if (names.empty())
    {
        // Which is every law that does not couple to another field, so this
        // costs them nothing
        return inputs;
    }

    if (tp.lawScalarInputs_.empty())
    {
        tp.lawScalarInputs_.setSize(laws_.size());
    }

    HashTable<autoPtr<scalarField>>& store = tp.lawScalarInputs_[lawI];

    forAll(names, i)
    {
        const word& name = names[i];

        if (!store.found(name))
        {
            store.insert(name, autoPtr<scalarField>(new scalarField()));
        }

        scalarField& fld = store[name]();
        fld.setSize(ipIDs.size(), 0.0);

        // The registry first, and only then the file. This is the opposite
        // order from a prescribed field, and deliberately so: a coupling input
        // is solved for, so the live field must always win and a file must
        // never shadow it.
        //
        // The file is still needed, because of when the first evaluation
        // happens. A solid model builds its implicit stiffness while it is
        // being constructed, and a derived model that solves for the coupling
        // field has not reached its own members yet - the base class runs
        // first. At that moment the field exists only as the initial condition
        // on disk, which is the right value to evaluate against anyway
        if (mesh_.foundObject<volScalarField>(name))
        {
            gatherToIntegrationPoints
            (
                mesh_.lookupObject<volScalarField>(name),
                lawI,
                topo,
                ipIDs,
                fld
            );
        }
        else
        {
            const tmp<volScalarField> tsrc(prescribedField<scalar>(name));

            if (!tsrc.valid())
            {
                FatalErrorInFunction
                    << "Mechanical constitutive law '" << laws_[lawI].type()
                    << "' reads the coupling input '" << name
                    << "', which is neither registered nor present as a field "
                    << "in the current or initial time directory." << nl
                    << "It is produced by another model, so the solid model "
                    << "has to solve for it, or the case has to supply it."
                    << exit(FatalError);
            }

            gatherToIntegrationPoints(tsrc(), lawI, topo, ipIDs, fld);
        }

        inputs.setScalar(name, fld);
    }

    return inputs;
}


void Foam::mechanicalConstitutiveLawManager::applyStateDefaults
(
    const mechanicalConstitutiveLawStateSpec& spec,
    mechanicalConstitutiveLawState& state
) const
{
    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        const mechanicalConstitutiveLawStateSpec::entry& e = es[i];

        if (e.typeName == "scalar")
        {
            state.scalarField(e.name) = e.scalarDefault;
            state.scalarField0(e.name) = e.scalarDefault;
        }
        else if (e.typeName == "vector")
        {
            state.vectorField(e.name) = e.vectorDefault;
            state.vectorField0(e.name) = e.vectorDefault;
        }
        else if (e.typeName == "tensor")
        {
            state.tensorField(e.name) = e.tensorDefault;
            state.tensorField0(e.name) = e.tensorDefault;
        }
        else if (e.typeName == "symmTensor")
        {
            state.symmTensorField(e.name) = e.symmTensorDefault;
            state.symmTensorField0(e.name) = e.symmTensorDefault;
        }
        else
        {
            FatalErrorInFunction
                << "State field '" << e.name << "' has unsupported type '"
                << e.typeName << "'." << exit(FatalError);
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::applyStateSpec
(
    const label lawI,
    const integrationPointTopology& topo,
    const labelList& ipIDs,
    mechanicalConstitutiveLawState& state
) const
{
    applyStateSpec(laws_[lawI], lawI, topo, ipIDs, state);
}


void Foam::mechanicalConstitutiveLawManager::applyStateSpec
(
    const mechanicalConstitutiveLaw& law,
    const label lawI,
    const integrationPointTopology& topo,
    const labelList& ipIDs,
    mechanicalConstitutiveLawState& state
) const
{
    // Order matters. The declared default goes down first, then the law is
    // given the chance to initialise state of its own over it, and only then
    // is a prescribed field read. A law that both declares a field and fills
    // it in initialiseState would otherwise have its work overwritten
    mechanicalConstitutiveLawStateSpec spec;
    law.declareState(spec);

    applyStateDefaults(spec, state);

    law.initialiseState(state);

    // A composite's sub-laws each get a child state of their own, prepared
    // exactly as this one was. Without this a sub-law's declared defaults and
    // prescribed fields would be missing and it would silently read zeros
    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        applyStateSpec
        (
            law.childLaw(childNames[i]),
            lawI,
            topo,
            ipIDs,
            state.child(childNames[i])
        );
    }

    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        const mechanicalConstitutiveLawStateSpec::entry& e = es[i];

        if (e.role != mechanicalConstitutiveLawStateSpec::stateRole::prescribed)
        {
            continue;
        }

        // A prescribed field is read from a volField the user supplied. Its
        // absence is not an error: the declared default stands, which is what
        // lets a law declare one without changing any existing case
        if (e.typeName == "scalar")
        {
            readPrescribed<scalar>
            (
                e.name, lawI, topo, ipIDs, state.scalarField(e.name)
            );
            state.scalarField0(e.name) = state.scalarField(e.name);
        }
        else if (e.typeName == "vector")
        {
            readPrescribed<vector>
            (
                e.name, lawI, topo, ipIDs, state.vectorField(e.name)
            );
            state.vectorField0(e.name) = state.vectorField(e.name);
        }
        else if (e.typeName == "tensor")
        {
            readPrescribed<tensor>
            (
                e.name, lawI, topo, ipIDs, state.tensorField(e.name)
            );
            state.tensorField0(e.name) = state.tensorField(e.name);
        }
        else
        {
            readPrescribed<symmTensor>
            (
                e.name, lawI, topo, ipIDs, state.symmTensorField(e.name)
            );
            state.symmTensorField0(e.name) = state.symmTensorField(e.name);
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::applyStateSpecScratch
(
    const label lawI,
    mechanicalConstitutiveLawState& state
) const
{
    applyStateSpecScratch(laws_[lawI], state);
}


void Foam::mechanicalConstitutiveLawManager::applyStateSpecScratch
(
    const mechanicalConstitutiveLaw& law,
    mechanicalConstitutiveLawState& state
) const
{
    mechanicalConstitutiveLawStateSpec spec;
    law.declareState(spec);

    applyStateDefaults(spec, state);

    law.initialiseState(state);

    // A composite's sub-laws each get a child state of their own, prepared
    // exactly as this one was. Without this a sub-law's declared defaults and
    // prescribed fields would be missing and it would silently read zeros
    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        applyStateSpecScratch
        (
            law.childLaw(childNames[i]),
            state.child(childNames[i])
        );
    }
}


void Foam::mechanicalConstitutiveLawManager::applyStateSpecPatch
(
    const label lawI,
    const label patchI,
    mechanicalConstitutiveLawState& state
) const
{
    applyStateSpecPatch(laws_[lawI], lawI, patchI, state);
}


void Foam::mechanicalConstitutiveLawManager::applyStateSpecPatch
(
    const mechanicalConstitutiveLaw& law,
    const label lawI,
    const label patchI,
    mechanicalConstitutiveLawState& state
) const
{
    mechanicalConstitutiveLawStateSpec spec;
    law.declareState(spec);

    applyStateDefaults(spec, state);

    law.initialiseState(state);

    // A composite's sub-laws each get a child state of their own, prepared
    // exactly as this one was. Without this a sub-law's declared defaults and
    // prescribed fields would be missing and it would silently read zeros
    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        applyStateSpecPatch
        (
            law.childLaw(childNames[i]),
            lawI,
            patchI,
            state.child(childNames[i])
        );
    }

    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        const mechanicalConstitutiveLawStateSpec::entry& e = es[i];

        if (e.role != mechanicalConstitutiveLawStateSpec::stateRole::prescribed)
        {
            continue;
        }

        if (e.typeName == "scalar")
        {
            readPrescribedPatch<scalar>
            (
                e.name, lawI, patchI, state.scalarField(e.name)
            );
            state.scalarField0(e.name) = state.scalarField(e.name);
        }
        else if (e.typeName == "vector")
        {
            readPrescribedPatch<vector>
            (
                e.name, lawI, patchI, state.vectorField(e.name)
            );
            state.vectorField0(e.name) = state.vectorField(e.name);
        }
        else if (e.typeName == "tensor")
        {
            readPrescribedPatch<tensor>
            (
                e.name, lawI, patchI, state.tensorField(e.name)
            );
            state.tensorField0(e.name) = state.tensorField(e.name);
        }
        else
        {
            readPrescribedPatch<symmTensor>
            (
                e.name, lawI, patchI, state.symmTensorField(e.name)
            );
            state.symmTensorField0(e.name) = state.symmTensorField(e.name);
        }
    }
}


void Foam::mechanicalConstitutiveLawManager::updateOldTimeIfNeeded()
{
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
        forAllIter
        (
            HashTable<autoPtr<topologyEntry>>, topologyEntries_, topoIter
        )
        {
            topologyEntry& entry = topoIter()();

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
    topologyEntries_()
{
    // Read the mechanical laws
    const PtrList<entry> lawEntries(dict.lookup("mechanical"));

    laws_.setSize(lawEntries.size());
    lawCells_.setSize(lawEntries.size());
    lawBoundaryFaces_.setSize(lawEntries.size());

    // Look list of law names
    wordList lawNames(lawEntries.size());
    forAll(lawNames, lawI)
    {
        lawNames[lawI] = lawEntries[lawI].keyword();
    }

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

        if (laws_.size() == 1)
        {
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

    // Set lawBoundaryFaces
    forAll(lawNames, lawI)
    {
        lawBoundaryFaces_[lawI].resize(mesh.boundary().size());

        forAll(lawBoundaryFaces_[lawI], patchI)
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            labelHashSet curFaceSet;

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];

                if (cellToLaw[cellID] == lawI)
                {
                    curFaceSet.insert(faceI);
                }
            }

            lawBoundaryFaces_[lawI][patchI] = curFaceSet.toc();
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawManager::~mechanicalConstitutiveLawManager()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const Foam::volScalarField& Foam::mechanicalConstitutiveLawManager::rho() const
{
    if (!rhoPtr_.valid())
    {
        rhoPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
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
                dimensionedScalar("kappa", dimDensity, 0.0),
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
        const integrationPointTopology& existing = topologyCache_[key]();

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

    return topology(topologyCache_[key]()).topology_;
}


void Foam::mechanicalConstitutiveLawManager::resetMaterialPropertyFields()
{
    rhoPtr_.clear();
    kappaPtr_.clear();
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
    const bool preserveState
)
{
    const word context = "updateStressSmallStrain (flat list)";
    const label nIP = topo.nIntegrationPoints();

    checkIntegrationPointListSize(nIP, gradD.size(), "gradD", context);
    checkIntegrationPointListSize(nIP, gradD0.size(), "gradD0", context);
    checkIntegrationPointListSize(nIP, stress.size(), "stress", context);

    checkTangentRequest(topo, tangentReq);

    checkTangentStorage
    (
        scalarTangentPtr != nullptr,
        fourthOrderTangentPtr != nullptr,
        tangentReq,
        context
    );

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
            << "material interface would be written more than once." << nl
            << "Use the surfaceField or pointField overload, which takes a "
            << "stressCollapseRule."
            << exit(FatalError);
    }

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    topologyEntry& tp = topology(topo);

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

        // A tangent query evaluates against a shadow of the law's state: the
        // shadow aliases the old-time fields, so history is read but never
        // written, and the law's outputs land where they are discarded
        autoPtr<mechanicalConstitutiveLawState> shadowPtr;
        if (preserveState)
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
            preserveState ? shadowPtr() : tp.states_[lawI];

        // Views into integration-point data (no copies)
        const UIndirectList<tensor> gradDView(gradD, ipIDs);
        const UIndirectList<tensor> gradD0View(gradD0, ipIDs);
        UIndirectList<symmTensor> stressView(stress, ipIDs);

        // Kinematics wrapper
        smallStrainMechanicalConstitutiveLawKinematics kin
        (
            gradDView, gradD0View
        );

        // Constitutive response
        if (needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView(*scalarTangentPtr, ipIDs);

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, lawState, response);
        }
        else if (needsFourthOrderTangent(tangentReq))
        {
            UIndirectList<mat66> tangentView(*fourthOrderTangentPtr, ipIDs);

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, lawState, response);
        }
        else
        {
            mechanicalConstitutiveLawResponse response(stressView, tangentReq);

            laws_[lawI].evaluate(kin, inputs, lawState, response);
        }
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
                const labelList patchIPs
                (
                    topo.boundaryIntegrationPointIDs(patchI)
                );

                if (patchIPs.empty())
                {
                    continue;
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

                labelList ipIDs;

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
                    continue;
                }

                // Live inputs for this law on this patch. The boundary points
                // are a different set from the internal ones, so the coupling
                // input has to be gathered for them rather than reused
                const mechanicalConstitutiveLawInputs inputs
                (
                    lawInputs(lawI, topo, ipIDs, dt, tp)
                );

                // An empty patch has no per-patch state: boundaryStates_ is
                // sized to the fvPatch, which has no faces there. Its faces
                // take no part in the discretisation, so give them a scratch
                // state of the right size rather than a shadow of a state
                // that is the wrong length. Nothing reads it afterwards
                autoPtr<mechanicalConstitutiveLawState> bScratchPtr;
                autoPtr<mechanicalConstitutiveLawState> bShadowPtr;

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

                mechanicalConstitutiveLawState& bState =
                    viaPolyPatch
                  ? bScratchPtr()
                  : preserveState
                  ? bShadowPtr()
                  : tp.boundaryStates_[lawI][patchI];

                const UIndirectList<tensor> gradDView(gradD, ipIDs);
                const UIndirectList<tensor> gradD0View(gradD0, ipIDs);
                UIndirectList<symmTensor> stressView(stress, ipIDs);

                smallStrainMechanicalConstitutiveLawKinematics kin
                (
                    gradDView, gradD0View
                );

                if (needsScalarTangent(tangentReq))
                {
                    UIndirectList<scalar> tanView(*scalarTangentPtr, ipIDs);

                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tanView, tangentReq
                    );

                    laws_[lawI].evaluate(kin, inputs, bState, response);
                }
                else if (needsFourthOrderTangent(tangentReq))
                {
                    UIndirectList<mat66> tanView
                    (
                        *fourthOrderTangentPtr, ipIDs
                    );

                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tanView, tangentReq
                    );

                    laws_[lawI].evaluate(kin, inputs, bState, response);
                }
                else
                {
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentReq
                    );

                    laws_[lawI].evaluate(kin, inputs, bState, response);
                }
            }
        }
    }
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
    const bool preserveState
)
{
    const word context = "updateStressFiniteStrain (flat list)";
    const label nIP = topo.nIntegrationPoints();

    checkIntegrationPointListSize(nIP, F.size(), "F", context);
    checkIntegrationPointListSize(nIP, F0.size(), "F0", context);
    checkIntegrationPointListSize(nIP, Finv.size(), "Finv", context);
    checkIntegrationPointListSize(nIP, Finv0.size(), "Finv0", context);
    checkIntegrationPointListSize(nIP, J.size(), "J", context);
    checkIntegrationPointListSize(nIP, J0.size(), "J0", context);
    checkIntegrationPointListSize(nIP, stress.size(), "stress", context);

    checkTangentRequest(topo, tangentReq);

    checkTangentStorage
    (
        scalarTangentPtr != nullptr,
        fourthOrderTangentPtr != nullptr,
        tangentReq,
        context
    );

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
            << exit(FatalError);
    }

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Live inputs for this evaluation. Built once and passed through
    // every evaluation, including each finite-difference perturbation,
    // so there is no per-call forwarding to get wrong
    const mechanicalConstitutiveLawInputs inputs
    (
        inputsWithoutCoupling(dt, "evaluateFiniteStrain")
    );

    topologyEntry& tp = topology(topo);

    // Loop over mechanical constitutive laws
    forAll(laws_, lawI)
    {
        const labelList& ipIDs = tp.lawIntegrationPointIDs_[lawI];

        if (ipIDs.empty())
        {
            continue;
        }

        // A tangent query evaluates against a shadow of the law's state: the
        // shadow aliases the old-time fields, so history is read but never
        // written, and the law's outputs land where they are discarded
        autoPtr<mechanicalConstitutiveLawState> shadowPtr;
        if (preserveState)
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
            preserveState ? shadowPtr() : tp.states_[lawI];

        // Views into integration-point data (no copies)
        const UIndirectList<tensor> FView(F, ipIDs);
        const UIndirectList<tensor> F0View(F0, ipIDs);
        const UIndirectList<tensor> FinvView(Finv, ipIDs);
        const UIndirectList<tensor> Finv0View(Finv0, ipIDs);
        const UIndirectList<scalar> JView(J, ipIDs);
        const UIndirectList<scalar> J0View(J0, ipIDs);
        UIndirectList<symmTensor> stressView(stress, ipIDs);

        // Kinematics wrapper
        finiteStrainMechanicalConstitutiveLawKinematics kin
        (
            FView,
            F0View,
            JView,
            J0View,
            FinvView,
            Finv0View
        );

        // Constitutive response
        if (needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView(*scalarTangentPtr, ipIDs);

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, lawState, response);
        }
        else if (needsFourthOrderTangent(tangentReq))
        {
            UIndirectList<mat66> tangentView(*fourthOrderTangentPtr, ipIDs);

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, lawState, response);
        }
        else
        {
            mechanicalConstitutiveLawResponse response(stressView, tangentReq);

            laws_[lawI].evaluate(kin, inputs, lawState, response);
        }
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
                const labelList patchIPs
                (
                    topo.boundaryIntegrationPointIDs(patchI)
                );

                if (patchIPs.empty())
                {
                    continue;
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

                labelList ipIDs;

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
                    continue;
                }

                // An empty patch has no per-patch state: boundaryStates_ is
                // sized to the fvPatch, which has no faces there. Its faces
                // take no part in the discretisation, so give them a scratch
                // state of the right size rather than a shadow of a state
                // that is the wrong length. Nothing reads it afterwards
                autoPtr<mechanicalConstitutiveLawState> bScratchPtr;
                autoPtr<mechanicalConstitutiveLawState> bShadowPtr;

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

                mechanicalConstitutiveLawState& bState =
                    viaPolyPatch
                  ? bScratchPtr()
                  : preserveState
                  ? bShadowPtr()
                  : tp.boundaryStates_[lawI][patchI];

                const UIndirectList<tensor> FView(F, ipIDs);
                const UIndirectList<tensor> F0View(F0, ipIDs);
                const UIndirectList<tensor> FinvView(Finv, ipIDs);
                const UIndirectList<tensor> Finv0View(Finv0, ipIDs);
                const UIndirectList<scalar> JView(J, ipIDs);
                const UIndirectList<scalar> J0View(J0, ipIDs);
                UIndirectList<symmTensor> stressView(stress, ipIDs);

                finiteStrainMechanicalConstitutiveLawKinematics kin
                (
                    FView, F0View, JView, J0View, FinvView, Finv0View
                );

                if (needsScalarTangent(tangentReq))
                {
                    UIndirectList<scalar> tanView(*scalarTangentPtr, ipIDs);

                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tanView, tangentReq
                    );

                    laws_[lawI].evaluate(kin, inputs, bState, response);
                }
                else if (needsFourthOrderTangent(tangentReq))
                {
                    UIndirectList<mat66> tanView
                    (
                        *fourthOrderTangentPtr, ipIDs
                    );

                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tanView, tangentReq
                    );

                    laws_[lawI].evaluate(kin, inputs, bState, response);
                }
                else
                {
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentReq
                    );

                    laws_[lawI].evaluate(kin, inputs, bState, response);
                }
            }
        }
    }
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
    const tangentRequest tangentReq
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
        true            // preserve the constitutive state
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
    const tangentRequest tangentReq
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
        tangentReq
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


void Foam::mechanicalConstitutiveLawManager::updateStressSmallStrain
(
    const volTensorField& gradD,
    const volTensorField& gradD0,
    const scalar dt,
    volSymmTensorField& stress,
    volScalarField* scalarTangentPtr,
    const tangentRequest tangentReq
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

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Look up the map and state for cell-based topologies
    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    // Update the internal field via the flat-list primitive: a cell-centred
    // topology has one integration point per cell, so the internal fields are
    // already in the flat form it expects
    updateStressSmallStrain
    (
        topo,
        Foam::primitiveField(gradD),
        Foam::primitiveField(gradD0),
        dt,
        Foam::primitiveFieldRef(stress),
        scalarTangentPtr ? &Foam::primitiveFieldRef(*scalarTangentPtr) : nullptr,
        nullptr,
        tangentReq
    );

    topologyEntry& tp = topology(topo);

    forAll(laws_, lawI)
    {
        // Optionally, update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        if (tp.boundaryAware_)
        {
            forAll(gradD.boundaryField(), patchI)
            {
                if (!gradD.boundaryField()[patchI].coupled())
                {
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

                    // Live inputs for this law on this patch. The patch
                    // values are what a boundary face sees, not the values
                    // in the cells behind it
                    const mechanicalConstitutiveLawInputs inputs
                    (
                        lawInputsPatch(lawI, patchI, faces, dt, tp)
                    );

                    // "View" into the kinematic and stress fields for this
                    // material => does not copy data
                    const UIndirectList<tensor> gradDView
                    (
                        gradD.boundaryField()[patchI], faces
                    );
                    const UIndirectList<tensor> gradD0View
                    (
                        gradD0.boundaryField()[patchI], faces
                    );
                    UIndirectList<symmTensor> stressView
                    (
                        Foam::boundaryFieldRef(stress)[patchI], faces
                    );

                    // Create wrapper for kinematic data: input to material law
                    // This does not copy data
                    smallStrainMechanicalConstitutiveLawKinematics kin
                    (
                        gradDView, gradD0View
                    );

                    // Create wrapper for output
                    if (scalarTangentPtr && needsScalarTangent(tangentReq))
                    {
                        // "View" into the tangent for this material => does not copy
                        // data
                        UIndirectList<scalar> tangentView
                        (
                            scalarTangentPtr->boundaryField()[patchI],
                            faces
                        );

                        // Create wrapper for material: output from material law
                        // This does not copy data
                        mechanicalConstitutiveLawResponse response
                        (
                            stressView, tangentView, tangentReq
                        );

                        // Update the material response, e.g. update the stress
                        laws_[lawI].evaluate
                        (
                            kin, inputs, tp.boundaryStates_[lawI][patchI], response
                        );
                    }
                    else
                    {
                        // Create wrapper for material: output from material law
                        // This does not copy data
                        mechanicalConstitutiveLawResponse response
                        (
                            stressView, tangentReq
                        );

                        // Update the material response, e.g. update the stress
                        laws_[lawI].evaluate
                        (
                            kin, inputs, tp.boundaryStates_[lawI][patchI], response
                        );
                    }
                }
            }
        }
    }

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
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
    else if (fourthOrderTangentPtr)
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

    // Look up the map and state for face-based topologies
    const integrationPointTopology& topo =
        topologyFor(faceCentredIntegrationPointTopology::typeName);

    topologyEntry& tp = topology(topo);

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Live inputs for this evaluation. Built once and passed through
    // every evaluation, including each finite-difference perturbation,
    // so there is no per-call forwarding to get wrong
    const mechanicalConstitutiveLawInputs inputs
    (
        inputsWithoutCoupling(dt, "updateStressSmallStrain")
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

        const UIndirectList<tensor> gradDView
        (
            gradD.internalField(), ipIDs
        );

        const UIndirectList<tensor> gradD0View
        (
            gradD0.internalField(), ipIDs
        );

        UIndirectList<symmTensor> stressView
        (
            stress.internalField(), ipIDs
        );

        smallStrainMechanicalConstitutiveLawKinematics kin
        (
            gradDView, gradD0View
        );

        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView
            (
                *scalarTangentPtr, ipIDs
            );

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, tp.states_[lawI], response);
        }
        else if
        (
            fourthOrderTangentPtr && needsFourthOrderTangent(tangentReq)
        )
        {
            UIndirectList<mat66> tangentView(*fourthOrderTangentPtr, ipIDs);

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, tp.states_[lawI], response);
        }
        else
        {
            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, tp.states_[lawI], response);
        }

        // Update stress accumulation fields used for stress collapse
        forAll(ipIDs, i)
        {
            const label faceI = ipIDs[i];

            stressSum[faceI] += stress[faceI];
            weightSum[faceI] += 1.0;

            if (scalarTangentPtr && needsScalarTangent(tangentReq))
            {
                const scalar K = (*scalarTangentPtr)[faceI];
                if (collapseRule == stressCollapseRule::average)
                {
                    (*tangentWeightPtr)[faceI] += K;
                }
                else if (collapseRule == stressCollapseRule::harmonic)
                {
                    (*tangentWeightPtr)[faceI] += 1.0/K;
                }
                else
                {
                    FatalErrorInFunction
                        << "Invalid stress collapse rule combination:\n"
                        << "collapseRule = "
                        << stressCollapseRuleName(collapseRule)
                        << nl
                        << "tangentReq   = " << tangentRequestName(tangentReq)
                        << nl
                        << "scalarTangentPtr = "
                        << (scalarTangentPtr ? "set" : "null")
                        << exit(FatalError);
                }
            }
        }

        // Optionally, update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        if (tp.boundaryAware_)
        {
            forAll(gradD.boundaryField(), patchI)
            {
                if (!gradD.boundaryField()[patchI].coupled())
                {
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

                    // "View" into the kinematic and stress fields for this
                    // material => does not copy data
                    const UIndirectList<tensor> gradDView
                    (
                        gradD.boundaryField()[patchI], faces
                    );
                    const UIndirectList<tensor> gradD0View
                    (
                        gradD0.boundaryField()[patchI], faces
                    );
                    UIndirectList<symmTensor> stressView
                    (
                        Foam::boundaryFieldRef(stress)[patchI], faces
                    );

                    // Create wrapper for kinematic data: input to material law
                    // This does not copy data
                    smallStrainMechanicalConstitutiveLawKinematics kin
                    (
                        gradDView, gradD0View
                    );

                    // Create wrapper for output
                    if (scalarTangentPtr && needsScalarTangent(tangentReq))
                    {
                        // "View" into the tangent for this material => does not copy
                        // data
                        UIndirectList<scalar> tangentView
                        (
                            scalarTangentPtr->boundaryField()[patchI],
                            faces
                        );

                        // Create wrapper for material: output from material law
                        // This does not copy data
                        mechanicalConstitutiveLawResponse response
                        (
                            stressView, tangentView, tangentReq
                        );

                        // Update the material response, e.g. update the stress
                        laws_[lawI].evaluate
                        (
                            kin, inputs, tp.boundaryStates_[lawI][patchI], response
                        );
                    }
                    else // Note: fourth order tangents not needed on the boundary
                    {
                        // Create wrapper for material: output from material law
                        // This does not copy data
                        mechanicalConstitutiveLawResponse response
                        (
                            stressView, tangentRequest::none
                        );

                        // Update the material response, e.g. update the stress
                        laws_[lawI].evaluate
                        (
                            kin, inputs, tp.boundaryStates_[lawI][patchI], response
                        );
                    }
                }
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

        // Stress collapse as arithmetic mean
        stress[faceI] = stressSum[faceI]/w;

        // Tangent collapse, if requested
        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            if (collapseRule == stressCollapseRule::average)
            {
                // Arithmetic mean of the contributing tangents
                (*scalarTangentPtr)[faceI] = (*tangentWeightPtr)[faceI]/w;
            }
            else if (collapseRule == stressCollapseRule::harmonic)
            {
                // Harmonicaly average the tangent
                (*scalarTangentPtr)[faceI] = w/(*tangentWeightPtr)[faceI];
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid stress collapse rule combination:\n"
                    << "collapseRule = " << stressCollapseRuleName(collapseRule)
                    << nl
                    << "tangentReq   = " << tangentRequestName(tangentReq) << nl
                    << "scalarTangentPtr = "
                    << (scalarTangentPtr ? "set" : "null")
                    << exit(FatalError);
            }
        }
    }

// This one guard is real, and only this one. OpenFOAM.org's fvsPatchField
// has no evaluate(), so correctBoundaryConditions() does not compile for a
// SURFACE field there. It compiles and is needed for volFields, where the
// guard was previously applied too and silently left the boundary values
// uncorrected on that fork
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

    // Live inputs for this evaluation, passed through unchanged
    const mechanicalConstitutiveLawInputs inputs
    (
        inputsWithoutCoupling(dt, "updateStressSmallStrain")
    );

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

        const UIndirectList<tensor> gradDView
        (
            gradD.internalField(), ipIDs
        );

        const UIndirectList<tensor> gradD0View
        (
            gradD0.internalField(), ipIDs
        );

        UIndirectList<symmTensor> stressView
        (
            stress.internalField(), ipIDs
        );

        smallStrainMechanicalConstitutiveLawKinematics kin
        (
            gradDView, gradD0View
        );

        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            UIndirectList<scalar> tangentView
            (
                scalarTangentPtr->internalField(), ipIDs
            );

            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, tp.states_[lawI], response);
        }
        else
        {
            mechanicalConstitutiveLawResponse response
            (
                stressView, tangentReq
            );

            laws_[lawI].evaluate(kin, inputs, tp.states_[lawI], response);
        }

        // Accumulate per point

        forAll(ipIDs, i)
        {
            const label pointI = ipIDs[i];

            stressSum[pointI] += stress[pointI];
            weightSum[pointI] += 1.0;

            if (tangentWeightPtr && scalarTangentPtr)
            {
                const scalar K = (*scalarTangentPtr)[pointI];

                if (collapseRule == stressCollapseRule::average)
                {
                    (*tangentWeightPtr)[pointI] += K;
                }
                else if (collapseRule == stressCollapseRule::harmonic)
                {
                    (*tangentWeightPtr)[pointI] += 1.0/max(K, SMALL);
                }
                else
                {
                    FatalErrorInFunction
                        << "Invalid stress collapse rule"
                        << exit(FatalError);
                }
            }
        }
    }

    // Collapse

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

        stress[pointI] = stressSum[pointI]/w;

        if (scalarTangentPtr && needsScalarTangent(tangentReq))
        {
            if (collapseRule == stressCollapseRule::average)
            {
                (*scalarTangentPtr)[pointI] =
                    (*tangentWeightPtr)[pointI]/w;
            }
            else if (collapseRule == stressCollapseRule::harmonic)
            {
                (*scalarTangentPtr)[pointI] =
                    w/max((*tangentWeightPtr)[pointI], SMALL);
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid stress collapse rule"
                    << exit(FatalError);
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
    const tangentRequest tangentReq
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

    // Update old time fields at the start of a new time step
    updateOldTimeIfNeeded();

    // Live inputs for this evaluation. Built once and passed through
    // every evaluation, including each finite-difference perturbation,
    // so there is no per-call forwarding to get wrong
    const mechanicalConstitutiveLawInputs inputs
    (
        inputsWithoutCoupling(dt, "updateStressFiniteStrain")
    );

    // Look up the map and state for cell-based topologies
    const integrationPointTopology& topo =
        topologyFor(cellCentredIntegrationPointTopology::typeName);

    // Update the internal field via the flat-list primitive: a cell-centred
    // topology has one integration point per cell, so the internal fields are
    // already in the flat form it expects
    updateStressFiniteStrain
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
        scalarTangentPtr ? &Foam::primitiveFieldRef(*scalarTangentPtr) : nullptr,
        nullptr,
        tangentReq
    );

    topologyEntry& tp = topology(topo);

    forAll(laws_, lawI)
    {
        // Update the boundary field
        // Boundary constitutive response uses independent state objects,
        // allowing history-dependent laws to operate correctly on boundary faces.
        forAll(F.boundaryField(), patchI)
        {
            if (!F.boundaryField()[patchI].coupled())
            {
                // Select all faces on the patch for which the adjacent cell is
                // in this material
                const labelList& faces = lawBoundaryFaces_[lawI][patchI];

                // "View" into the J for this material => does not copy data
                const UIndirectList<scalar> JView
                (
                    J.boundaryField()[patchI], faces
                );

                // "View" into the J0 for this material => does not copy data
                const UIndirectList<scalar> J0View
                (
                    J0.boundaryField()[patchI], faces
                );

                // "View" into the F for this material => does not copy data
                const UIndirectList<tensor> FView
                (
                    F.boundaryField()[patchI], faces
                );

                // "View" into the F0 for this material => does not copy data
                const UIndirectList<tensor> F0View
                (
                    F0.boundaryField()[patchI], faces
                );

                // "View" into the Finv for this material => does not copy data
                const UIndirectList<tensor> FinvView
                (
                    Finv.boundaryField()[patchI], faces
                );

                // "View" into the Finv0 for this material => does not copy data
                const UIndirectList<tensor> Finv0View
                (
                    Finv0.boundaryField()[patchI], faces
                );

                // "View" into the stress for this material => does not copy data
                UIndirectList<symmTensor> stressView
                (
                    Foam::boundaryFieldRef(stress)[patchI], faces
                );

                // Create wrapper for kinematic data: input to material law
                // This does not copy data
                finiteStrainMechanicalConstitutiveLawKinematics kin
                (
                    FView, F0View, JView, J0View, FinvView, Finv0View
                );

                // Create wrapper for output
                if (scalarTangentPtr && needsScalarTangent(tangentReq))
                {
                    // "View" into the tangent for this material => does not copy
                    // data
                    UIndirectList<scalar> tangentView
                    (
                        scalarTangentPtr->boundaryField()[patchI],
                        faces
                    );

                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentView, tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, inputs, tp.boundaryStates_[lawI][patchI], response
                    );
                }
                else
                {
                    // Create wrapper for material: output from material law
                    // This does not copy data
                    mechanicalConstitutiveLawResponse response
                    (
                        stressView, tangentReq
                    );

                    // Update the material response, e.g. update the stress
                    laws_[lawI].evaluate
                    (
                        kin, inputs, tp.boundaryStates_[lawI][patchI], response
                    );
                }
            }
        }
    }

    // Update boundaries including syncing coupled boundaries
    stress.correctBoundaryConditions();

    if (scalarTangentPtr && needsScalarTangent(tangentReq))
    {
        scalarTangentPtr->correctBoundaryConditions();
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
    checkCompactLayoutConsistency
    (
        F,
        F0,
        stress,
        scalarTangentPtr,
        "updateStressFiniteStrain (CompactListList)"
    );

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


void Foam::mechanicalConstitutiveLawManager::endTimeStep()
{
    const scalar time = mesh_.time().value();
    const label timeIndex = mesh_.time().timeIndex();

    // Loop over all topology entries
    forAllIter
    (
        HashTable<autoPtr<topologyEntry>>, topologyEntries_, topoIter
    )
    {
        topologyEntry& tp = topoIter()();

        forAll(laws_, lawI)
        {
            // Call endTimeStep with each topo type for each law
            laws_[lawI].endTimeStep(tp.states_[lawI], time, timeIndex);
        }
    }
}


// ************************************************************************* //
