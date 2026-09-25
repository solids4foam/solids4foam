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

#include "mechanicalConstitutiveLawStateSetup.H"
#include "mechanicalConstitutiveLawStateIO.H"
#include "labelIOList.H"
#include "mechanicalConstitutiveLawGather.H"
#include "calculatedFvPatchFields.H"
#include "compatibilityFunctions.H"

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
            part.internal->getField<Type>(plan.variable);
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
                part.patchStates[patchI]->getField<Type>(plan.variable);
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


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::mechanicalConstitutiveLawStateSetup::prescribedField
(
    const word& name
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    // The file comes first, and the registry only after it.
    //
    // Reading the file first means what the case says on disk is what
    // the law gets, even if something has registered a field of the
    // same name, and the registry is left to serve the fields that
    // only ever exist in memory because another model computes them
    IOobject io
    (
        name,
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!mechanicalConstitutiveLawHeaderIsA<VolFieldType>(io))
    {
        // A prescribed field describes the case, not the state at a
        // particular time, so it is written once into 0 and a restart
        // will not find it beside the fields it restarts from. Look
        // there before giving up
        io.instance() = "0";
    }

    if (mechanicalConstitutiveLawHeaderIsA<VolFieldType>(io))
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


template<class Type>
void Foam::mechanicalConstitutiveLawStateSetup::readPrescribed
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

    gatherToIntegrationPoints
    (
        mesh_, lawCells_[lawI], src, topo, ipIDs, fld
    );
}


template<class Type>
void Foam::mechanicalConstitutiveLawStateSetup::readPrescribedPatch
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalConstitutiveLawStateSetup::mechanicalConstitutiveLawStateSetup
(
    const fvMesh& mesh,
    const PtrList<mechanicalConstitutiveLaw>& laws,
    const wordList& lawNames,
    const List<labelList>& lawCells,
    const List<List<labelList>>& lawBoundaryFaces,
    PtrList<regIOobject>& stateWriteProxies,
    const bool& restartKinematicsAvailable
)
:
    mesh_(mesh),
    laws_(laws),
    lawNames_(lawNames),
    lawCells_(lawCells),
    lawBoundaryFaces_(lawBoundaryFaces),
    stateWriteProxies_(stateWriteProxies),
    restartKinematicsAvailable_(restartKinematicsAvailable)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mechanicalConstitutiveLawStateSetup::readPrescribedFields
(
    const mechanicalConstitutiveLawStateSpec& spec,
    mechanicalConstitutiveLawState& state,
    const prescribedSource source,
    const label lawI,
    const integrationPointTopology* topoPtr,
    const labelList* ipIDsPtr,
    const label patchI
) const
{
    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        const mechanicalConstitutiveLawStateSpec::entry& e = es[i];

        if (e.role != mechanicalConstitutiveLawStateSpec::stateRole::prescribed)
        {
            continue;
        }

        // A prescribed field is read from a field the case supplied. Its
        // absence is not an error: the declared default stands, which is what
        // lets a law declare one without changing any existing case.
        //
        // The old time is set to match, because a prescribed field is never
        // written and a law reads it at old time so that a tangent query
        // evaluated into a shadow state sees it
        if (e.typeName == "scalar")
        {
            Field<scalar>& f = state.scalarField(e.name);

            if (source == prescribedSource::internalPoints)
            {
                readPrescribed<scalar>
                (
                    e.name, lawI, *topoPtr, *ipIDsPtr, f
                );
            }
            else
            {
                readPrescribedPatch<scalar>(e.name, lawI, patchI, f);
            }

            state.scalarField0(e.name) = f;
        }
        else if (e.typeName == "vector")
        {
            Field<vector>& f = state.vectorField(e.name);

            if (source == prescribedSource::internalPoints)
            {
                readPrescribed<vector>
                (
                    e.name, lawI, *topoPtr, *ipIDsPtr, f
                );
            }
            else
            {
                readPrescribedPatch<vector>(e.name, lawI, patchI, f);
            }

            state.vectorField0(e.name) = f;
        }
        else if (e.typeName == "tensor")
        {
            Field<tensor>& f = state.tensorField(e.name);

            if (source == prescribedSource::internalPoints)
            {
                readPrescribed<tensor>
                (
                    e.name, lawI, *topoPtr, *ipIDsPtr, f
                );
            }
            else
            {
                readPrescribedPatch<tensor>(e.name, lawI, patchI, f);
            }

            state.tensorField0(e.name) = f;
        }
        else if (e.typeName == "symmTensor")
        {
            Field<symmTensor>& f = state.symmTensorField(e.name);

            if (source == prescribedSource::internalPoints)
            {
                readPrescribed<symmTensor>
                (
                    e.name, lawI, *topoPtr, *ipIDsPtr, f
                );
            }
            else
            {
                readPrescribedPatch<symmTensor>(e.name, lawI, patchI, f);
            }

            state.symmTensorField0(e.name) = f;
        }
        else
        {
            // As at the other places a state type is dispatched on. No law
            // can reach this, since only the typed adders set a type name,
            // but a type added to one place and not here must not be read as
            // a symmTensor
            FatalErrorInFunction
                << "State field '" << e.name << "' has unsupported type '"
                << e.typeName << "'." << exit(FatalError);
        }
    }
}


void Foam::mechanicalConstitutiveLawStateSetup::applyStateDefaults
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


bool Foam::mechanicalConstitutiveLawStateSetup::declaresPersistentState
(
    const mechanicalConstitutiveLaw& law
) const
{
    mechanicalConstitutiveLawStateSpec spec;
    law.declareState(spec);

    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        if
        (
            es[i].role
         == mechanicalConstitutiveLawStateSpec::stateRole::persistent
        )
        {
            return true;
        }
    }

    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        if (declaresPersistentState(law.childLaw(childNames[i])))
        {
            return true;
        }
    }

    return false;
}


void Foam::mechanicalConstitutiveLawStateSetup::checkRestartKinematics() const
{
    bool anyPersistent = false;

    forAll(laws_, lawI)
    {
        if (declaresPersistentState(laws_[lawI]))
        {
            anyPersistent = true;
            break;
        }
    }

    if (!anyPersistent)
    {
        return;
    }

    // History is measured against the displacement gradient of the step it
    // began in, so restoring one without the other is half a restart. An
    // incremental law then takes the strain accumulated since the beginning of
    // the run as though it happened in one step, which is not a small error
    // and does not announce itself: the run continues and the answer is wrong.
    //
    // Whether the kinematic history was written is the solid model's to say,
    // since it is the solid model that writes it
    const bool present = restartKinematicsAvailable_;

    if (!present)
    {
        FatalErrorInFunction
            << "Restarting from time " << mesh_.time().timeName()
            << " with a material that keeps history, but the displacement "
            << "gradient at old time was never written." << nl << nl
            << "Set" << nl << nl
            << "    restart yes;" << nl << nl
            << "in the solidModel's coefficients dictionary and run the case "
            << "again from the start. It makes the solid model write the "
            << "kinematic history a constitutive history is measured against."
            << nl << nl
            << "Without it the constitutive state would come back and the "
            << "strain it is measured against would not, and an incremental "
            << "law would read the whole run's strain as one step's."
            << exit(FatalError);
    }
}


void Foam::mechanicalConstitutiveLawStateSetup::setupStateRestart
(
    const List<labelList>& lawIntegrationPointIDs,
    PtrList<mechanicalConstitutiveLawState>& states,
    PtrList<PtrList<mechanicalConstitutiveLawState>>& boundaryStates,
    const bool boundaryAware,
    const integrationPointTopology& topo,
    const word& topologyKey
) const
{
    // A run that begins at a time other than the first is continuing, and a
    // history dependent law must pick its history back up. A run beginning at
    // the start has none to pick up and writes from its first output.
    //
    // Asked of the time the run *started* at, not the time it has reached. A
    // topology is created the first time something asks for it, which can be
    // well after the run has advanced, and judging by the current index would
    // call an ordinary run a restart the moment it took a step and then refuse
    // to find files it never wrote
    const bool isRestart = mesh_.time().startTimeIndex() > 0;

    if (isRestart)
    {
        checkRestartKinematics();
    }

    forAll(laws_, lawI)
    {
        // The pieces of this law's state, in the order they are written: the
        // law's own integration points, then its points on each patch. A
        // topology with no boundary points contributes the first alone
        List<mechanicalConstitutiveLawState*> parts(1);
        parts[0] = &states[lawI];

        if (boundaryAware)
        {
            parts.setSize(1 + boundaryStates[lawI].size());

            forAll(boundaryStates[lawI], patchI)
            {
                parts[1 + patchI] = &boundaryStates[lawI][patchI];
            }
        }

        // Where each of this law's points sits on the mesh, in the order the
        // state is written. Written alongside the state so that a run on a
        // different decomposition can match entry to entry by mesh entity
        // instead of by position - the state's own order comes from a hash set
        // for the boundary part and means nothing across decompositions.
        //
        // Only a topology whose points are the cells for now: a cell maps
        // through decomposePar's addressing directly, whereas a point- or
        // dual-based one needs its own translation that is not written yet.
        // The others write no locations, and a restart on a changed
        // decomposition then refuses rather than guessing
        labelList entities;

        const bool topologyRecordsLocations =
            topo.integrationPointsAreCells();

        if (topologyRecordsLocations)
        {
            const labelList& ipIDs = lawIntegrationPointIDs[lawI];

            label nEntities = ipIDs.size();

            if (boundaryAware)
            {
                forAll(mesh_.boundary(), patchI)
                {
                    nEntities += lawBoundaryFaces_[lawI][patchI].size();
                }
            }

            entities.setSize(nEntities);

            label k = 0;

            forAll(ipIDs, i)
            {
                entities[k++] =
                    mechanicalConstitutiveLawStateIO::cellEntity(ipIDs[i]);
            }

            if (boundaryAware)
            {
                forAll(mesh_.boundary(), patchI)
                {
                    const labelList& pf = lawBoundaryFaces_[lawI][patchI];
                    const label start = mesh_.boundaryMesh()[patchI].start();

                    forAll(pf, i)
                    {
                        entities[k++] =
                            mechanicalConstitutiveLawStateIO::faceEntity
                            (
                                start + pf[i]
                            );
                    }
                }
            }
        }

        const word entityName
        (
            mechanicalConstitutiveLawStateIO::fieldName
            (
                lawNames_[lawI], topologyKey, wordList(), "integrationPoints"
            )
        );

        // On whether this topology records locations, not on whether this
        // rank has any. A rank holding no cells of this material has an empty
        // list, and that is a statement about the decomposition rather than an
        // absence of information: the file still has to be written, because
        // reconstructing a decomposed state requires the pair of files from
        // every processor directory and refuses the lot if one is missing.
        // A topology that records no locations at all writes neither file, and
        // a restart on a changed decomposition then refuses rather than
        // guessing, which is the intended behaviour
        if (topologyRecordsLocations)
        {
            // Written in the numbering of the undecomposed mesh, whichever way
            // this run is being run. A serial run's locations are what a
            // decomposed one needs to distribute the state; a decomposed run's
            // are what a reconstructed one needs to put it back together. Said
            // in the same numbering, one list serves both
            const labelList written
            (
                mechanicalConstitutiveLawStateIO::toSerialEntities
                (
                    mesh_, entities, entityName
                )
            );

            const label proxyI = stateWriteProxies_.size();
            stateWriteProxies_.setSize(proxyI + 1);

            stateWriteProxies_.set
            (
                proxyI,
                new labelIOList
                (
                    IOobject
                    (
                        entityName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    written
                )
            );

            stateWriteProxies_[proxyI].note() =
                mechanicalConstitutiveLawStateIO::decompositionIdentity(mesh_);
        }

        setupStateRestartLaw
        (
            laws_[lawI],
            lawNames_[lawI],
            topologyKey,
            wordList(),
            entityName,
            entities,
            parts,
            isRestart
        );
    }
}


void Foam::mechanicalConstitutiveLawStateSetup::setupStateRestartLaw
(
    const mechanicalConstitutiveLaw& law,
    const word& lawName,
    const word& topologyName,
    const wordList& childPath,
    const word& entityName,
    const labelUList& entities,
    const List<mechanicalConstitutiveLawState*>& parts,
    const bool isRestart
) const
{
    mechanicalConstitutiveLawStateSpec spec;
    law.declareState(spec);

    const UList<mechanicalConstitutiveLawStateSpec::entry>& es = spec.entries();

    forAll(es, i)
    {
        const mechanicalConstitutiveLawStateSpec::entry& e = es[i];

        // Only history is written. A prescribed field is the user's input and
        // is read again from where it came; a fixed one cannot change
        if (e.role != mechanicalConstitutiveLawStateSpec::stateRole::persistent)
        {
            continue;
        }

        const word name
        (
            mechanicalConstitutiveLawStateIO::fieldName
            (
                lawName, topologyName, childPath, e.name
            )
        );

        if (e.typeName == "scalar")
        {
            mechanicalConstitutiveLawStateIO::restartField<scalar>
            (
                mesh_, name, entityName, entities, parts, e.name, isRestart,
                stateWriteProxies_
            );
        }
        else if (e.typeName == "vector")
        {
            mechanicalConstitutiveLawStateIO::restartField<vector>
            (
                mesh_, name, entityName, entities, parts, e.name, isRestart,
                stateWriteProxies_
            );
        }
        else if (e.typeName == "tensor")
        {
            mechanicalConstitutiveLawStateIO::restartField<tensor>
            (
                mesh_, name, entityName, entities, parts, e.name, isRestart,
                stateWriteProxies_
            );
        }
        else if (e.typeName == "symmTensor")
        {
            mechanicalConstitutiveLawStateIO::restartField<symmTensor>
            (
                mesh_, name, entityName, entities, parts, e.name, isRestart,
                stateWriteProxies_
            );
        }
        else
        {
            FatalErrorInFunction
                << "State field '" << e.name << "' has unsupported type '"
                << e.typeName << "'." << exit(FatalError);
        }
    }

    // A composite keeps its sub-laws' history in child states. Without this
    // the composite would restart and its sub-laws would not, which is the
    // failure that is hardest to see: the run continues and only the part of
    // the answer the sub-law owns is wrong
    const wordList childNames(law.childStateNames());

    forAll(childNames, i)
    {
        wordList path(childPath.size() + 1);

        forAll(childPath, j)
        {
            path[j] = childPath[j];
        }

        path[childPath.size()] = childNames[i];

        List<mechanicalConstitutiveLawState*> childParts(parts.size());

        forAll(parts, partI)
        {
            childParts[partI] = &parts[partI]->child(childNames[i]);
        }

        setupStateRestartLaw
        (
            law.childLaw(childNames[i]),
            lawName,
            topologyName,
            path,
            entityName,
            entities,
            childParts,
            isRestart
        );
    }
}


void Foam::mechanicalConstitutiveLawStateSetup::applyStateSpec
(
    const label lawI,
    const integrationPointTopology& topo,
    const labelList& ipIDs,
    mechanicalConstitutiveLawState& state
) const
{
    applyStateSpec(laws_[lawI], lawI, topo, ipIDs, state);
}


void Foam::mechanicalConstitutiveLawStateSetup::applyStateSpec
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

    readPrescribedFields
    (
        spec,
        state,
        prescribedSource::internalPoints,
        lawI,
        &topo,
        &ipIDs,
        -1
    );
}


void Foam::mechanicalConstitutiveLawStateSetup::applyStateSpecScratch
(
    const label lawI,
    mechanicalConstitutiveLawState& state
) const
{
    applyStateSpecScratch(laws_[lawI], state);
}


void Foam::mechanicalConstitutiveLawStateSetup::applyStateSpecScratch
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


void Foam::mechanicalConstitutiveLawStateSetup::applyStateSpecPatch
(
    const label lawI,
    const label patchI,
    mechanicalConstitutiveLawState& state
) const
{
    applyStateSpecPatch(laws_[lawI], lawI, patchI, state);
}


void Foam::mechanicalConstitutiveLawStateSetup::applyStateSpecPatch
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

    readPrescribedFields
    (
        spec,
        state,
        prescribedSource::patchFaces,
        lawI,
        nullptr,
        nullptr,
        patchI
    );
}


void Foam::mechanicalConstitutiveLawStateSetup::writeStateFields
(
    const List<labelList>& lawIntegrationPointIDs,
    const PtrList<mechanicalConstitutiveLawState>& states,
    const PtrList<PtrList<mechanicalConstitutiveLawState>>& boundaryStates,
    const bool boundaryAware,
    HashSet<word>& reportedStateFields
) const
{
    HashTable<stateFieldPlan> plans;

    forAll(laws_, lawI)
    {
        if (lawI >= states.size() || !states.set(lawI))
        {
            continue;
        }

        stateFieldPart part;
        part.internal = &states[lawI];
        part.cells = &lawIntegrationPointIDs[lawI];
        part.patchStates.setSize(mesh_.boundary().size(), nullptr);
        part.patchFaces.setSize(mesh_.boundary().size(), nullptr);

        if
        (
            boundaryAware
         && lawI < boundaryStates.size()
         && lawI < lawBoundaryFaces_.size()
        )
        {
            const PtrList<mechanicalConstitutiveLawState>& bStates =
                boundaryStates[lawI];

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
            if (!reportedStateFields.found("state:" + fieldName))
            {
                reportedStateFields.insert("state:" + fieldName);

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


// ************************************************************************* //
