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

#include "mechanicalConstitutiveLawStateIO.H"
#include "IFstream.H"
#include "Pstream.H"
#include "labelIOList.H"
#include "OSspecific.H"
#include "fileName.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::stateIOFieldProxy<Type>::writeData(Ostream& os) const
{
    // Gather the parts into one list in write order. The state is the
    // authoritative copy and is read here rather than mirrored, so there is no
    // second copy that could be stale by the time the registry writes
    Field<Type> all(mechanicalConstitutiveLawStateIO::totalSize(parts_));

    label n = 0;
    forAll(parts_, partI)
    {
        const Field<Type>& f =
            stateFieldAccess<Type>::get(*parts_[partI], variableName_);

        forAll(f, i)
        {
            all[n++] = f[i];
        }
    }

    // Written at full round-trip precision, whatever the case asks of its
    // other output. writePrecision is a choice about reading numbers and
    // looking at them; nobody reads this file, and at the default of 6 the
    // rounding alone moves a continued run by 6e-6 - a thousand times the
    // 1e-8 the restart otherwise achieves. 17 significant digits is what an
    // IEEE double survives a decimal round trip with
    const int oldPrecision = os.precision(17);

    os << all;

    os.precision(oldPrecision);

    return os.good();
}


//- Read a list written by a serial run, from the undecomposed case.
//  decomposePar does not copy these files into the processor directories, so
//  a decomposed restart goes back to the serial one and maps it. Read by hand
//  rather than through the registry because the object wanted is in another
//  case directory entirely
template<class ListType>
bool Foam::mechanicalConstitutiveLawStateIO::readSerialList
(
    const fvMesh& mesh,
    const word& name,
    ListType& data,
    string& note
)
{
    return readListFrom
    (
        globalCasePath(mesh)/mesh.time().timeName(),
        mesh,
        name,
        data,
        note
    );
}


template<class ListType>
bool Foam::mechanicalConstitutiveLawStateIO::readListFrom
(
    const fileName& dir,
    const fvMesh& mesh,
    const word& name,
    ListType& data,
    string& note
)
{
    const fileName path(dir/name);

    IFstream is(path);

    if (!is.good())
    {
        return false;
    }

    // Not registered: this is a file from another case directory being read
    // for its contents, and it must not appear in this run's registry
    IOobject io
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!io.readHeader(is))
    {
        return false;
    }

    is >> data;
    note = io.note();

    return is.good() || is.eof();
}


//- Register one state variable for writing, and read it back on a restart.
//
//  The state is written as it is held, a flat list, so that the value read is
//  the value written. A geometric field would have been viewable and would
//  have survived a change of decomposition, and would also have needed a
//  dimensionSet no law declares and a boundary field whose evaluation would
//  overwrite the history being restored. Visualisation is a separate,
//  write-only projection; this path is for restarting, and it is exact
template<class Type>
bool Foam::mechanicalConstitutiveLawStateIO::readOwnFile
(
    const fvMesh& mesh,
    const word& name,
    const stateParts& parts,
    const word& variableName
)
{
    IOobject readIO
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    // Both OpenFOAM.com and OpenFOAM.org spell this typeHeaderOk; foam-extend
    // has only the untyped headerOk
#ifdef OPENFOAM_NOT_EXTEND
    bool present = readIO.typeHeaderOk<IOField<Type>>(false);
#else
    bool present = readIO.headerOk();
#endif

    // Every rank takes the same route. They should agree already - the state
    // is written by all of them at once - but if a directory were half deleted
    // they would not, and then some would go on to read the processor
    // addressing while others did not. Under a collated file handler that read
    // can be collective, and a collective some ranks skip is a hang rather
    // than an error
    if (Pstream::parRun())
    {
        bool allPresent = present;
        reduce(allPresent, andOp<bool>());
        present = allPresent;
    }

    if (!present)
    {
        return false;
    }

    const IOField<Type> stored(readIO);

    // The header says which decomposition wrote this. A file of the right
    // length from the wrong one reads as plausible values in the wrong cells,
    // which is the single way this could still be silently wrong
    const string expected(decompositionIdentity(mesh));

    if (stored.note() != expected)
    {
        FatalErrorInFunction
            << "Constitutive state field '" << name << "' was written for "
            << "a different decomposition." << nl
            << "  the file says: " << stored.note() << nl
            << "  this run is:   " << expected << nl
            << "This usually means processor directories were left behind "
            << "by an earlier decomposePar. The values would be the right "
            << "count and the wrong cells, so they are refused rather "
            << "than read."
            << exit(FatalError);
    }

    if (stored.size() != totalSize(parts))
    {
        FatalErrorInFunction
            << "Constitutive state field '" << name << "' holds "
            << stored.size() << " values but this run needs "
            << totalSize(parts) << '.'
            << exit(FatalError);
    }

    // Scatter back in write order, and set the old time to match. The old time
    // is deliberately not written: at the instant a time directory is written
    // the two are equal, so copying is what reading a second file would have
    // produced anyway
    label k = 0;

    forAll(parts, partI)
    {
        Field<Type>& f =
            stateFieldAccess<Type>::ref(*parts[partI], variableName);
        Field<Type>& f0 =
            stateFieldAccess<Type>::ref0(*parts[partI], variableName);

        forAll(f, i)
        {
            f[i] = stored[k];
            f0[i] = stored[k];
            k++;
        }
    }

    return true;
}


template<class Type>
bool Foam::mechanicalConstitutiveLawStateIO::distributeFromSerial
(
    const fvMesh& mesh,
    const word& name,
    const word& entityName,
    const labelUList& entities,
    const stateParts& parts,
    const word& variableName
)
{
    // Nothing in this processor's directory. decomposePar does not
    // copy these files, so a case decomposed after it was run has its
    // state sitting in the undecomposed directory - which is readable,
    // and which decomposePar's own addressing says how to distribute
    Field<Type> serialData;
    labelList serialEntities;
    string dataNote, entityNote;

    bool haveData =
        readSerialList(mesh, name, serialData, dataNote);
    bool haveEntities =
        readSerialList(mesh, entityName, serialEntities, entityNote);

    // Agreed across the ranks before any of them acts on it. These are
    // direct reads that do not go through the file handler, so a rank
    // can fail one on its own; the next step reads the processor
    // addressing through the registry, which the handler may make
    // collective, and a collective some ranks skip is a hang
    bool haveBoth = haveData && haveEntities;
    reduce(haveBoth, andOp<bool>());
    haveData = haveBoth;
    haveEntities = haveBoth;

    if (haveData && haveEntities)
    {
        // The values and the locations they belong to are two files,
        // and they only mean anything as a pair. One of them left over
        // from an earlier run would map real values through the wrong
        // locations, which is a silent wrong answer rather than a
        // missing one, so they have to agree on which run wrote them
        if (dataNote != entityNote)
        {
            FatalErrorInFunction
                << "The serial state '" << name << "' and the "
                << "integration point locations it is read with were "
                << "written by different runs." << nl
                << "  the state says:     " << dataNote << nl
                << "  the locations say:  " << entityNote << nl
                << "They are written together and are only meaningful "
                << "together."
                << exit(FatalError);
        }

        if (serialData.size() != serialEntities.size())
        {
            FatalErrorInFunction
                << "The serial state '" << name << "' holds "
                << serialData.size() << " values and "
                << serialEntities.size() << " integration point "
                << "locations. They are written together and must "
                << "agree." << exit(FatalError);
        }

        // The pair is self-consistent; that does not make it this
        // case. decomposePar preserves cells, so the undecomposed cell
        // count the file records has to be the number of cells this
        // run holds between them - a collective, but one every rank
        // reaches, having just agreed to be here
        const label serialNCells = serialCellCount(dataNote);

        if (serialNCells < 0)
        {
            FatalErrorInFunction
                << "The state file '" << name << "' in the "
                << "undecomposed case was not written by a serial run."
                << nl << "  its header says: " << dataNote
                << exit(FatalError);
        }

        const label totalCells =
            returnReduce(mesh.nCells(), sumOp<label>());

        if (serialNCells != totalCells)
        {
            FatalErrorInFunction
                << "The state in the undecomposed case was written for "
                << "a mesh of " << serialNCells << " cells, and this "
                << "run holds " << totalCells << " between its ranks."
                << nl
                << "It is not the same mesh, so its integration point "
                << "locations mean nothing here."
                << exit(FatalError);
        }

        if (entities.empty())
        {
            FatalErrorInFunction
                << "Restarting '" << name << "' from the undecomposed case, "
                << "but this topology records no integration point locations, "
                << "so its entries cannot be matched to the ones in that file."
                << nl
                << "Only the cell-centred topology writes locations, so a case "
                << "on any other cannot change its decomposition between runs."
                << exit(FatalError);
        }

        const labelList positions
        (
            serialPositions(mesh, entities, serialEntities, name)
        );

        if (positions.size() != totalSize(parts))
        {
            FatalErrorInFunction
                << "Restarting '" << name << "': " << positions.size()
                << " locations were matched for " << totalSize(parts)
                << " values." << nl
                << "These have to agree; a mismatch means the locations do not "
                << "describe this run's state."
                << exit(FatalError);
        }

        label k = 0;
        forAll(parts, partI)
        {
            Field<Type>& f =
                stateFieldAccess<Type>::ref(*parts[partI], variableName);
            Field<Type>& f0 =
                stateFieldAccess<Type>::ref0(*parts[partI], variableName);

            forAll(f, i)
            {
                f[i] = serialData[positions[k]];
                f0[i] = serialData[positions[k]];
                k++;
            }
        }

        Info<< "    Mapped '" << name << "' from the undecomposed case"
            << endl;

        return true;
    }

    return false;
}




template<class Type>
void Foam::mechanicalConstitutiveLawStateIO::restartField
(
    const fvMesh& mesh,
    const word& name,
    const word& entityName,
    const labelUList& entities,
    const mechanicalConstitutiveLawStateIO::stateParts& parts,
    const word& variableName,
    const bool isRestart,
    PtrList<regIOobject>& proxies
)
{
    if (isRestart)
    {
        // Three places the state can be, in the order they are worth looking.
        // Its own directory, if this run has the shape the last one had; the
        // undecomposed case, if this run is decomposed and that one was not;
        // the processor directories, if the reverse. Each returns whether it
        // found it, and none of them guesses
        bool restored = readOwnFile<Type>(mesh, name, parts, variableName);

        if (!restored)
        {
            if (Pstream::parRun())
            {
                restored = distributeFromSerial<Type>
                (
                    mesh, name, entityName, entities, parts, variableName
                );
            }
            else
            {
                restored = gatherFromProcessors<Type>
                (
                    mesh, name, entityName, entities, parts, variableName
                );
            }
        }

        if (!restored)
        {
            FatalErrorInFunction
                << "Restarting from time " << mesh.time().timeName()
                << " but the constitutive state field '" << name
                << "' is not there." << nl
                << "A history dependent material cannot continue without it."
                << nl
                << "It is written from the first time step of a run that uses "
                << "the mechanical constitutive law framework. A case "
                << "decomposed or reconstructed after it was run keeps its "
                << "state where the run that produced it left it, and both "
                << "are looked for, so this means none of them is present."
                << exit(FatalError);
        }
    }

    // Registered so that the state is written whenever the run writes, by
    // whatever triggered it. The proxy gathers from the state at write time
    // rather than holding a copy that could be stale
    // set() rather than append(): foam-extend's PtrList has no append
    const label proxyI = proxies.size();
    proxies.setSize(proxyI + 1);

    proxies.set
    (
        proxyI,
        new stateIOFieldProxy<Type>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            parts,
            variableName,
            mesh
        )
    );
}


template<class Type>
bool Foam::mechanicalConstitutiveLawStateIO::gatherFromProcessors
(
    const fvMesh& mesh,
    const word& name,
    const word& entityName,
    const labelUList& serialEntities,
    const stateParts& parts,
    const word& variableName
)
{
    // A case run in parallel and then reconstructed. reconstructPar knows no
    // more about these files than decomposePar did, so the state is still
    // sitting in the processor directories, in as many pieces as there were
    // processors. Each piece says where its points were in the undecomposed
    // mesh, so the pieces can be put back together without knowing how the
    // mesh was cut
    const fileName casePath(globalCasePath(mesh));
    const word timeName(mesh.time().timeName());

    // What a value at each undecomposed mesh entity is
    HashTable<Type, label, Hash<label>> valueOf(4*serialEntities.size());

    label nProcs = 0;
    label totalProcCells = 0;

    while (true)
    {
        const fileName procDir
        (
            casePath/("processor" + Foam::name(nProcs))/timeName
        );

        if (!isDir(casePath/("processor" + Foam::name(nProcs))))
        {
            break;
        }

        Field<Type> procData;
        labelList procEntities;
        string dataNote, entityNote;

        const bool haveData =
            readListFrom(procDir, mesh, name, procData, dataNote);
        const bool haveEntities =
            readListFrom(procDir, mesh, entityName, procEntities, entityNote);

        if (!haveData || !haveEntities)
        {
            if (nProcs == 0)
            {
                return false;
            }

            FatalErrorInFunction
                << "Restarting '" << name << "' from the processor "
                << "directories, but processor" << nProcs << " has no state "
                << "at time " << timeName << '.' << nl
                << "The pieces of a decomposed state are only meaningful "
                << "together, so a missing one is refused rather than filled "
                << "in."
                << exit(FatalError);
        }

        if (dataNote != entityNote)
        {
            FatalErrorInFunction
                << "In processor" << nProcs << ", the state '" << name
                << "' and the integration point locations it is read with "
                << "were written by different runs." << nl
                << "  the state says:     " << dataNote << nl
                << "  the locations say:  " << entityNote
                << exit(FatalError);
        }

        if (procData.size() != procEntities.size())
        {
            FatalErrorInFunction
                << "In processor" << nProcs << ", the state '" << name
                << "' holds " << procData.size() << " values and "
                << procEntities.size() << " locations. They are written "
                << "together and must agree."
                << exit(FatalError);
        }

        // Each piece says how many cells it held. decomposePar preserves
        // cells, so they have to add up to this mesh, and if they do not then
        // these directories belong to some other case
        const label procCells = cellCountFromNote(dataNote);
        totalProcCells += (procCells > 0 ? procCells : 0);

        forAll(procEntities, i)
        {
            // set rather than insert: a face that decomposition cut is a
            // boundary on two processors and interior in the undecomposed
            // mesh, so both offer a value for it and neither is ever asked
            // for one
            valueOf.set(procEntities[i], procData[i]);
        }

        nProcs++;
    }

    if (nProcs > 0 && totalProcCells != mesh.nCells())
    {
        FatalErrorInFunction
            << "The state in " << nProcs << " processor directories was "
            << "written for a mesh of " << totalProcCells << " cells, and "
            << "this one has " << mesh.nCells() << '.' << nl
            << "Those directories belong to a different case, or to a "
            << "different mesh of the same case."
            << exit(FatalError);
    }

    if (nProcs == 0)
    {
        // Collated output puts every rank's data in one processors<N>
        // directory instead of a processorN each, and this reads neither the
        // layout nor the format. Saying nothing here would look exactly like a
        // case that was never run in parallel
        // OpenFOAM.org spells the file type as a scoped enumeration where
        // the others use a member constant
#ifdef OPENFOAM_ORG
        fileNameList dirs(readDir(casePath, fileType::directory));
#else
        fileNameList dirs(readDir(casePath, fileName::DIRECTORY));
#endif

        forAll(dirs, i)
        {
            if (dirs[i].find("processors") == 0)
            {
                FatalErrorInFunction
                    << "The constitutive state was written with collated file "
                    << "handling, into " << dirs[i] << '.' << nl
                    << "This reads the state from a processorN directory per "
                    << "rank, which is what the uncollated handler writes." << nl
                    << "Re-run the parallel leg with"
                    << " -fileHandler uncollated, or reconstruct it with a "
                    << "handler that writes the state where this can find it."
                    << exit(FatalError);
            }
        }

        return false;
    }

    if (serialEntities.size() != totalSize(parts))
    {
        FatalErrorInFunction
            << "Restarting '" << name << "' from the processor directories, "
            << "but this run records " << serialEntities.size()
            << " integration point locations for " << totalSize(parts)
            << " values." << nl
            << "Only the cell-centred topology records locations, so a case on "
            << "any other cannot be reassembled from a decomposed run."
            << exit(FatalError);
    }

    label k = 0;

    forAll(parts, partI)
    {
        Field<Type>& f =
            stateFieldAccess<Type>::ref(*parts[partI], variableName);
        Field<Type>& f0 =
            stateFieldAccess<Type>::ref0(*parts[partI], variableName);

        forAll(f, i)
        {
            typename HashTable<Type, label, Hash<label>>::const_iterator iter =
                valueOf.find(serialEntities[k]);

            if (iter == valueOf.end())
            {
                FatalErrorInFunction
                    << "Restarting '" << name << "' from " << nProcs
                    << " processor directories: this mesh holds a point on "
                    << "entity " << serialEntities[k] << " that none of them "
                    << "wrote." << nl
                    << "The processor directories and this mesh do not "
                    << "describe the same case."
                    << exit(FatalError);
            }

            f[i] = iter();
            f0[i] = iter();
            k++;
        }
    }

    Info<< "    Mapped '" << name << "' from " << nProcs
        << " processor directories" << endl;

    return true;
}


// ************************************************************************* //
