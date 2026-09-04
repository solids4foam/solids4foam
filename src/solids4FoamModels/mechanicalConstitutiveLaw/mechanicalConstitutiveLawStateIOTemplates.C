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
    const fileName path
    (
        mechanicalConstitutiveLawStateIO::globalCasePath(mesh)
      / mesh.time().timeName()
      / name
    );

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
    const label n = totalSize(parts);

    if (isRestart)
    {
        IOobject readIO
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        // A run continuing from a time directory has history to restore. If it
        // is not there then the state would silently fall back to its cold
        // start defaults and the run would continue a different calculation
        // that looks plausible, which is the failure this whole path exists to
        // remove. Refuse instead
        // Both OpenFOAM.com and OpenFOAM.org spell this typeHeaderOk;
        // foam-extend has only the untyped headerOk
#ifdef OPENFOAM_NOT_EXTEND
        bool present = readIO.typeHeaderOk<IOField<Type>>(false);
#else
        bool present = readIO.headerOk();
#endif

        // Set when the values came from the undecomposed case rather than
        // from this processor's own directory
        bool mapped = false;

        // Every rank takes the same route. They should agree already - the
        // state is written by all of them at once - but if a directory were
        // half deleted they would not, and then some would go on to read the
        // processor addressing while others did not. Under a collated file
        // handler that read can be collective, and a collective some ranks
        // skip is a hang rather than an error
        if (Pstream::parRun())
        {
            bool allPresent = present;
            reduce(allPresent, andOp<bool>());
            present = allPresent;
        }

        if (!present && Pstream::parRun())
        {
            // Nothing in this processor's directory. decomposePar does not
            // copy these files, so a case decomposed after it was run has its
            // state sitting in the undecomposed directory - which is readable,
            // and which decomposePar's own addressing says how to distribute
            Field<Type> serialData;
            labelList serialEntities;
            string dataNote, entityNote;

            const bool haveData =
                readSerialList(mesh, name, serialData, dataNote);
            const bool haveEntities =
                readSerialList(mesh, entityName, serialEntities, entityNote);

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

                const label serialNCells =
                    serialCellCount(dataNote);

                if (serialNCells < 0)
                {
                    FatalErrorInFunction
                        << "The state file '" << name << "' in the "
                        << "undecomposed case was not written by a serial run."
                        << nl << "  its header says: " << dataNote
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

                const labelList positions
                (
                    serialPositions
                    (
                        mesh, entities, serialEntities, serialNCells, name
                    )
                );

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

                present = true;
                mapped = true;
            }
        }

        if (!present)
        {
            FatalErrorInFunction
                << "Restarting from time " << mesh.time().timeName()
                << " but the constitutive state field '" << name
                << "' is not there." << nl
                << "A history dependent material cannot continue without it."
                << nl
                << "It is written from the first time step of a run that uses "
                << "the mechanical constitutive law framework. A case "
                << "decomposed after it was run keeps its state in the "
                << "undecomposed directory, and that is looked for too, so "
                << "this means neither is present."
                << exit(FatalError);
        }

        if (!mapped)
        {
        const IOField<Type> stored(readIO);

        // The header says which decomposition wrote this. A file of the right
        // length from the wrong one reads as plausible values in the wrong
        // cells, which is the single way this could still be silently wrong
        const string expected
        (
            decompositionIdentity(mesh)
        );

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

        if (stored.size() != n)
        {
            FatalErrorInFunction
                << "Constitutive state field '" << name << "' holds "
                << stored.size() << " values but this run needs " << n << '.'
                << exit(FatalError);
        }

        // Scatter back in write order, and set the old time to match. The old
        // time is deliberately not written: at the instant a time directory is
        // written the two are equal, so copying is what reading a second file
        // would have produced anyway
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


// ************************************************************************* //
