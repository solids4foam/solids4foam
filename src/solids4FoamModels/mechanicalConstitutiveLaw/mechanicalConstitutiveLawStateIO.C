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
#include "labelIOList.H"
#include "polyMesh.H"
#include "Pstream.H"
#include "OStringStream.H"
#include "HashTable.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::mechanicalConstitutiveLawStateIO::fieldName
(
    const word& lawName,
    const word& topologyName,
    const wordList& childPath,
    const word& variableName
)
{
    // Every part is needed to make the name unique. The law, because two
    // materials may each declare epsilonP and mean different things by it. The
    // topology, because the same law may hold state at cell centres and at
    // faces in the same run, and the two are different lengths. The child
    // path, because a composite law and its sub-law each have their own
    // unqualified namespace and may both use the same name
    word name(lawName + ':' + topologyName);

    forAll(childPath, i)
    {
        name += ':' + childPath[i];
    }

    return name + ':' + variableName;
}


Foam::label Foam::mechanicalConstitutiveLawStateIO::totalSize
(
    const stateParts& parts
)
{
    label n = 0;

    forAll(parts, partI)
    {
        n += parts[partI]->size();
    }

    return n;
}

Foam::labelList Foam::mechanicalConstitutiveLawStateIO::cellProcAddressing
(
    const fvMesh& mesh
)
{
    if (!Pstream::parRun())
    {
        return labelList();
    }

    // Written by decomposePar into each processor's mesh directory, against
    // the instance the mesh faces came from rather than constant, which is
    // what a moving or topology-changing mesh needs
    IOobject io
    (
        "cellProcAddressing",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

#ifdef OPENFOAM_NOT_EXTEND
    const bool present = io.typeHeaderOk<labelIOList>(false);
#else
    const bool present = io.headerOk();
#endif

    if (!present)
    {
        return labelList();
    }

    return labelList(labelIOList(io));
}


Foam::string Foam::mechanicalConstitutiveLawStateIO::decompositionIdentity
(
    const fvMesh& mesh
)
{
    OStringStream os;

    if (!Pstream::parRun())
    {
        // Cells, faces and points rather than cells alone. A cell count is a
        // weak thing to identify a mesh by - two unrelated meshes can share
        // one - and the entity numbering this file is written in means
        // nothing except against the mesh that produced it
        os  << "serial nCells=" << mesh.nCells()
            << " nFaces=" << mesh.nFaces()
            << " nPoints=" << mesh.nPoints();

        return os.str();
    }

    // Which cells this rank holds, not merely how many. Two decompositions
    // into the same number of processors can agree on the count and disagree
    // on everything else, and it is the disagreement that has to be caught
    const labelList addr(cellProcAddressing(mesh));

    // Order matters as much as membership: the state is written in cell
    // order, so a permutation is as wrong as a different set.
    //
    // Two hashes rather than one, mixed differently and reported side by side.
    // This is the only check standing between a stale processor directory and
    // a silently wrong restart, so the cost of a collision is a wrong answer
    // that looks right, and one 31-bit checksum is a thin thing to hang that
    // on
    unsigned long hashA = 5381;
    unsigned long hashB = 0;

    forAll(addr, i)
    {
        const unsigned long v = (unsigned long)(addr[i]);

        hashA = ((hashA << 5) + hashA) + v;
        hashB = hashB*1000003UL + (v ^ (unsigned long)(i));
    }

    const label hash = label(hashA & 0x7FFFFFFF);
    const label hash2 = label((hashB >> 7) & 0x7FFFFFFF);

    os  << "procs=" << Pstream::nProcs()
        << " rank=" << Pstream::myProcNo()
        << " nCells=" << mesh.nCells()
        << " addr=" << hash << ',' << hash2;

    return os.str();
}


Foam::labelList Foam::mechanicalConstitutiveLawStateIO::faceProcAddressing
(
    const fvMesh& mesh
)
{
    if (!Pstream::parRun())
    {
        return labelList();
    }

    IOobject io
    (
        "faceProcAddressing",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

#ifdef OPENFOAM_NOT_EXTEND
    const bool present = io.typeHeaderOk<labelIOList>(false);
#else
    const bool present = io.headerOk();
#endif

    if (!present)
    {
        return labelList();
    }

    return labelList(labelIOList(io));
}


Foam::fileName Foam::mechanicalConstitutiveLawStateIO::globalCasePath
(
    const fvMesh& mesh
)
{
    return mesh.time().rootPath()/mesh.time().globalCaseName();
}


Foam::label Foam::mechanicalConstitutiveLawStateIO::serialCellCount
(
    const string& note
)
{
    // The note reads "serial nCells=<n>". Anything else is not a serial file
    const std::string s(note);
    const std::string key("serial nCells=");

    if (s.compare(0, key.size(), key) != 0)
    {
        return -1;
    }

    // The count runs to the next space; the rest of the note is the other
    // half of the mesh's fingerprint and is compared as a whole elsewhere
    return readLabel(IStringStream(s.substr(key.size()))());
}


//- Where each of this processor's state entries sat in the serial run.
//  The serial file lists, for every entry it holds, the piece of mesh that
//  entry belongs to. This processor knows the same thing about its own
//  entries, in its own numbering, and decomposePar's addressing translates
//  between the two - so the two lists are matched on mesh entity rather than
//  on position, and the order either was written in stops mattering
Foam::labelList Foam::mechanicalConstitutiveLawStateIO::serialPositions
(
    const fvMesh& mesh,
    const labelUList& localEntities,
    const labelUList& serialEntities,
    const label serialNCells,
    const word& name
)
{
    const labelList cellAddr
    (
        cellProcAddressing(mesh)
    );
    const labelList faceAddr
    (
        faceProcAddressing(mesh)
    );

    if
    (
        !cellAddr.empty()
     && !faceAddr.empty()
     && (cellAddr.size() != mesh.nCells() || faceAddr.size() != mesh.nFaces())
    )
    {
        FatalErrorInFunction
            << "Restarting '" << name << "': the processor addressing does "
            << "not describe this mesh." << nl
            << "  cellProcAddressing holds " << cellAddr.size()
            << " entries for " << mesh.nCells() << " cells" << nl
            << "  faceProcAddressing holds " << faceAddr.size()
            << " entries for " << mesh.nFaces() << " faces" << nl
            << "This mesh was not the one decomposePar wrote that addressing "
            << "for."
            << exit(FatalError);
    }

    if (cellAddr.empty() || faceAddr.empty())
    {
        FatalErrorInFunction
            << "Restarting '" << name << "' from the undecomposed case needs "
            << "cellProcAddressing and faceProcAddressing, which decomposePar "
            << "writes into each processor's mesh directory." << nl
            << "They are not there, so this decomposition cannot be related "
            << "to the serial one."
            << exit(FatalError);
    }

    // Serial entity to its position in the serial file
    HashTable<label, label, Hash<label>> positionOf(2*serialEntities.size());

    forAll(serialEntities, i)
    {
        positionOf.insert(serialEntities[i], i);
    }

    labelList positions(localEntities.size(), -1);

    forAll(localEntities, i)
    {
        const label local = localEntities[i];

        // Translate this processor's entity into the serial one
        label serial = -1;

        if (local < mesh.nCells())
        {
            serial =
                cellEntity(cellAddr[local]);
        }
        else
        {
            const label localFace = local - mesh.nCells();

            // Offset by one and signed, so that face zero can carry a sign
            const label serialFace = mag(faceAddr[localFace]) - 1;

            serial = faceEntity
            (
                serialNCells, serialFace
            );
        }

        HashTable<label, label, Hash<label>>::const_iterator iter =
            positionOf.find(serial);

        // Only a face decomposition itself created may fall back. The
        // fallback is right for those and wrong for every other face: on a
        // real boundary a miss means the serial state and this mesh disagree
        // about something, and quietly substituting a nearby value would hide
        // exactly the mistake this mapping is supposed to refuse
        bool decompositionFace = false;

        if (local >= mesh.nCells())
        {
            const label patchI =
                mesh.boundaryMesh().whichPatch(local - mesh.nCells());

            decompositionFace =
                patchI >= 0 && mesh.boundary()[patchI].coupled();
        }

        if (iter == positionOf.end() && decompositionFace)
        {
            // A face this processor calls a boundary that the serial mesh
            // called interior: decomposition made it one, by cutting the mesh
            // there. The serial run kept no history on it, because a
            // cell-centred topology keeps none on interior faces, so there is
            // nothing to look up and nothing has been lost.
            //
            // Its owner cell is the nearest thing that does have history, and
            // is what the face's own state grew from, so that is what it gets.
            // A decomposed run continued from a serial one is therefore not
            // identical to a decomposed run continued from a decomposed one at
            // these faces, and cannot be: the history it would need was never
            // written because it never existed
            const label localFace = local - mesh.nCells();
            const label ownerCell = mesh.faceOwner()[localFace];

            iter = positionOf.find
            (
                cellEntity
                (
                    cellAddr[ownerCell]
                )
            );
        }

        if (iter == positionOf.end())
        {
            FatalErrorInFunction
                << "Restarting '" << name << "': this processor holds a "
                << "point on mesh entity " << serial << " of the serial mesh, "
                << "and the serial state does not have one there." << nl
                << "The serial state and this mesh do not describe the same "
                << "case, or the material a cell belongs to has changed."
                << exit(FatalError);
        }

        positions[i] = iter();
    }

    return positions;
}


// ************************************************************************* //
