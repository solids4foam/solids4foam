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
        os << "serial nCells=" << mesh.nCells();

        return os.str();
    }

    // Which cells this rank holds, not merely how many. Two decompositions
    // into the same number of processors can agree on the count and disagree
    // on everything else, and it is the disagreement that has to be caught
    const labelList addr(cellProcAddressing(mesh));

    label hash = 0;

    forAll(addr, i)
    {
        // Order matters as much as membership: the state is written in cell
        // order, so a permutation is as wrong as a different set
        hash = (hash*31 + addr[i] + i) & 0x7FFFFFFF;
    }

    os  << "procs=" << Pstream::nProcs()
        << " rank=" << Pstream::myProcNo()
        << " nCells=" << mesh.nCells()
        << " addr=" << hash;

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

    return readLabel(IStringStream(s.substr(key.size()))());
}


// ************************************************************************* //
