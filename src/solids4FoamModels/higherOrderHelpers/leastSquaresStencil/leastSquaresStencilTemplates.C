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

#include "leastSquaresStencil.H"
#include "fvMesh.H"

// * * * * * * * * * * *  Public Member Functions  * * * * * * * * * * * * * //


template<class Type>
Foam::Map<Type>
Foam::leastSquaresStencil::remoteFieldMap
(
    const UList<Type>& fld
) const
{
    const List<labelList>& remoteCells = remoteCellsPerProc();
    const List<Field<Type>> remoteFld = remoteFieldPerProc(fld);

    label nRemote = 0;
    forAll(remoteCells, procI)
    {
        nRemote += remoteCells[procI].size();
    }

    Map<Type> remoteFldMap(2*nRemote);

    // Fill map
    forAll(remoteCells, procI)
    {
        const labelList& globalIDs = remoteCells[procI];
        const Field<Type>& values = remoteFld[procI];

        if (globalIDs.size() != values.size())
        {
            FatalErrorInFunction
                << "Mismatch between remote cell IDs and received field values "
                << "for proc " << procI
                << abort(FatalError);
        }

        forAll(globalIDs, i)
        {
            remoteFldMap.insert(globalIDs[i], values[i]);
        }
    }

    return remoteFldMap;
}


template<class Type>
Foam::List<Foam::Field<Type>>
Foam::leastSquaresStencil::remoteFieldPerProc
(
    const UList<Type>& fld
) const
{
    const List<labelList>& remoteCells = remoteCellsPerProc();

    List<Field<Type>> remoteFld(Pstream::nProcs());

    if (!Pstream::parRun())
    {
        return remoteFld;
    }

    PstreamBuffers reqBufs(Pstream::commsTypes::nonBlocking);

    // Phase 1: send requested global cell IDs to each processor
    forAll(remoteCells, procI)
    {
        if (procI == Pstream::myProcNo())
        {
            continue;
        }

        if (!remoteCells[procI].empty())
        {
            UOPstream os(procI, reqBufs);
            os << remoteCells[procI];
        }
    }

    reqBufs.finishedSends();

    // Phase 2: reply with field values for requested cells
    PstreamBuffers repBufs(Pstream::commsTypes::nonBlocking);

    forAll(remoteCells, procI)
    {
        if (procI == Pstream::myProcNo())
        {
            continue;
        }

#ifdef OPENFOAM_COM
        if (!reqBufs.recvDataCount(procI))
        {
            continue;
        }
#endif

        UIPstream is(procI, reqBufs);

        labelList requestedGlobalIDs;
        is >> requestedGlobalIDs;

        Field<Type> values(requestedGlobalIDs.size());

        forAll(requestedGlobalIDs, i)
        {
            const label globalCellID = requestedGlobalIDs[i];
            const label localCellID = globalCells_.toLocal(globalCellID);

            if (localCellID < 0 || localCellID >= mesh_.nCells())
            {
                FatalErrorInFunction
                    << "Invalid global->local mapping: globalCellID="
                    << globalCellID << ", localCellID=" << localCellID
                    << ", on proc " << Pstream::myProcNo()
                    << abort(FatalError);
            }

            values[i] = fld[localCellID];
        }

        UOPstream os(procI, repBufs);
        os << values;
    }

    repBufs.finishedSends();

    // Phase 3: receive field values
    forAll(remoteCells, procI)
    {
        if (procI == Pstream::myProcNo())
        {
            continue;
        }

        if (remoteCells[procI].empty())
        {
            continue;
        }

#ifdef OPENFOAM_COM
        if (!repBufs.recvDataCount(procI))
        {
            FatalErrorInFunction
                << "Did not receive field values from proc " << procI
                << abort(FatalError);
        }
#endif

        UIPstream is(procI, repBufs);

        Field<Type> values;
        is >> values;

        if (values.size() != remoteCells[procI].size())
        {
            FatalErrorInFunction
                << "Field reply size mismatch from proc " << procI
                << ": got " << values.size()
                << ", expected " << remoteCells[procI].size()
                << abort(FatalError);
        }

        remoteFld[procI].transfer(values);
    }

    return remoteFld;
}


// ************************************************************************* //
