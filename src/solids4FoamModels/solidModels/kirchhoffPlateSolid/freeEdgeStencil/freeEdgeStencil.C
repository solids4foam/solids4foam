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

#ifndef OPENFOAM_ORG

#include "freeEdgeStencil.H"
#include "faMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelListList Foam::freeEdgeStencil::edgeNeighbours(const faPatch& p)
{
    const labelListList& pointEdges = p.pointEdges();

    List<DynamicList<label> > dynNbrs(p.size());

    forAll(pointEdges, pointI)
    {
        const labelList& pe = pointEdges[pointI];

        if (pe.size() == 2)
        {
            dynNbrs[pe[0]].append(pe[1]);
            dynNbrs[pe[1]].append(pe[0]);
        }
    }

    labelListList nbrs(p.size());

    forAll(nbrs, edgeI)
    {
        nbrs[edgeI].transfer(dynNbrs[edgeI]);
    }

    return nbrs;
}


Foam::labelListList Foam::freeEdgeStencil::edgeEndPoints(const faPatch& p)
{
    const labelList& pointLabels = p.pointLabels();
    const labelListList& pointEdges = p.pointEdges();

    List<DynamicList<label> > dynEnds(p.size());

    forAll(pointEdges, pointI)
    {
        const labelList& pe = pointEdges[pointI];

        if (pe.size() == 1)
        {
            dynEnds[pe[0]].append(pointLabels[pointI]);
        }
    }

    labelListList ends(p.size());

    forAll(ends, edgeI)
    {
        ends[edgeI].transfer(dynEnds[edgeI]);
    }

    return ends;
}


bool Foam::freeEdgeStencil::fixesValueAtPoint
(
    const areaScalarField& fld,
    const faPatch& p,
    const label pointI
)
{
    const faBoundaryMesh& boundary = p.boundaryMesh();

    forAll(boundary, patchI)
    {
        if (patchI == p.index() || !fld.boundaryField()[patchI].fixesValue())
        {
            continue;
        }

        const labelList& pointLabels = boundary[patchI].pointLabels();

        forAll(pointLabels, pI)
        {
            if (pointLabels[pI] == pointI)
            {
                return true;
            }
        }
    }

    return false;
}


Foam::vectorField Foam::freeEdgeStencil::edgeTangents
(
    const faPatch& p,
    const labelListList& nbrs
)
{
    const vectorField& Ce = p.edgeCentres();

    vectorField t(p.size(), vector::zero);
    boolList visited(p.size(), false);

    forAll(nbrs, startI)
    {
        if (visited[startI])
        {
            continue;
        }

        // Walk back to the start of the chain holding this edge; a closed
        // loop is walked from this edge
        label prev = -1;
        label cur = startI;

        while (nbrs[cur].size() == 2)
        {
            const label next =
                (nbrs[cur][0] == prev) ? nbrs[cur][1] : nbrs[cur][0];

            if (next == startI)
            {
                break;
            }

            prev = cur;
            cur = next;
        }

        // The edge behind the first edge of a closed loop is its other
        // neighbour
        prev = (nbrs[cur].size() == 2) ? nbrs[cur][1] : -1;

        // Walk forward along the chain assigning the tangents
        while (true)
        {
            visited[cur] = true;

            const labelList& nbr = nbrs[cur];
            label next = -1;

            if (nbr.size() == 2)
            {
                next = (nbr[0] == prev) ? nbr[1] : nbr[0];
            }
            else if (nbr.size() == 1 && nbr[0] != prev)
            {
                next = nbr[0];
            }

            const vector from = (prev == -1) ? Ce[cur] : Ce[prev];
            const vector to = (next == -1) ? Ce[cur] : Ce[next];

            t[cur] = (to - from)/(mag(to - from) + VSMALL);

            if (next == -1 || visited[next])
            {
                break;
            }

            prev = cur;
            cur = next;
        }
    }

    return t;
}

#endif // OPENFOAM_ORG

// ************************************************************************* //
