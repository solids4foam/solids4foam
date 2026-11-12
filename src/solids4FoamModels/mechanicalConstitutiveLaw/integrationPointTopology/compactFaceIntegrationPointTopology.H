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

Description

Class
    Foam::compactFaceIntegrationPointTopology

Description
    Integration points packed several per mesh FACE, as produced by a
    high-order discretisation that evaluates the stress at face quadrature
    points.

    This is the face-based counterpart of compactCellIntegrationPointTopology,
    and the two differ in more than which entity indexes the rows. A face is
    shared by two cells, so a face quadrature point is reachable from both and
    the same point is aggregated twice when the manager gathers the points of
    a material. That is why requiresUniqueIntegrationPointsPerMaterial() is
    true here and false there, and it is what restricts a flat-list update
    over this topology to a single material.

    With more than one material the points of a face are instead given to
    one of its cells only (cellAssociation): the owner, or the neighbour.
    Each point then belongs to exactly one material, and a face between two
    materials is evaluated once from each side, by two topologies: the
    owner-side one, holding every face, and a neighbour-side one, holding
    the material-interface faces only.

SourceFiles
    None: header only

\*---------------------------------------------------------------------------*/

#ifndef compactFaceIntegrationPointTopology_H
#define compactFaceIntegrationPointTopology_H

#include "integrationPointTopology.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class compactFaceIntegrationPointTopology Declaration
\*---------------------------------------------------------------------------*/

class compactFaceIntegrationPointTopology
:
    public integrationPointTopology
{
public:

    //- Which cells a face's points are associated with
    enum class cellAssociation
    {
        both,       // the owner and the neighbour: points are shared
        owner,      // the owner only
        neighbour   // the neighbour only (internal faces)
    };


private:

    //- Which cells a face's points are associated with
    cellAssociation association_;

    //- Flattened integration-point indices of every face, in face order
    labelList integrationPointIDs_;

    //- Start of each face's block in integrationPointIDs_.
    //  Stored as explicit starts rather than a CompactListList, since the
    //  forks disagree on whether an offset marks the start or the end of a row
    labelList faceStartOffsets_;

    //- Flattened integration-point indices of every cell, gathered from the
    //  faces of that cell
    labelList cellIntegrationPointIDs_;

    //- Start of each cell's block in cellIntegrationPointIDs_
    labelList cellStartOffsets_;

    //- The face each integration point lies on
    labelList integrationPointFaces_;

public:

    //- Runtime type information
    TypeName("compactFaceIntegrationPointTopology");

    //- Construct from mesh and face-to-integration-point addressing.
    //  With cellAssociation::neighbour every boundary face must have no
    //  points, as it has no neighbour
    compactFaceIntegrationPointTopology
    (
        const fvMesh& mesh,
        CompactListList<label>&& faceToIntegrationPointIDs,
        const cellAssociation association = cellAssociation::both
    )
    :
        integrationPointTopology(mesh),
        association_(association),
        integrationPointIDs_(),
        faceStartOffsets_(),
        cellIntegrationPointIDs_(),
        cellStartOffsets_(),
        integrationPointFaces_()
    {
        if (faceToIntegrationPointIDs.size() != mesh_.nFaces())
        {
            FatalErrorInFunction
                << "Size of faceToIntegrationPointIDs ("
                << faceToIntegrationPointIDs.size()
                << ") does not match the number of mesh faces ("
                << mesh_.nFaces() << ")." << nl
                << "A layout with one row per cell is a "
                << "compactCellIntegrationPointTopology instead."
                << exit(FatalError);
        }

        // Flatten the face rows
        const labelList rowSizes(faceToIntegrationPointIDs.sizes());

        faceStartOffsets_.setSize(rowSizes.size() + 1);
        faceStartOffsets_[0] = 0;
        forAll(rowSizes, faceI)
        {
            faceStartOffsets_[faceI + 1] =
                faceStartOffsets_[faceI] + rowSizes[faceI];
        }

        integrationPointIDs_ = faceToIntegrationPointIDs.m();

        // Each point's face, which is the row it appears in
        integrationPointFaces_.setSize(integrationPointIDs_.size(), -1);
        forAll(rowSizes, faceI)
        {
            for
            (
                label k = faceStartOffsets_[faceI];
                k < faceStartOffsets_[faceI + 1];
                ++k
            )
            {
                integrationPointFaces_[integrationPointIDs_[k]] = faceI;
            }
        }

        // Gather each cell's points from its faces. A cell's faces are its
        // own, so this is the natural association; an internal face
        // contributes to both of its cells, which is what makes the points
        // shared. With an owner or neighbour association it contributes to
        // that cell only
        const cellList& cells = mesh_.cells();
        const labelUList& own = mesh_.faceOwner();
        const labelUList& nei = mesh_.faceNeighbour();
        const label nInternalFaces = mesh_.nInternalFaces();

        if (association_ == cellAssociation::neighbour)
        {
            for (label faceI = nInternalFaces; faceI < mesh_.nFaces(); ++faceI)
            {
                if (rowSizes[faceI] > 0)
                {
                    FatalErrorInFunction
                        << "Boundary face " << faceI << " has integration "
                        << "points, but the points are associated with the "
                        << "face neighbour, which a boundary face does not "
                        << "have" << exit(FatalError);
                }
            }
        }

        // Whether cellI takes the points of its face faceI
        const auto takesFace = [&](const label cellI, const label faceI)
        {
            if (association_ == cellAssociation::owner)
            {
                return own[faceI] == cellI;
            }
            else if (association_ == cellAssociation::neighbour)
            {
                return faceI < nInternalFaces && nei[faceI] == cellI;
            }

            return true;
        };

        labelList cellSizes(mesh_.nCells(), 0);
        forAll(cells, cellI)
        {
            const labelList& cFaces = cells[cellI];
            forAll(cFaces, i)
            {
                if (takesFace(cellI, cFaces[i]))
                {
                    cellSizes[cellI] += rowSizes[cFaces[i]];
                }
            }
        }

        cellStartOffsets_.setSize(mesh_.nCells() + 1);
        cellStartOffsets_[0] = 0;
        forAll(cellSizes, cellI)
        {
            cellStartOffsets_[cellI + 1] =
                cellStartOffsets_[cellI] + cellSizes[cellI];
        }

        cellIntegrationPointIDs_.setSize(cellStartOffsets_[mesh_.nCells()]);

        label pos = 0;
        forAll(cells, cellI)
        {
            const labelList& cFaces = cells[cellI];
            forAll(cFaces, i)
            {
                const label faceI = cFaces[i];

                if (!takesFace(cellI, faceI))
                {
                    continue;
                }

                for
                (
                    label k = faceStartOffsets_[faceI];
                    k < faceStartOffsets_[faceI + 1];
                    ++k
                )
                {
                    cellIntegrationPointIDs_[pos++] = integrationPointIDs_[k];
                }
            }
        }
    }

    //- Destructor
    virtual ~compactFaceIntegrationPointTopology() = default;

    //- Total number of integration points
    virtual label nIntegrationPoints() const override
    {
        return integrationPointIDs_.size();
    }

    //- Integration-point indices associated with a given cell
    virtual const labelUList cellIntegrationPointIDs
    (
        const label cellI
    ) const override
    {
        return labelUList
        (
            const_cast<label*>(cellIntegrationPointIDs_.begin())
          + cellStartOffsets_[cellI],
            cellStartOffsets_[cellI + 1] - cellStartOffsets_[cellI]
        );
    }

    //- Every point lies on a face, the one whose row holds it
    virtual label integrationPointFace(const label ipI) const override
    {
        return integrationPointFaces_[ipI];
    }

    //- Does this topology store boundary integration points per patch
    virtual bool boundaryAware() const override
    {
        return false;
    }

    //- Which cells a face's points are associated with
    cellAssociation association() const
    {
        return association_;
    }

    //- With both cells associated, a face is shared by two cells, so the
    //  same quadrature point is reached from both and must not be aggregated
    //  twice. With one cell, every point belongs to that cell alone
    virtual bool requiresUniqueIntegrationPointsPerMaterial() const override
    {
        return association_ == cellAssociation::both;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
