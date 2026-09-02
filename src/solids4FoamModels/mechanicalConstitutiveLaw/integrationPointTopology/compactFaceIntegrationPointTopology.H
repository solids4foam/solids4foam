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

public:

    //- Runtime type information
    TypeName("compactFaceIntegrationPointTopology");

    //- Construct from mesh and face-to-integration-point addressing
    compactFaceIntegrationPointTopology
    (
        const fvMesh& mesh,
        CompactListList<label>&& faceToIntegrationPointIDs
    )
    :
        integrationPointTopology(mesh),
        integrationPointIDs_(),
        faceStartOffsets_(),
        cellIntegrationPointIDs_(),
        cellStartOffsets_()
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

        // Gather each cell's points from its faces. A cell's faces are its
        // own, so this is the natural association; an internal face
        // contributes to both of its cells, which is what makes the points
        // shared
        const cellList& cells = mesh_.cells();

        labelList cellSizes(mesh_.nCells(), 0);
        forAll(cells, cellI)
        {
            const labelList& cFaces = cells[cellI];
            forAll(cFaces, i)
            {
                cellSizes[cellI] += rowSizes[cFaces[i]];
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

    //- Does this topology store boundary integration points per patch
    virtual bool boundaryAware() const override
    {
        return false;
    }

    //- A face is shared by two cells, so the same quadrature point is reached
    //  from both and must not be aggregated twice
    virtual bool requiresUniqueIntegrationPointsPerMaterial() const override
    {
        return true;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
