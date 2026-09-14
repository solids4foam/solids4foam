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
    Integration-point topology for higher-order solid mechanics formulations
    using compact per-cell storage.

    This class defines a mapping from mesh cells to integration-point indices
    when integration-point data are stored in a compact form, such as a
    CompactListList<Type>. Each cell is associated with one or more
    integration-point indices, corresponding to entries in the packed
    constitutive working set.

    The integration-point indices returned by this topology refer to positions
    in the flat storage returned by CompactListList<Type>::values(). The
    topology itself is independent of material definitions and constitutive
    behaviour, and only describes how integration points are associated with
    mesh cells.

    This class is intended for use with higher-order solid models that evaluate
    constitutive laws at multiple integration points per cell.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef compactCellIntegrationPointTopology_H
#define compactCellIntegrationPointTopology_H

#include "integrationPointTopology.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class compactCellIntegrationPointTopology
\*---------------------------------------------------------------------------*/

class compactCellIntegrationPointTopology
:
    public integrationPointTopology
{
    //- Integration-point indices for every cell, concatenated.
    //  The indices refer to positions in the packed constitutive working set
    //  (e.g. CompactListList<Type>::values())
    labelList integrationPointIDs_;

    //- Start of each cell's entries in integrationPointIDs_, of size
    //  nCells + 1.
    //  Held flat rather than as a CompactListList because that container's
    //  const operator[] does not compile on foam-extend, and because the two
    //  forks disagree on whether an offset marks the start or the end of a row
    labelList startOffsets_;

public:

    //- Runtime type information
    TypeName("compactCellIntegrationPointTopology");

    //- Construct from mesh and cell-to-integration-point addressing
    compactCellIntegrationPointTopology
    (
        const fvMesh& mesh,
        CompactListList<label>&& cellToIntegrationPointIDs
    )
    :
        integrationPointTopology(mesh),
        integrationPointIDs_(),
        startOffsets_()
    {
        // Basic consistency check: one entry per mesh cell
        if (cellToIntegrationPointIDs.size() != mesh_.nCells())
        {
            FatalErrorInFunction
                << "Size of cellToIntegrationPointIDs ("
                << cellToIntegrationPointIDs.size()
                << ") does not match number of mesh cells ("
                << mesh_.nCells() << ")."
                << exit(FatalError);
        }

        // Flatten. sizes() and m() mean the same thing on every supported
        // fork, which const operator[] does not
        const labelList rowSizes(cellToIntegrationPointIDs.sizes());

        startOffsets_.setSize(rowSizes.size() + 1);
        startOffsets_[0] = 0;
        forAll(rowSizes, cellI)
        {
            startOffsets_[cellI + 1] = startOffsets_[cellI] + rowSizes[cellI];
        }

        integrationPointIDs_ = cellToIntegrationPointIDs.m();
    }

    //- Destructor
    virtual ~compactCellIntegrationPointTopology() = default;

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
            const_cast<label*>(integrationPointIDs_.begin())
          + startOffsets_[cellI],
            startOffsets_[cellI + 1] - startOffsets_[cellI]
        );
    }

    //- Does this topology store boundary integration points per patch
    virtual bool boundaryAware() const override
    {
        return false;
    }

    //- Should integration-point indices be unique when aggregating across
    //  multiple cells belonging to the same material.
    virtual bool requiresUniqueIntegrationPointsPerMaterial() const override
    {
        return false;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
