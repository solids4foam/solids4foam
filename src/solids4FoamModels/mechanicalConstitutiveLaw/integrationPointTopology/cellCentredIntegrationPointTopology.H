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
    Cell-centred integration-point topology for low-order solid mechanics
    formulations.

    This class provides a trivial integration-point topology in which each
    mesh cell is associated with exactly one integration point, identified by
    the cell index itself. The constitutive working set therefore coincides
    with the mesh cell set, and integration-point indices are identical to
    cell indices.

    This topology is suitable for standard cell-centred finite-volume solid
    models, where constitutive quantities (e.g. stress, strain) are evaluated
    once per cell.

    The class contains no material information and does not perform any
    constitutive evaluation; it only describes the mapping between mesh cells
    and integration-point indices.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef cellCentredIntegrationPointTopology_H
#define cellCentredIntegrationPointTopology_H

#include "integrationPointTopology.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class cellCentredIntegrationPointTopology
\*---------------------------------------------------------------------------*/

class cellCentredIntegrationPointTopology
:
    public integrationPointTopology
{
    //- Per-cell integration-point addressing.
    //  Each cell is associated with a single integration point, whose
    //  index is equal to the cell index.
    List<labelList> cellToIntegrationPointIDs_;

public:

    //- Runtime type information
    TypeName("cellCentredIntegrationPointTopology");

    //- Construct from mesh
    explicit cellCentredIntegrationPointTopology(const fvMesh& mesh)
    :
        integrationPointTopology(mesh),
        cellToIntegrationPointIDs_(mesh.nCells())
    {
        forAll(cellToIntegrationPointIDs_, cellI)
        {
            cellToIntegrationPointIDs_[cellI].setSize(1);
            cellToIntegrationPointIDs_[cellI][0] = cellI;
        }
    }

    //- Destructor
    virtual ~cellCentredIntegrationPointTopology() = default;

    //- Total number of integration points
    virtual label nIntegrationPoints() const override
    {
        return mesh_.nCells();
    }

    //- Integration-point indices associated with a given cell
    virtual const labelUList cellIntegrationPointIDs
    (
        const label cellI
    ) const override
    {
        return cellToIntegrationPointIDs_[cellI];
    }

    //- Does this topology store boundary integration points per patch
    virtual bool boundaryAware() const override
    {
        return true;
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
