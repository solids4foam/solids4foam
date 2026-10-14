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
    Integration-point topology for point-centred (vertex-centred solid
    mechanics formulations.

    This topology defines one integration point per mesh point (vertex).
    Each cell maps to the integration points associated with its points.

    Mesh points are shared by two or more cells and therefore appear in the
    integration-point lists of both adjacent cells. As a result,
    integration-point indices must be deduplicated when aggregating
    per-material working sets.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pointCentredIntegrationPointTopology_H
#define pointCentredIntegrationPointTopology_H

#include "integrationPointTopology.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class pointCentredIntegrationPointTopology
\*---------------------------------------------------------------------------*/

class pointCentredIntegrationPointTopology
:
    public integrationPointTopology
{
    //- Mapping from cells to point-based integration-point indices
    //  Built once at construction
    labelListList cellToIntegrationPointIDs_;

public:

    //- Runtime type information
    TypeName("pointCentredIntegrationPointTopology");

    //- Construct from mesh
    explicit pointCentredIntegrationPointTopology(const fvMesh& mesh)
    :
        integrationPointTopology(mesh),
        cellToIntegrationPointIDs_(mesh.nCells())
    {
        const labelListList& cellPoints = mesh.cellPoints();

        forAll(cellPoints, cellI)
        {
            const labelList& pts = cellPoints[cellI];
            cellToIntegrationPointIDs_[cellI] = pts;
        }

        // Boundary faces are dealt separately and should not be included here
    }

    //- Destructor
    virtual ~pointCentredIntegrationPointTopology() = default;

    //- Total number of integration points (one per point)
    virtual label nIntegrationPoints() const override
    {
        return mesh_.nPoints();
    }

    //- Integration-point indices associated with a given cell
    virtual const labelUList cellIntegrationPointIDs
    (
        const label cellI
    ) const override
    {
        return cellToIntegrationPointIDs_[cellI];
    }

    //- Boundary integration patches not dealt with
    virtual bool boundaryAware() const override
    {
        return false;
    }

    //- Points are shared between cells → uniqueness required
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
