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
    Integration-point topology for face-centred solid mechanics formulations.

    This topology defines one integration point per mesh face. Each cell
    maps to the integration points associated with its faces.

    Internal faces are shared by two cells and therefore appear in the
    integration-point lists of both adjacent cells. As a result,
    integration-point indices must be deduplicated when aggregating
    per-material working sets.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef faceCentredIntegrationPointTopology_H
#define faceCentredIntegrationPointTopology_H

#include "integrationPointTopology.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class faceCentredIntegrationPointTopology
\*---------------------------------------------------------------------------*/

class faceCentredIntegrationPointTopology
:
    public integrationPointTopology
{
    //- Mapping from cell index to face-based integration-point indices
    //  The indices refer to global face indices:
    //    - internal faces: [0 .. nInternalFaces-1]
    //    - boundary faces: offset internally per patch
    List<labelList> cellToIntegrationPointIDs_;

public:

    //- Runtime type information
    TypeName("faceCentredIntegrationPointTopology");

    //- Construct from mesh
    explicit faceCentredIntegrationPointTopology(const fvMesh& mesh)
    :
        integrationPointTopology(mesh),
        cellToIntegrationPointIDs_(mesh.nCells())
    {
        // Internal faces
        forAll(mesh_.owner(), faceI)
        {
            const label owner = mesh_.owner()[faceI];
            const label neigh = mesh_.neighbour()[faceI];

            cellToIntegrationPointIDs_[owner].append(faceI);
            cellToIntegrationPointIDs_[neigh].append(faceI);
        }

        // Boundary faces are dealt separately and should not be included here
    }

    //- Destructor
    virtual ~faceCentredIntegrationPointTopology() = default;

    //- Total number of integration points (one per face)
    virtual label nIntegrationPoints() const override
    {
        return mesh_.nFaces();
    }

    //- Integration-point indices associated with a given cell
    virtual const labelUList cellIntegrationPointIDs
    (
        const label cellI
    ) const override
    {
        return cellToIntegrationPointIDs_[cellI];
    }

    //- Boundary integration points are supported
    virtual bool boundaryAware() const override
    {
        return true;
    }

    //- Faces are shared between cells → uniqueness required
    virtual bool requiresUniqueIntegrationPointsPerMaterial() const override
    {
        return true;
    }

    //- Faces are flux-evaluation points, so a fourth-order tangent is
    //  meaningful here. Note that requiresUniqueIntegrationPointsPerMaterial()
    //  is true, so a fourth-order tangent is only usable with a single
    //  mechanical constitutive law: a face on a material interface has no
    //  single tangent
    virtual bool supportsFourthOrderTangent() const override
    {
        return true;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
