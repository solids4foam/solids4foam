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
    Base class defining the topology of integration points used for mechanical
    constitutive-law evaluation.

    The integrationPointTopology provides a mapping from mesh entities
    (currently cells) to integration-point indices belonging to the
    constitutive working set. These indices are used by the
    mechanicalConstitutiveLawManager to determine where each mechanical
    constitutive law is evaluated.

    The constitutive working set is owned by the solid model and may represent:
      - cell-centred values for low-order finite-volume formulations, or
      - integration-point values stored in a compact form (e.g.
        CompactListList<Type>::values()) for higher-order formulations.

    This class deliberately contains no material or constitutive information.
    It only describes how integration points are associated with mesh entities,
    allowing the mechanical constitutive law manager to remain agnostic to the
    discretisation order, quadrature scheme, and data storage layout.

    Overlap or multiplicity of integration-point indices across entities is
    permitted and depends on the chosen discretisation (e.g. cell-, face-, or
    vertex-based formulations).

    See also:
    cellCentredIntegrationPointTopology
    mechanicalConstitutiveLawManager

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef integrationPointTopology_H
#define integrationPointTopology_H

#include "labelList.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class integrationPointTopology
\*---------------------------------------------------------------------------*/

class integrationPointTopology
{
protected:

    //- Reference to mesh (needed for size consistency checks etc.)
    const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("integrationPointTopology");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        integrationPointTopology,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

    //- Selector
    static autoPtr<integrationPointTopology> New
    (
        const word& topologyType,
        const fvMesh& mesh
    );

    //- Constructor
    explicit integrationPointTopology(const fvMesh& mesh)
    :
        mesh_(mesh)
    {}

    //- Destructor
    virtual ~integrationPointTopology() = default;

    //- Total number of integration points in the constitutive working set
    virtual label nIntegrationPoints() const = 0;

    //- Integration-point indices associated with a given cell
    //  The returned object is a lightweight, non-owning view (labelUList)
    //  into integration-point addressing stored by the topology.
    //  No data is copied; the underlying storage remains owned by the
    //  topology object and must outlive the use of the returned view.
    virtual const labelUList cellIntegrationPointIDs
    (
        const label cellI
    ) const = 0;

    //- Does this topology store boundary integration points per patch
    virtual bool boundaryAware() const = 0;

    //- Integration-point indices for the faces of a boundary patch, in patch
    //  face order, or an empty list where the topology has none.
    //
    //  This is a separate question from boundaryAware(), and the two must not
    //  be confused. boundaryAware() asks whether the topology keeps a
    //  constitutive state per patch, which a cell-centred topology does, for
    //  the GeometricField updates that evaluate a law on patch faces. This
    //  asks whether the topology's *flat index space* contains slots for
    //  boundary faces, which a cell-centred topology's does not: it has one
    //  integration point per cell and nothing else.
    //
    //  Returning empty therefore means "no boundary slots to fill", and is the
    //  right answer for every topology whose integration points are interior
    //  quantities. A face-centred topology overrides it, because its index
    //  space is the mesh faces, boundary faces included.
    //
    //  Returned by value: a view onto storage owned by the topology would be
    //  invalidated by the next call, which is a trap for any caller holding
    //  two patches at once
    virtual labelList boundaryIntegrationPointIDs(const label patchI) const
    {
        return labelList();
    }

    //- Can this topology carry a fourth-order (mat66) material tangent?
    //  This is true only where the integration points are the locations at
    //  which a Jacobian operator evaluates fluxes, i.e. face-like rather than
    //  cell-like topologies. Cell-, point- and cell-quadrature-based
    //  topologies never carry a fourth-order tangent, because the operators
    //  that consume one are assembled from face fluxes.
    //  Note that this is a structural property of the topology alone. The
    //  separate requirement that each integration point belong to exactly one
    //  mechanical constitutive law is expressed by
    //  requiresUniqueIntegrationPointsPerMaterial().
    virtual bool supportsFourthOrderTangent() const
    {
        return false;
    }

    //- Should integration-point indices be unique when aggregating across
    //  multiple cells belonging to the same material.
    //  For cell-centred topologies, each cell maps to a unique integration
    //  point and duplicates cannot arise, so this should return false.
    //  For face- or point-centred topologies, the same integration point may
    //  be associated with multiple cells. In such cases, this function should
    //  return true to ensure that each integration point is evaluated exactly
    //  once per material law.
    virtual bool requiresUniqueIntegrationPointsPerMaterial() const = 0;

    //- Whether an integration point is a mesh face, with the integration
    //  point index being the global face index.
    //
    //  A value the user gives per cell reaches a face by interpolation, not
    //  by taking one of the two adjacent cell values, and a topology that
    //  answers true is asking to be filled that way. Answering false means
    //  a cell's value is broadcast to the integration points of that cell,
    //  which is right for a topology whose points sit inside one cell.
    //
    //  Note the index requirement: this is only true of a topology whose
    //  integration point index *is* the face index. A topology holding
    //  several points per face, as the compact face one does, indexes into a
    //  per-face block instead and must answer false
    virtual bool integrationPointsAreFaces() const
    {
        return false;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
