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
    Integration points located at the internal faces of a dual mesh.

    Used by the vertex-centred solid models, where stresses and material
    tangents are evaluated at dual-mesh faces while materials remain defined
    per primary-mesh cell.

    Each dual face lies within (or, for a boundary dual face, on the boundary
    of) exactly one primary-mesh cell, so each integration point belongs to
    exactly one mechanical constitutive law. No stress or tangent collapse is
    therefore required, and a fourth-order tangent is well defined even for a
    multi-material case.

    Only internal dual faces are integration points. The vertex-centred
    assembly routines (vfvm::divSigma and vfvm::laplacian) iterate
    dualMesh.owner(), i.e. internal dual faces only, and the primary-mesh cell
    that owns a boundary dual face already carries its own integration points.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef dualFaceIntegrationPointTopology_H
#define dualFaceIntegrationPointTopology_H

#include "integrationPointTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class dualFaceIntegrationPointTopology
\*---------------------------------------------------------------------------*/

class dualFaceIntegrationPointTopology
:
    public integrationPointTopology
{
    // Private Data

        //- Number of internal faces in the dual mesh
        const label nInternalDualFaces_;

        //- Internal dual faces lying within each primary-mesh cell.
        //  This is the inverse of the dualFaceToCell map, restricted to
        //  internal dual faces
        labelListList cellToIntegrationPointIDs_;


public:

    //- Runtime type information
    TypeName("dualFaceIntegrationPointTopology");

    //- Construct from the primary mesh, the dual-face-to-cell map and the
    //  number of internal dual faces.
    //  The map is indexed by dual face and is not retained: only its
    //  inverse is stored.
    //  Note that this topology cannot be built by the run-time selection
    //  table, which supplies the mesh alone; it must be handed to the
    //  mechanicalConstitutiveLawManager via registerTopology()
    dualFaceIntegrationPointTopology
    (
        const fvMesh& mesh,
        const labelUList& dualFaceToCell,
        const label nInternalDualFaces
    );

    //- Destructor
    virtual ~dualFaceIntegrationPointTopology() = default;

    //- Total number of integration points (one per internal dual face)
    virtual label nIntegrationPoints() const override
    {
        return nInternalDualFaces_;
    }

    //- Integration-point indices associated with a given primary-mesh cell
    virtual const labelUList cellIntegrationPointIDs
    (
        const label cellI
    ) const override
    {
        return cellToIntegrationPointIDs_[cellI];
    }

    //- Boundary integration points are not stored: a boundary dual face is
    //  already covered by the primary cell that owns it
    virtual bool boundaryAware() const override
    {
        return false;
    }

    //- Each internal dual face lies in exactly one primary cell, so an
    //  integration point cannot be reached from two cells of the same
    //  material and no de-duplication is needed
    virtual bool requiresUniqueIntegrationPointsPerMaterial() const override
    {
        return false;
    }

    //- Dual faces are the flux-evaluation points of the vertex-centred
    //  operators, and the map from a dual face to a mechanical constitutive
    //  law is single-valued, so a fourth-order tangent is well defined here
    //  even with more than one law
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
