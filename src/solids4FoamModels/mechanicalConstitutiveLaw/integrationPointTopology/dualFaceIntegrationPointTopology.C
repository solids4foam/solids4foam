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

#include "dualFaceIntegrationPointTopology.H"
#include "ListOps.H"
#include "SubList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dualFaceIntegrationPointTopology, 0);

    // Note: this topology is not added to the run-time selection table, as it
    // cannot be constructed from the mesh alone: it also requires the
    // dual-face-to-cell map
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dualFaceIntegrationPointTopology::dualFaceIntegrationPointTopology
(
    const fvMesh& mesh,
    const labelUList& dualFaceToCell,
    const label nInternalDualFaces
)
:
    integrationPointTopology(mesh),
    nInternalDualFaces_(nInternalDualFaces),
    cellToIntegrationPointIDs_()
{
    if (nInternalDualFaces_ < 0 || nInternalDualFaces_ > dualFaceToCell.size())
    {
        FatalErrorInFunction
            << "The number of internal dual faces (" << nInternalDualFaces_
            << ") is not consistent with the size of the dualFaceToCell map ("
            << dualFaceToCell.size() << ")."
            << exit(FatalError);
    }

    // Every internal dual face must map to a primary cell, otherwise its
    // mechanical constitutive law is undefined. dualMeshToMeshMap already
    // enforces this, but the map is passed here as a plain list so check again
    for (label dualFaceI = 0; dualFaceI < nInternalDualFaces_; ++dualFaceI)
    {
        if (dualFaceToCell[dualFaceI] < 0)
        {
            FatalErrorInFunction
                << "Internal dual face " << dualFaceI << " does not map to a "
                << "primary mesh cell, so its mechanical constitutive law "
                << "cannot be determined."
                << exit(FatalError);
        }
    }

    // Invert the map to give the internal dual faces within each primary cell
    cellToIntegrationPointIDs_ =
        invertOneToMany
        (
            mesh_.nCells(),
            SubList<label>(dualFaceToCell, nInternalDualFaces_)
        );
}


// ************************************************************************* //
