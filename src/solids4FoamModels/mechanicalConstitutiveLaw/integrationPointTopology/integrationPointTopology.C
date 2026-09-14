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

#include "integrationPointTopology.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(integrationPointTopology, 0);
    defineRunTimeSelectionTable
    (
        integrationPointTopology, fvMesh
    );
}

// * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::integrationPointTopology>
Foam::integrationPointTopology::New
(
    const word& topologyType,
    const fvMesh& mesh
)
{
#if (OPENFOAM >= 2112)
    auto* ctorPtr = fvMeshConstructorTable(topologyType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown integrationPointTopology type: "
            << topologyType << nl << nl
            << "Valid integrationPointTopology types are:" << nl
            << fvMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
#else
    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(topologyType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown integrationPointTopology type: "
            << topologyType << nl << nl
            << "Valid integrationPointTopology types are:" << nl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    auto* ctorPtr = cstrIter();
#endif

    return autoPtr<integrationPointTopology>(ctorPtr(mesh));
}


// ************************************************************************* //
