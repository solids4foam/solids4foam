/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

InClass
    solidPolyMeshLduAddressing

\*---------------------------------------------------------------------------*/

#include "solidPolyMeshLduAddressing.H"
#include "globalMeshData.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPolyMeshLduAddressing::solidPolyMeshLduAddressing
(
    const solidPolyMesh& mesh
)
:
    lduAddressing(mesh.nVariables()),
    lowerAddr_(),
    upperAddr_(),
    patchAddr_(mesh.boundary().boundaryMesh().size()),
    patchSchedule_(mesh.globalData().patchSchedule())
{
    const edgeList& implicitBonds = mesh.implicitBonds();

    upperAddr_.setSize(implicitBonds.size());
    lowerAddr_.setSize(implicitBonds.size());

    forAll(lowerAddr_, implicitBondI)
    {
        lowerAddr_[implicitBondI] = implicitBonds[implicitBondI].start();
        upperAddr_[implicitBondI] = implicitBonds[implicitBondI].end();
    }

    forAll(patchAddr_, patchI)
    {
        patchAddr_[patchI] = mesh.boundary().boundaryMesh()[patchI].faceCells();
    }
}


// ************************************************************************* //
