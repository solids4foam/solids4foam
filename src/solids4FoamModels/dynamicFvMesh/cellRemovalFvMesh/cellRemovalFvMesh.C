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

\*---------------------------------------------------------------------------*/

#include "cellRemovalFvMesh.H"
#include "foamTime.H"
#include "regionSplit.H"
#include "mapPolyMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellRemovalFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, cellRemovalFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::cellRemovalFvMesh::calcMeshMod() const
{
    if (meshModPtr_)
    {
        FatalErrorIn("void Foam::cellRemovalFvMesh::calcMeshMod() const")
            << "pointe already set" << abort(FatalError);
    }

    meshModPtr_ = new directTopoChange(*this);
}


Foam::directTopoChange& Foam::cellRemovalFvMesh::meshModifier()
{
    if (!meshModPtr_)
    {
        calcMeshMod();
    }

    return *meshModPtr_;
}


void Foam::cellRemovalFvMesh::clearOut() const
{
    deleteDemandDrivenData(meshModPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellRemovalFvMesh::cellRemovalFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    writeMesh_(dict_.lookupOrDefault<Switch>("writeMesh", false)),
    removeDeadCells_(dict_.lookupOrDefault<Switch>("removeDeadCells", true)),
    meshModPtr_(NULL),
    cellRemover_(*this, true),
    lawPtr_(cellRemovalLaw::New("law", *this, dict_.subDict("law")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellRemovalFvMesh::~cellRemovalFvMesh()
{
    deleteDemandDrivenData(meshModPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellRemovalFvMesh::update()
{
    // Check if there are cells to remove
    const labelField cellsToRemove = lawPtr_->cellsToRemove();

    const label nCellsToRemove =
        returnReduce(cellsToRemove.size(), sumOp<int>());

    if (nCellsToRemove)
    {
        // Exposed faces will be inserted into the open patch

        Info<< nl << "There are " << nCellsToRemove
            << " cells to be removed" << endl;

        // Find faces that will be exposed
        // These faces will become boundary faces

        const labelList facesToExpose
        (
            cellRemover_.getExposedFaces(cellsToRemove)
        );

        Info<< "There are " << facesToExpose.size()
            << " internal faces that will be exposed" << endl;

        // Clear the directTopoChanger before giving it to the cell remover
        clearOut();

        // Set actions in cell remover
        cellRemover_.setRefinement
        (
            cellsToRemove,
            facesToExpose,
            labelList(facesToExpose.size(), lawPtr_->exposedFacesPatchID()),
            meshModifier()
        );

        // Change the mesh
        Info<< "Performing mesh change" << endl;
        autoPtr<mapPolyMesh> map = meshModifier().changeMesh(*this, true);

        // Update cell remover map
        cellRemover_.updateMesh(map());

        // Update mesh fields e.g. U, sigma, etc.
        updateMesh(map);

        // Update fields on newly exposed faces
        Info<< "Updating field values on newly exposed faces" << endl;
        updateVolFieldsExposedFaces<scalar>(map, facesToExpose);
        updateVolFieldsExposedFaces<vector>(map, facesToExpose);
        updateVolFieldsExposedFaces<tensor>(map, facesToExpose);
        updateVolFieldsExposedFaces<symmTensor>(map, facesToExpose);
        updateVolFieldsExposedFaces<diagTensor>(map, facesToExpose);
        updateVolFieldsExposedFaces<sphericalTensor>(map, facesToExpose);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            fvMesh::movePoints(map().preMotionPoints());
        }

        // Remove any dead cells
        if (removeDeadCells_)
        {
            Info<< "Checking for dead cells" << endl;

            regionSplit regions(*this);

            // Make nCellsInRegion list

            labelList nCellsInRegion(regions.nRegions(), 0);

            forAll(regions, cellI)
            {
                nCellsInRegion[regions[cellI]]++;
            }

            // Find regions with less cells than the dead-cell limit
            Info<< "The number of ncells is " << nCells()<< endl;

            const label deadCellLimit = 0.1*nCells();
            labelHashSet deadCellsSet;

            forAll(regions, cellI)
            {
                if (nCellsInRegion[regions[cellI]] < deadCellLimit)
                {
                    deadCellsSet.insert(cellI);
                }
            }

            const labelList deadCells = deadCellsSet.toc();

            label ndeadCells = deadCells.size();
            reduce(ndeadCells, maxOp<label>());

            if (ndeadCells)
            {
                // Remove dead cells from the mesh

                // Clear the directTopoChanger before giving it to the cell
                // remover
                clearOut();

                // Set actions in cell remover
                cellRemover_.setRefinement
                (
                    deadCells,
                    labelList(0),    // no faces to expose
                    labelList(0),    // open patch not needed
                    meshModifier()
                );

                // Change the mesh
                Info<< "    Removing " << deadCells.size() << " dead cells"
                    << endl;

                autoPtr<mapPolyMesh> dmap =
                    meshModifier().changeMesh(*this, true);

                // Update cell remover map
                cellRemover_.updateMesh(dmap());

                // Update mesh fields e.g. U, sigma, etc.
                updateMesh(dmap);

                // Move mesh (since morphing does not do this)
                if (dmap().hasMotionPoints())
                {
                    fvMesh::movePoints(dmap().preMotionPoints());
                }
            }
            else
            {
                Info<< "    There are no dead cells" << endl;
            }
        }
    }

    // We will write the mesh if there is a topological change
    if (writeMesh_)
    {
        Info<< "Writing the mesh" << endl;
        write();
    }

    return bool(nCellsToRemove > 0);
}


// ************************************************************************* //
