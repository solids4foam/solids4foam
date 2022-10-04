/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "globalPointIndices.H"
#include "globalMeshData.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "processorPointPatchFields.H"
#include "ListMaxOp.H"
#include "ListSumOp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPointIndices, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalPointIndices::globalPointIndices(const polyMesh& mesh)
:
    mesh_(mesh),
    ownedByThisProc_(mesh.nPoints(), true),
    localToGlobalPointMap_(mesh.nPoints(), 0),
    stencilSizeOwned_(mesh.nPoints(), 0),
    stencilSizeNotOwned_(mesh.nPoints(), 0)
{
    // There are three categories of points on each proc:
    // 1. global points: shared by more than two procs
    // 2. processor boundary points: shared by exactly two procs
    // 3. normal points: not 1 and not 2
    // For points 1 and 2, the processor with the lowest proc ID is the global
    // owner

    // Construct the ownedByThisProc list:
    // - set to true for all points
    // - for all proc patch, if the neighProcID is less than myProcID
    //   then set all patch points to false
    // - for all global points, set to true for lowest procID

    const labelList& sharedPointLabels = mesh_.globalData().sharedPointLabels();
    const pointMesh& pMesh(pointMesh::New(mesh_));

    // Correct points on processor patches
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& ppatch = mesh_.boundaryMesh()[patchI];

        if (ppatch.type() == processorPolyPatch::typeName)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(ppatch);

            const label neiProcNo = procPatch.neighbProcNo();

            // If neiProcNo is lower than myProcNo then neiProc owns the points
            if (neiProcNo < Pstream::myProcNo())
            {
                const labelList& meshPoints = ppatch.meshPoints();

                forAll(meshPoints, pI)
                {
                    const label pointID = meshPoints[pI];
                    ownedByThisProc_[pointID] = false;
                }
            }
        }
    }

    // Correct global points

    // Construct a list of all global points with the lowest processor index
    // which owns each point
    labelList globalPointOwner
    (
        mesh_.globalData().nGlobalPoints(), Pstream::nProcs()
    );
    const labelList& sharedPointAddr = mesh_.globalData().sharedPointAddr();
    forAll(sharedPointAddr, gpI)
    {
        const label globalPointID = sharedPointAddr[gpI];
        globalPointOwner[globalPointID] = Pstream::myProcNo();
    }
    reduce(globalPointOwner, minOp<labelList>());

    // Correct ownedByThisProc for global points

    forAll(sharedPointLabels, gpI)
    {
        const label globaPointID = sharedPointAddr[gpI];
        const label pointID = sharedPointLabels[gpI];

        if (globalPointOwner[globaPointID] == Pstream::myProcNo())
        {
            if (debug)
            {
                Pout<< "Proc " << Pstream::myProcNo() << " owns global point "
                    << mesh_.points()[pointID] << endl;
            }

            ownedByThisProc_[pointID] = true; 
        }
        else
        {
            ownedByThisProc_[pointID] = false;
        }
    }


    // Construct localToGlobalPointMap
    // - count the number of point owned by each proc
    // - scatter this list to all procs
    // - determine the starting global point ID for each proc
    // - create the map on each proc using this starting global ID 

    // Count the number of point owned by each proc
    label nPointsOwnByThisProc = 0;
    forAll(ownedByThisProc_, pointI)
    {
        if (ownedByThisProc_[pointI])
        {
            nPointsOwnByThisProc++;
        }
    }

    // Scater the list of points on each proc to all procs
    labelList nPointsOwnByThisProcList(Pstream::nProcs(), 0);
    nPointsOwnByThisProcList[Pstream::myProcNo()] = nPointsOwnByThisProc;
    Pstream::gatherList(nPointsOwnByThisProcList);
    Pstream::scatterList(nPointsOwnByThisProcList);

    // Determine the starting global index for the current proc by summing the
    // number of points owned by processors with a lower index
    label globalIndexStart = 0;
    forAll(nPointsOwnByThisProcList, procI)
    {
        if (procI == Pstream::myProcNo())
        {
            break;
        }

        globalIndexStart += nPointsOwnByThisProcList[procI];
    }

    if (debug)
    {
        Pout<< "globalIndexStart = " << int(globalIndexStart) << endl;
    }

    // Create map starting from globalIndexStart
    int localPointI = 0;
    forAll(localToGlobalPointMap_, pointI)
    {
        if (ownedByThisProc_[pointI])
        {
            localToGlobalPointMap_[pointI] = globalIndexStart + localPointI;
            localPointI++;
        }
    }

    // Synchronise the off-proc global indices
    // - sync the processor patch points first
    // - then sync the global points

    // Send processor point fields
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].type() == "processor")
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            const label neiProcNo = procPatch.neighbProcNo();

            // The processor with the lower ID owns the points and sends the
            // global point indices to the processor with the higher ID
            if (Pstream::myProcNo() < neiProcNo)
            {
                // Prepare field to send
                labelList toSend(procPatch.nPoints());

                // Insert global point indices
                const labelList& meshPoints = procPatch.meshPoints();
                forAll(meshPoints, pI)
                {
                    const label pointID = meshPoints[pI];
                    toSend[pI] = localToGlobalPointMap_[pointID];
                }

                // Send
                OPstream::write
                (
#ifdef OPENFOAMESIORFOUNDATION
                    Pstream::commsTypes::blocking,
#else
                    Pstream::blocking,
#endif
                    procPatch.neighbProcNo(),
                    reinterpret_cast<const char*>(toSend.begin()),
                    toSend.byteSize()
                );
            }
        }
    }

    // Receive processor point fields
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].type() == "processor")
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);
            const processorPointPatch& procPointPatch =
                refCast<const processorPointPatch>(pMesh.boundary()[patchI]);

            const label neiProcNo = procPatch.neighbProcNo();

            // The processor with the higher ID receives
            if (Pstream::myProcNo() > neiProcNo)
            {
                // Prepare field to receive
                labelList toReceive(procPatch.nPoints());

                // Receive
                IPstream::read
                ( 
#ifdef OPENFOAMESIORFOUNDATION
                    Pstream::commsTypes::blocking,
#else
                    Pstream::blocking,
#endif
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>(toReceive.begin()),
                    toReceive.byteSize()
                );

                // Insert indices into local map
                // Take care as the point ordering is reversed on the receiving
                // side
                #ifdef FOAMEXTEND
                    const labelList& meshPoints = procPatch.meshPoints();
                #else
                    const labelList& meshPoints =
                        procPointPatch.reverseMeshPoints();
                #endif
                forAll(meshPoints, pI)
                {
                    const label pointID = meshPoints[pI];
                    localToGlobalPointMap_[pointID] = toReceive[pI];
                }
            }
        }
    }


    // Sync global points

    // Correct global points

    if (debug)
    {
        Pout<< "sharedPointAddr = " << sharedPointAddr << nl
            << "sharedPointLabels = " << sharedPointLabels << endl;
    }

    {
        // Prepare a list of all global points, i.e. once synced, this list will
        // be the same on all procs
        labelList gpf(mesh_.globalData().nGlobalPoints(), -1);

        // Loop through all global points on this proc
        forAll(sharedPointAddr, globalPointI)
        {
            // Local point ID
            const label pointID = sharedPointLabels[globalPointI];

            // Insert global point ID, if this proc owns the point 
            if (ownedByThisProc_[pointID])
            {
                gpf[sharedPointAddr[globalPointI]] =
                    localToGlobalPointMap_[pointID];
            }
        }

        // Take the maximum across all procs, i.e. the point global index
        reduce(gpf, ListMaxOp<label>());

        if (debug)
        {
            Pout<< "global point indices = " << gpf << endl;
        }

        // Extract local data
        forAll(sharedPointAddr, globalPointI)
        {
            // Local point ID
            const label pointID = sharedPointLabels[globalPointI];

            // Only the procs that do not own the point need to update their map
            if (!ownedByThisProc_[pointID])
            {
                localToGlobalPointMap_[pointID] =
                    gpf[sharedPointAddr[globalPointI]];

                if (debug)
                {
                    Pout<< "Updating localToGlobalPointMap[" << pointID << "]:"
                        << nl << "    globalPointI = " << globalPointI << nl
                        << "    sharedPointAddr = "
                        << sharedPointAddr[globalPointI] << nl
                        << "    gpf = " << gpf[sharedPointAddr[globalPointI]]
                        << endl;
                }
            }
        }
    }


    if (debug)
    {
        // Write ownedByThisProc field
        pointScalarField ownedByThisProc
        (
            IOobject
            (
                "ownedByThisProc",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", dimless, -1)
        );

        forAll(ownedByThisProc, i)
        {
            ownedByThisProc[i] = scalar(int(ownedByThisProc_[i]));
        }
        ownedByThisProc.write();

        // Write global point index field
        pointScalarField t
        (
            IOobject
            (
                "globalIDs",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", dimless, -1)
        );

        forAll(t, i)
        {
            t[i] = scalar(int(localToGlobalPointMap_[i]));
        }
        t.write();
    }


    // Create the stencilSize lists
    // For each point:
    // - on-core: count all points which share a cell with the point, and
    //   add 1 for the point itself
    // - off-core: receive on-core counts across processor boundaries and via
    //   global points. Note: the stencilSizeNotOwned will not be exact, instead,
    //   it will be a conservative estimate (some points counted twice or more)

    // Check local points
    const labelListList& pointCells = mesh_.pointCells();
    const labelListList& cellPoints = mesh_.cellPoints();
    forAll(stencilSizeOwned_, pointI)
    {
        const labelList& curPointCells = pointCells[pointI];
        labelHashSet stencilOwned;
        labelHashSet stencilNotOwned;

        forAll(curPointCells, pcI)
        {
            const label cellID = curPointCells[pcI];
            const labelList& curCellPoints = cellPoints[cellID];
            forAll(curCellPoints, cpI)
            {
                const label neiPointID = curCellPoints[cpI];

                if (ownedByThisProc_[neiPointID])
                {
                    if (!stencilOwned.found(neiPointID))
                    {
                        stencilOwned.insert(neiPointID);
                    }
                }
                else
                {
                    if (!stencilNotOwned.found(neiPointID))
                    {
                        stencilNotOwned.insert(neiPointID);
                    }
                }                        
            }
        }

        stencilSizeOwned_[pointI] += stencilOwned.size();
        stencilSizeNotOwned_[pointI] += stencilNotOwned.size();
    }

    // Update stencilSizeNotOwned for points on processer boundaries

    // Send local fields
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].type() == "processor")
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            // Prepare field to send
            labelList toSend(procPatch.nPoints());

            // Insert local stencilSizes
            const labelList& meshPoints = procPatch.meshPoints();
            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                toSend[pI] =
                    stencilSizeOwned_[pointID] + stencilSizeNotOwned_[pointID];
            }

            // Send
            OPstream::write
            (
#ifdef OPENFOAMESIORFOUNDATION
                Pstream::commsTypes::blocking,
#else
                Pstream::blocking,
#endif
                procPatch.neighbProcNo(),
                reinterpret_cast<const char*>(toSend.begin()),
                toSend.byteSize()
            );
        }
    }

    // Receive processor fields
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].type() == "processor")
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            // Prepare field to receive
            labelList toReceive(procPatch.nPoints());

            // Receive
            IPstream::read
            (
#ifdef OPENFOAMESIORFOUNDATION
                Pstream::commsTypes::blocking,
#else
                Pstream::blocking,
#endif
                procPatch.neighbProcNo(),
                reinterpret_cast<char*>(toReceive.begin()),
                toReceive.byteSize()
            );

            // Add neiProc sizes to the local sizes: this will mean some points
            // are doubly counted, and so it will be a conservative estimate
            const labelList& meshPoints = procPatch.meshPoints();
            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                // Increment owned and not-owned fields
                stencilSizeOwned_[pointID] += toReceive[pI];
                stencilSizeNotOwned_[pointID] += toReceive[pI];
            }
        }
    }

    // Correct global points for stencilSizeNotOwned
    {
        // Prepare a list of all global points, i.e. once synced, this list will
        // be the same on all procs
        labelList gpf(mesh_.globalData().nGlobalPoints(), 0);

        // Loop through all global points on this proc
        forAll(sharedPointAddr, globalPointI)
        {
            // Local point ID
            const label pointID = sharedPointLabels[globalPointI];

            // Insert owned and not-owned stencil sizes
            gpf[sharedPointAddr[globalPointI]] =
                stencilSizeOwned_[pointID] + stencilSizeNotOwned_[pointID];
        }

        // Sum across all procs. Once again, this will be conservative as it
        // will include some points more than once
        reduce(gpf, ListSumOp<label>());

        // Extract local data
        forAll(sharedPointAddr, globalPointI)
        {
            // Local point ID
            const label pointID = sharedPointLabels[globalPointI];

            // Update the off-core stencil sizes
            stencilSizeNotOwned_[pointID] += gpf[sharedPointAddr[globalPointI]];
        }
    }

    if (debug)
    {
        // Write global point index field
        const pointMesh& pMesh(pointMesh::New(mesh_));
        pointScalarField a
        (
            IOobject
            (
                "stencilSizeOwned",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", dimless, -1)
        );
        pointScalarField b
        (
            IOobject
            (
                "stencilSizeNotOwned",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedScalar("0", dimless, -1)
        );

        forAll(a, i)
        {
            a[i] = scalar(int(stencilSizeOwned_[i]));
            b[i] = scalar(int(stencilSizeNotOwned_[i]));
        }
        a.write();
        b.write();
    }

    // if (debug)
    // {
    //     label maxStencilSize = 0;
    //     label minStencilSize = 1000;
    //     label avStencilSize = 0;
    //     forAll(stencilSizeOwned_, pointI)
    //     {
    //         maxStencilSize = max(maxStencilSize, stencilSizeOwned_[pointI]);
    //         minStencilSize = min(minStencilSize, stencilSizeOwned_[pointI]);
    //         avStencilSize += stencilSizeOwned_[pointI];
    //     }
    //     avStencilSize /= mesh_.nPoints();
    //     Info<< "maxStencilSize = " << maxStencilSize
    //         << ", minStencilSize = " << minStencilSize
    //         << ", avStencilSize = " << avStencilSize << endl;
    // }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::globalPointIndices::stencilSize() const
{
    labelList result(stencilSizeOwned_.size(), 0);

    forAll(result, i)
    {
        result[i] = stencilSizeOwned_[i] + stencilSizeNotOwned_[i];
    }

    return result;
}

// ************************************************************************* //
