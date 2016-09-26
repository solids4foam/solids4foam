/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fluidModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidModel, 0);
    defineRunTimeSelectionTable(fluidModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

 void Foam::fluidModel::calcGlobalFaceZones() const
 {
     // Find global face zones
     if (globalFaceZonesPtr_)
     {
         FatalErrorIn
         (
             "void fluidModel::calcGlobalFaceZones() const"
         )
             << "Global face zones already fonud"
                 << abort(FatalError);
     }

     SLList<label> globalFaceZonesSet;

     const faceZoneMesh& faceZones = mesh().faceZones();

     forAll(faceZones, zoneI)
     {
         const faceZone& curFaceZone = faceZones[zoneI];

         forAll(curFaceZone, faceI)
         {
//              if unused face exist
             if (curFaceZone[faceI] >= mesh().nFaces())
             {
                 globalFaceZonesSet.insert(zoneI);
                 break;
             }
         }
     }

     globalFaceZonesPtr_ = new labelList(globalFaceZonesSet);
 }


 void Foam::fluidModel::calcGlobalToLocalFaceZonePointMap() const
 {
     // Find global face zones
     if (globalToLocalFaceZonePointMapPtr_)
     {
         FatalErrorIn
         (
             "void fluidModel::calcGlobalToLocalFaceZonePointMap() const"
         )
             << "Global to local face zones point map already exists"
                 << abort(FatalError);
     }

     globalToLocalFaceZonePointMapPtr_ =
         new labelListList(globalFaceZones().size());

     labelListList& globalToLocalFaceZonePointMap =
         *globalToLocalFaceZonePointMapPtr_;

     forAll(globalFaceZones(), zoneI)
     {
         label curZoneID = globalFaceZones()[zoneI];

         labelList curMap(mesh().faceZones()[curZoneID]().nPoints(), -1);

         vectorField fzGlobalPoints =
             mesh().faceZones()[curZoneID]().localPoints();

         //- set all slave points to zero because only the master order is used
         if(!Pstream::master())
         {
             fzGlobalPoints *= 0.0;
         }

         //- pass points to all procs
         reduce(fzGlobalPoints, sumOp<vectorField>());

         //- now every proc has the master's list of FZ points
         //- every proc must now find the mapping from their local FZ points to
         //- the global FZ points

         const vectorField& fzLocalPoints =
             mesh().faceZones()[curZoneID]().localPoints();

         const edgeList& fzLocalEdges =
             mesh().faceZones()[curZoneID]().edges();

         const labelListList& fzPointEdges =
             mesh().faceZones()[curZoneID]().pointEdges();

         scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

         forAll(minEdgeLength, pI)
         {
             const labelList& curPointEdges = fzPointEdges[pI];

             forAll(curPointEdges, eI)
             {
                 scalar Le = fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
                 if (Le < minEdgeLength[pI])
                 {
                     minEdgeLength[pI] = Le;
                 }
             }
         }

         forAll(fzGlobalPoints, globalPointI)
         {
             boolList visited(fzLocalPoints.size(), false);

             forAll(fzLocalPoints, procPointI)
             {
                 if (!visited[procPointI])
                 {
                     visited[procPointI] = true;

                     label nextPoint = procPointI;

                     scalar curDist =
                         mag
                         (
                             fzLocalPoints[nextPoint]
                           - fzGlobalPoints[globalPointI]
                         );

                     if (curDist < 1e-4*minEdgeLength[nextPoint])
                     {
                         curMap[globalPointI] = nextPoint;
                         break;
                     }

                     label found = false;

                     while (nextPoint != -1)
                     {
                         const labelList& nextPointEdges =
                             fzPointEdges[nextPoint];

                         scalar minDist = GREAT;
                         label index = -1;
                         forAll(nextPointEdges, edgeI)
                         {
                             label curNgbPoint =
                                 fzLocalEdges[nextPointEdges[edgeI]]
                                .otherVertex(nextPoint);

                             if (!visited[curNgbPoint])
                             {
                                 visited[curNgbPoint] = true;

                                 scalar curDist =
                                     mag
                                     (
                                         fzLocalPoints[curNgbPoint]
                                       - fzGlobalPoints[globalPointI]
                                     );

                                 if (curDist < 1e-4*minEdgeLength[curNgbPoint])
                                 {
                                     curMap[globalPointI] = curNgbPoint;
                                     found = true;
                                     break;
                                 }
                                 else if (curDist < minDist)
                                 {
                                     minDist = curDist;
                                     index = curNgbPoint;
                                 }
                             }
                         }

                         nextPoint = index;
                     }

                     if (found)
                     {
                         break;
                     }
                 }
             }
         }


         forAll(curMap, globalPointI)
         {
             if (curMap[globalPointI] == -1)
             {
                 FatalErrorIn
                 (
                     "fluidModel::calcGlobalToLocalFaceZonePointMap()"
                 )
                     << "local to global face zone point map is not correct"
                         << abort(FatalError);
             }
         }

         globalToLocalFaceZonePointMap[zoneI] = curMap;
     }
 }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidModel::fluidModel
(
    const word& type,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fluidProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    fluidProperties_(subDict(type + "Coeffs")),
     globalFaceZonesPtr_(NULL),
     globalToLocalFaceZonePointMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidModel::~fluidModel()
{
     deleteDemandDrivenData(globalFaceZonesPtr_);
     deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 const Foam::labelList&
 Foam::fluidModel::globalFaceZones() const
 {
     if (!globalFaceZonesPtr_)
     {
         calcGlobalFaceZones();
     }

     return *globalFaceZonesPtr_;
 }

 const Foam::labelListList&
 Foam::fluidModel::globalToLocalFaceZonePointMap() const
 {
     if (!globalToLocalFaceZonePointMapPtr_)
     {
         calcGlobalToLocalFaceZonePointMap();
     }

     return *globalToLocalFaceZonePointMapPtr_;
 }

bool Foam::fluidModel::read()
{
    if (regIOobject::read())
    {
        fluidProperties_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
