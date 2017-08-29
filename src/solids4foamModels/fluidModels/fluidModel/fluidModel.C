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
#include "fv.H"
#include "elasticWallPressureFvPatchScalarField.H"
#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

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
        )   << "Global face zones already fonud"
            << abort(FatalError);
    }

    SLList<label> globalFaceZonesSet;

    const faceZoneMesh& faceZones = mesh().faceZones();

    forAll(faceZones, zoneI)
    {
        const faceZone& curFaceZone = faceZones[zoneI];

        bool globalFaceZone = false;

        forAll(curFaceZone, faceI)
        {
            // If unused face exist, then this must be a global face zone
            if (curFaceZone[faceI] >= mesh().nFaces())
            {
                globalFaceZonesSet.insert(zoneI);
                break;
            }
        }

        reduce(globalFaceZone, orOp<bool>());

        if (globalFaceZone)
        {
            globalFaceZonesSet.insert(zoneI);
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
        )   << "Global to local face zones point map already exists"
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
            fzGlobalPoints = vector::zero;
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fluidModel::updateRobinFsiInterface
(
    const volScalarField& p,
    const volVectorField& U,
    surfaceScalarField& phi,
    surfaceScalarField& rAUf
)
{
    forAll(p.boundaryField(), patchI)
    {
        if
        (
            (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
             && isA<elasticSlipWallVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
         || (
                isA<elasticWallPressureFvPatchScalarField>
                (
                    p.boundaryField()[patchI]
                )
             && isA<elasticWallVelocityFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
        )
        {
            const word ddtScheme =
                mesh().schemesDict().ddtScheme("ddt(" + U.name() +')');

            if (ddtScheme == fv::EulerDdtScheme<vector>::typeName)
            {
                phi.boundaryField()[patchI] =
                    phi.oldTime().boundaryField()[patchI];
                rAUf.boundaryField()[patchI] = runTime().deltaT().value();
            }
            else if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
            {
                if (runTime().timeIndex() == 1)
                {
                    phi.boundaryField()[patchI] =
                        phi.oldTime().boundaryField()[patchI];
                    rAUf.boundaryField()[patchI] = runTime().deltaT().value();

                    phi.oldTime().oldTime();
                }
                else
                {
                    scalar deltaT = runTime().deltaT().value();
                    scalar deltaT0 = runTime().deltaT0().value();

                    scalar Cn = 1 + deltaT/(deltaT + deltaT0);
                    scalar Coo = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
                    scalar Co = Cn + Coo;

                    phi.boundaryField()[patchI] =
                        (Co/Cn)*phi.oldTime().boundaryField()[patchI]
                      - (Coo/Cn)
                       *phi.oldTime().oldTime().boundaryField()[patchI];

                    rAUf.boundaryField()[patchI] =
                        runTime().deltaT().value()/Cn;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidModel::fluidModel
(
    const word& type,
    Time& runTime,
    const word& region
)
:
    physicsModel(type, runTime),
    IOdictionary
    (
        IOobject
        (
            "fluidProperties",
            // If region == "region0" then read from the main case
            // Otherwise, read from the region/sub-mesh directory e.g.
            // constant/fluid or constant/solid
            bool(region == dynamicFvMesh::defaultRegion)
          ? fileName(runTime.constant())
          : fileName(runTime.constant()/region),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    meshPtr_
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                region,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    ),
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

Foam::autoPtr<Foam::fluidModel> Foam::fluidModel::New
(
    Time& runTime,
    const word& region
)
{
    word fluidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary fluidProperties
        (
            IOobject
            (
                "fluidProperties",
                //runTime.constant(),
                fileName("constant"/region),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        fluidProperties.lookup("fluidModel")
            >> fluidModelTypeName;
    }

    Info<< nl << "Selecting fluidModel " << fluidModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fluidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fluidModel::New(Time&, const word&)"
        )   << "Unknown fluidModel type " << fluidModelTypeName
            << endl << endl
            << "Valid fluidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<fluidModel>(cstrIter()(runTime, region));
}


const Foam::labelList& Foam::fluidModel::globalFaceZones() const
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
