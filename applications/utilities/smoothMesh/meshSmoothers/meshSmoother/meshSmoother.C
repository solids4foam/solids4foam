/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "meshSmoother.H"
#include "volFields.H"
#include "fvc.H"
#include "fvMesh.H"
#include "pointFields.H"
#include "globalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(meshSmoother, 0);
defineRunTimeSelectionTable(meshSmoother, dictionary);


// * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

void meshSmoother::correctBoundaryMotion
(
    pointVectorField& pointMotionD
)
{
    if (debug)
    {
        Info<< "    Correcting the boundary" << endl;
    }

    // Take references
    vectorField& pointMotionDI = pointMotionD.internalField();
    const polyMesh& mesh = pointMotionD.mesh().mesh();


    // Points on the boundary are only allowed to slide/slip along the boundary
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (!mesh.boundaryMesh()[patchI].coupled())
        {
            const vectorField& pointNormals =
                mesh.boundaryMesh()[patchI].pointNormals();
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                pointMotionDI[pointID] =
                    (I - sqr(pointNormals[pI])) & pointMotionDI[pointID];
            }
        }
    }


    //
    // Correct feature edges
    //


    // Finally, we will look for feature edges based on the angle between faces
    // For points on feature edges, the points are along allowed to slide along
    // the edges

    if (dict().lookupOrDefault<Switch>("correctFeatureEdges", false))
    {
        if (debug)
        {
            Info<< "    Correcting feature lines" << endl;
        }

        const scalar featureAngle =
            dict().lookupOrDefault<scalar>("featureAngle", 45);

        if (debug)
        {
            Info<< "        feature angle: " << featureAngle << endl;
        }

        // Multiple the angle by 0.5 because we use the point normal below
        const scalar featureAngleDot =
            Foam::cos(0.5*featureAngle*mathematicalConstant::pi/180.0);

        // Method
        // Create a globalPolyPatch for either all the boundary or
        // individually for each patch, then we can treat the patch as we would
        // in serial: then a serial algorithm could be: for each point on the
        // boundary, for each pointFace, compare the angle between the face
        // normal and the point normal; if greater than the specified feature
        // angle (or feature angle divided by two since we are using the point
        // normal) then we restrict the point motion to the face plane i.e.
        // remove the out of plane component of motion

        const label nPatches =
            returnReduce(mesh.boundaryMesh().size(), minOp<int>());

        for (label patchI = 0; patchI < nPatches; patchI++)
        {
            // Skip empty patches (this will happen because we are using a
            // sub-mesh) and skip coupled patches
            if
            (
                (
                    returnReduce
                    (
                        mesh.boundaryMesh()[patchI].size(), sumOp<int>()
                    ) == 0
                )
             || mesh.boundaryMesh()[patchI].coupled()
            )
            {
                continue;
            }

            // Contruct the local patch point motion
            vectorField localPointMotionD
                (
                    mesh.boundaryMesh()[patchI].nPoints(),
                    vector::zero
                );

            // Create local scope so we don't subsequently confuse the local
            // meshPoints with the global patch addressing
            {
                const labelList& meshPoints =
                    mesh.boundaryMesh()[patchI].meshPoints();

                forAll(localPointMotionD, pI)
                {
                    localPointMotionD[pI] = pointMotionDI[meshPoints[pI]];
                }
            }

            // Create a global patch from this local patch
            globalPolyPatch gpatch
                (
                    mesh.boundaryMesh()[patchI].name(),
                    mesh
                );

            // Map the local patch motion to the global patch
            vectorField globalPointMotionD =
                gpatch.patchPointToGlobal(localPointMotionD);

            // Take references
            const polyPatch& ppatch = gpatch.patch();
            const vectorField& faceNormals = ppatch.faceNormals();
            const vectorField& pointNormals = ppatch.pointNormals();
            const labelListList& pointFaces = ppatch.pointFaces();

            forAll(pointNormals, pI)
            {
                vector& curPointMotionD = globalPointMotionD[pI];
                const labelList& curPointFaces = pointFaces[pI];
                const vector& curPointNormal = pointNormals[pI];

                forAll(curPointFaces, fI)
                {
                    const label faceID = curPointFaces[fI];
                    const vector& curFaceNormal = faceNormals[faceID];

                    if ((curPointNormal & curFaceNormal) < featureAngleDot)
                    {
                        // Set the point motion to zero: this is a simple
                        // solution. It would be better to allow sliding along
                        // the feature edge but that would require a bit more
                        // work
                        curPointMotionD = vector::zero;
                    }
                }
            }

            // Map the global patch motion back to the local patch
            localPointMotionD = gpatch.globalPointToPatch(globalPointMotionD);

            // Map this patch field back into the pointMotionD
            {
                const labelList& meshPoints =
                    mesh.boundaryMesh()[patchI].meshPoints();

                forAll(localPointMotionD, pI)
                {
                    pointMotionDI[meshPoints[pI]] = localPointMotionD[pI];
                }
            }
        }
    }

    if (dict().found("fixedPatches"))
    {
        if (debug)
        {
            Info<< "    Correcting fixed patches" << endl;
        }

        const wordList fixedPatches
        (
            dict().lookup("fixedPatches")
        );

        forAll(fixedPatches, fpI)
        {
            const word& patchName = fixedPatches[fpI];

            const label patchID = mesh.boundaryMesh().findPatchID(patchName);

            if (patchID == -1)
            {
                FatalErrorIn
                (
                    "void implicitVolumeSmoother::mapVolToPoint\n"
                    "(\n"
                    "    const volVectorField& cellMotionD,\n"
                    "    pointVectorField& pointMotionD\n"
                    ")"
                )   << "Fixed patch " << patchName << " does not exist!"
                    << abort(FatalError);
            }

            // Set pointMotion to zero on the patch

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchID].meshPoints();

            forAll(meshPoints, pI)
            {
                pointMotionDI[meshPoints[pI]] = vector::zero;
            }
        }
    }

    if (dict().found("flatPatches"))
    {
        if (debug)
        {
            Info<< "    Correcting flat patches" << endl;
        }

        const wordList flatPatches
        (
            dict().lookup("flatPatches")
        );

        forAll(flatPatches, fpI)
        {
            const word& patchName = flatPatches[fpI];

            const label patchID = mesh.boundaryMesh().findPatchID(patchName);

            if (patchID == -1)
            {
                FatalErrorIn
                (
                    "void implicitVolumeSmoother::mapVolToPoint\n"
                    "(\n"
                    "    const volVectorField& cellMotionD,\n"
                    "    pointVectorField& pointMotionD\n"
                    ")"
                )   << "Flat patch " << patchName << " does not exist!"
                    << abort(FatalError);
            }

            // Calculate the average patch normal
            const vector avgN =
                gAverage(mesh.boundaryMesh()[patchID].faceNormals());

            // Remove the patch normal component of the point motion

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchID].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                pointMotionDI[pointID] =
                    ((I - sqr(avgN)) & pointMotionDI[pointID]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSmoother::meshSmoother
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    oldInstance_("constant"),
    oldPoints_(mesh.allPoints()),
    fieldsAdvected_(false)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<meshSmoother> meshSmoother::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    // Lookup the type of mesh smoother from the dict
    const word smootherType = word(dict.lookup("type"));

    Info<< "    Mesh Smoother: " << smootherType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(smootherType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "meshSmoother::New(fvMesh& mesh, const dictionary& dict)"
        )   << "Unknown meshSmoother type "
            << smootherType << endl << endl
            << "Valid  meshSmoothers are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<meshSmoother>
    (
        cstrIter()(mesh, dict.subDict(smootherType + "Coeffs"))
    );
}


void meshSmoother::storeOldPoints(const pointField& oldPoints)
{
    oldPoints_ = oldPoints;
}


// void meshSmoother::mapFields()
// {
//     if (fieldsAdvected_)
//     {
//         FatalErrorIn("void meshSmoother::mapFields()")
//             << "Either mapFields OR advect fields can be called but not both!"
//             << abort(FatalError);
//     }

//     Info<< nl << "Performing conservative mapping of the fields" << nl
//         << "    internal field: 2nd order" << nl
//         << "    boundary field: 1st order"
//         << nl << endl;

//     // We assume that the current mesh is in the new position, so we will make
//     // the original mesh using the previously stored oldPoints
//     // This function maps the geometric fields from the previous mesh to the new
//     // mesh, where only point motion has occurred i.e. not topological changes

//     // Make a copy of the current mesh components as they will be transferred to
//     // the mesh
//     pointField pointsCopy = oldPoints_;
//     faceList facesCopy = mesh_.faces();
//     labelList allOwnerCopy = mesh_.faceOwner();
//     labelList allNeighbourCopy = mesh_.faceNeighbour();

//     // Create the old mesh
//     fvMesh oldMesh
//     (
//         IOobject
//         (
//             "oldMesh",
//             mesh_.time().timeName(),
//             mesh_.time(),
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         xferMove(pointsCopy),
//         xferMove(facesCopy),
//         xferMove(allOwnerCopy),
//         xferMove(allNeighbourCopy)
//     );

//     // Add the boundary patches by copy the current mesh boundary
//     List<polyPatch*> meshBoundary(mesh_.boundaryMesh().size());
//     forAll(mesh_.boundaryMesh(), patchI)
//     {
//         meshBoundary[patchI] =
//             mesh_.boundaryMesh()[patchI].clone
//             (
//                 oldMesh.boundaryMesh(),
//                 patchI,
//                 mesh_.boundaryMesh()[patchI].size(),
//                 mesh_.boundaryMesh()[patchI].start()
//             ).ptr();
//     }
//     oldMesh.addFvPatches(meshBoundary);

//     // Create the interpolator from the old mesh to the new mesh
//     conservativeMeshToMesh meshToMeshInterp
//     (
//         oldMesh,              // meshSource
//         mesh_,                // meshTarget
//         1,                    // nThreads: is this OpenMP?
//         false,                // forceRecalc
//         false,                // writeAddr
//         true                  // rescale weights for inconsistent cells
//     );

//     // Map volFields
//     // Mapping options:
//     // conservativeMeshToMesh::CONSERVATIVE
//     // conservativeMeshToMesh::CONSERVATIVE_FIRST_ORDER
//     // conservativeMeshToMesh::INVERSE_DISTANCE

//     MapConservativeVolFields<scalar>
//     (
//         meshToMeshInterp,
//         conservativeMeshToMesh::CONSERVATIVE
//     );

//     MapConservativeVolFields<vector>
//     (
//         meshToMeshInterp,
//         conservativeMeshToMesh::CONSERVATIVE
//     );

//     MapConservativeVolFields<tensor>
//     (
//         meshToMeshInterp,
//         conservativeMeshToMesh::CONSERVATIVE
//     );

//     MapConservativeVolFields<symmTensor>
//     (
//         meshToMeshInterp,
//         conservativeMeshToMesh::CONSERVATIVE
//     );

//     MapConservativeVolFields<sphericalTensor>
//     (
//         meshToMeshInterp,
//         conservativeMeshToMesh::CONSERVATIVE
//     );
// }


void meshSmoother::advectFields
(
    const scalarField& sweptVol,
    const scalarField& volOld,
    const scalarField& volNew
)
{
    Info<< nl << "Advecting the fields" << nl << endl;

    fieldsAdvected_ = true;

    // Create sweptVol surfaceScalarField
    surfaceScalarField sweptVolField
    (
        IOobject
        (
            "sweptVolField",
            mesh().time().timeName(),
            mesh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVolume, 0.0)
    );

    // Populate the internal field
    sweptVolField.internalField() =
        scalarField::subField(sweptVol, mesh().nInternalFaces());

    // Populate the boundary
    const fvPatchList& patches = mesh().boundary();

    forAll(patches, patchI)
    {
        sweptVolField.boundaryField()[patchI] =
            patches[patchI].patchSlice(sweptVol);
    }

    // Create volFields for the cell volumes fields
    volScalarField volOldField
    (
        IOobject
        (
            "oldVolField",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVolume, 0.0),
        "zeroGradient"
    );

    // Set the internal field
    volOldField.internalField() = volOld;

    // Update the boundaries
    volOldField.correctBoundaryConditions();

    // Create volFields for the cell volumes fields
    volScalarField volNewField
    (
        IOobject
        (
            "newVolField",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimVolume, 0.0),
        "zeroGradient"
    );

    // Set the internal field
    volNewField.internalField() = volNew;

    // Update the boundaries
    volNewField.correctBoundaryConditions();

    // List of fiels that should not be advected
    wordList ignoreFields(17);
    ignoreFields[0] = volOldField.name();
    ignoreFields[1] = volNewField.name();
    ignoreFields[2] = "mu";
    ignoreFields[3] = "lambda";
    ignoreFields[4] = "twoMuLambda";
    ignoreFields[5] = "materials";
    ignoreFields[6] = "curMaterial_0";
    ignoreFields[7] = "curMaterial_1";
    ignoreFields[8] = "curMaterial_2";
    ignoreFields[9] = "curMaterial_3";
    ignoreFields[10] = "curMaterial_4";
    ignoreFields[11] = "curMaterial_5";
    ignoreFields[12] = "curMaterial_6";
    ignoreFields[13] = "curMaterial_7";
    ignoreFields[14] = "curMaterial_8";
    ignoreFields[15] = "curMaterial_9";
    ignoreFields[16] = "curMaterial_10";
    Info<< "Currently the following fields are not advected: " << ignoreFields
        << endl;

    // Advect the fields
    AdvectVolFields<scalar>
    (
        sweptVolField, volOldField, volNewField, ignoreFields
    );
    AdvectVolFields<vector>
    (
        sweptVolField, volOldField, volNewField, ignoreFields
    );
    AdvectVolFields<tensor>
    (
        sweptVolField, volOldField, volNewField, ignoreFields
    );
    AdvectVolFields<symmTensor>
    (
        sweptVolField, volOldField, volNewField, ignoreFields
    );
    AdvectVolFields<sphericalTensor>
    (
        sweptVolField, volOldField, volNewField, ignoreFields
    );
}


void meshSmoother::writeMesh()
{
    Info<< nl << "Writing the mesh" << nl << endl;

    if (overwrite())
    {
        Info<< "    Overwriting the mesh in " << oldInstance() << endl;
        mesh().setInstance(oldInstance());
    }

    mesh().write();
}

} // End namespace Foam

// ************************************************************************* //
