/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Application
    decomposePar

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of FOAM.

Usage

    - decomposePar [OPTION]

    @param -cellDist \n
    Write the cell distribution as a labelList for use with 'manual'
    decomposition method and as a volScalarField for post-processing.

    @param -region regionName \n
    Decompose named region. Does not check for existence of processor*.

    @param -copyUniform \n
    Copy any @a uniform directories too.

    @param -fields \n
    Use existing geometry decomposition and convert fields only.

    @param -filterPatches \n
    Remove empty patches when decomposing the geometry.

    @param -force \n
    Remove any existing @a processor subdirectories before decomposing the
    geometry.

    @param -ifRequired \n
    Only decompose the geometry if the number of domains has changed from a
    previous decomposition. No @a processor subdirectories will be removed
    unless the @a -force option is also specified. This option can be used
    to avoid redundant geometry decomposition (eg, in scripts), but should
    be used with caution when the underlying (serial) geometry or the
    decomposition method etc. have been changed between decompositions.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "IOobjectList.H"
#include "processorFvPatchFields.H"
#include "domainDecomposition.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"

#include "tetPointFields.H"
#include "elementFields.H"
#include "tetFemMatrices.H"
#include "tetPointFieldDecomposer.H"

#include "pointFields.H"

#include "readFields.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "lagrangianFieldDecomposer.H"

#include "faCFD.H"
#include "emptyFaPatch.H"
#include "faMeshDecomposition.H"
#include "faFieldDecomposer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"
    argList::validOptions.insert("cellDist", "");
    argList::validOptions.insert("copyUniform", "");
    argList::validOptions.insert("fields", "");
    argList::validOptions.insert("filterPatches", "");
    argList::validOptions.insert("force", "");
    argList::validOptions.insert("ifRequired", "");

#   include "setRootCase.H"

    word regionName = fvMesh::defaultRegion;
    word regionDir = word::null;

    if (args.optionFound("region"))
    {
        regionName = args.option("region");
        regionDir = regionName;
        Info<< "Decomposing mesh " << regionName << nl << endl;
    }


    bool writeCellDist = args.optionFound("cellDist");
    bool copyUniform = args.optionFound("copyUniform");
    bool decomposeFieldsOnly = args.optionFound("fields");
    bool filterPatches = args.optionFound("filterPatches");
    bool forceOverwrite = args.optionFound("force");
    bool ifRequiredDecomposition = args.optionFound("ifRequired");

#   include "createTime.H"

    Info<< "Time = " << runTime.timeName() << endl;

    // Determine the existing processor count directly
    label nProcs = 0;
    while
    (
        isDir
        (
            runTime.path()
           /(word("processor") + name(nProcs))
           /runTime.constant()
           /regionDir
           /polyMesh::meshSubDir
        )
    )
    {
        ++nProcs;
    }

    // get requested numberOfSubdomains
    label nDomains = 0;
    {
        IOdictionary decompDict
        (
            IOobject
            (
                "decomposeParDict",
                runTime.time().system(),
                regionDir,          // use region if non-standard
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        decompDict.lookup("numberOfSubdomains") >> nDomains;
    }

    if (decomposeFieldsOnly)
    {
        // Sanity check on previously decomposed case
        if (nProcs != nDomains)
        {
            FatalErrorIn(args.executable())
                << "Specified -fields, but the case was decomposed with "
                << nProcs << " domains"
                << nl
                << "instead of " << nDomains
                << " domains as specified in decomposeParDict"
                << nl
                << exit(FatalError);
        }
    }
    else if (nProcs)
    {
        bool procDirsProblem = true;

        if (regionName != fvMesh::defaultRegion)
        {
            decomposeFieldsOnly = false;
            procDirsProblem = false;
        }


        if (ifRequiredDecomposition && nProcs == nDomains)
        {
            // we can reuse the decomposition
            decomposeFieldsOnly = true;
            procDirsProblem = false;
            forceOverwrite = false;

            Info<< "Using existing processor directories" << nl;
        }

        if (forceOverwrite)
        {
            Info<< "Removing " << nProcs
                << " existing processor directories" << endl;

            // remove existing processor dirs
            // reverse order to avoid gaps if someone interrupts the process
            for (label procI = nProcs-1; procI >= 0; --procI)
            {
                fileName procDir
                (
                    runTime.path()/(word("processor") + name(procI))
                );

                rmDir(procDir);
            }

            procDirsProblem = false;
        }

        if (procDirsProblem)
        {
            FatalErrorIn(args.executable())
                << "Case is already decomposed with " << nProcs
                << " domains, use the -force option or manually" << nl
                << "remove processor directories before decomposing. e.g.,"
                << nl
                << "    rm -rf " << runTime.path().c_str() << "/processor*"
                << nl
                << exit(FatalError);
        }
    }

    Info<< "Create mesh for region " << regionName << endl;
    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    domainDecomposition meshDecomp
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );


    // Decompose the mesh
    if (!decomposeFieldsOnly)
    {
        meshDecomp.decomposeMesh(filterPatches);

        meshDecomp.writeDecomposition();

        if (writeCellDist)
        {
            const labelList& procIds = meshDecomp.cellToProc();

            // Write the decomposition as labelList for use with 'manual'
            // decomposition method.
            labelIOList cellDecomposition
            (
                IOobject
                (
                    "cellDecomposition",
                    mesh.facesInstance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procIds
            );
            cellDecomposition.write();

            Info<< nl << "Wrote decomposition to "
                << cellDecomposition.objectPath()
                << " for use in manual decomposition." << endl;

            // Write as volScalarField for post-processing
            volScalarField cellDist
            (
                IOobject
                (
                    "cellDist",
                    runTime.timeName(),
                    mesh.dbDir(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("cellDist", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            );

            forAll(procIds, celli)
            {
               cellDist[celli] = procIds[celli];
            }

            cellDist.write();

            Info<< nl << "Wrote decomposition as volScalarField to "
                << cellDist.name() << " for use in post-processing."
                << endl;
        }
    }


    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    // Construct the vol fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    PtrList<volScalarField> volScalarFields;
    readFields(mesh, objects, volScalarFields);

    PtrList<volVectorField> volVectorFields;
    readFields(mesh, objects, volVectorFields);

    PtrList<volSphericalTensorField> volSphericalTensorFields;
    readFields(mesh, objects, volSphericalTensorFields);

    PtrList<volSymmTensorField> volSymmTensorFields;
    readFields(mesh, objects, volSymmTensorFields);

    PtrList<volTensorField> volTensorFields;
    readFields(mesh, objects, volTensorFields);


    // Construct the surface fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PtrList<surfaceScalarField> surfaceScalarFields;
    readFields(mesh, objects, surfaceScalarFields);
    PtrList<surfaceVectorField> surfaceVectorFields;
    readFields(mesh, objects, surfaceVectorFields);
    PtrList<surfaceSphericalTensorField> surfaceSphericalTensorFields;
    readFields(mesh, objects, surfaceSphericalTensorFields);
    PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
    readFields(mesh, objects, surfaceSymmTensorFields);
    PtrList<surfaceTensorField> surfaceTensorFields;
    readFields(mesh, objects, surfaceTensorFields);


    // Construct the point fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Fix ZT/PC, 29-Aug-18
    const pointMesh& pMesh = pointMesh::New(mesh, true);

    PtrList<pointScalarField> pointScalarFields;
    readFields(pMesh, objects, pointScalarFields);

    PtrList<pointVectorField> pointVectorFields;
    readFields(pMesh, objects, pointVectorFields);

    PtrList<pointSphericalTensorField> pointSphericalTensorFields;
    readFields(pMesh, objects, pointSphericalTensorFields);

    PtrList<pointSymmTensorField> pointSymmTensorFields;
    readFields(pMesh, objects, pointSymmTensorFields);

    PtrList<pointTensorField> pointTensorFields;
    readFields(pMesh, objects, pointTensorFields);


    // Construct the tetPoint fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    tetPolyMesh* tetMeshPtr = NULL;

    PtrList<tetPointScalarField> tetPointScalarFields;
    PtrList<tetPointVectorField> tetPointVectorFields;
    PtrList<tetPointSphericalTensorField> tetPointSphericalTensorFields;
    PtrList<tetPointSymmTensorField> tetPointSymmTensorFields;
    PtrList<tetPointTensorField> tetPointTensorFields;

    PtrList<elementScalarField> elementScalarFields;
    PtrList<elementVectorField> elementVectorFields;

    if
    (
        objects.lookupClass("tetPointScalarField").size() > 0
     || objects.lookupClass("tetPointVectorField").size() > 0
     || objects.lookupClass("tetPointSphericalTensorField").size() > 0
     || objects.lookupClass("tetPointSymmTensorField").size() > 0
     || objects.lookupClass("tetPointTensorField").size() > 0

     || objects.lookupClass("elementScalarField").size() > 0
     || objects.lookupClass("elementVectorField").size() > 0
    )
    {
        tetMeshPtr = new tetPolyMesh(mesh);
        tetPolyMesh& tetMesh = *tetMeshPtr;

        readFields(tetMesh, objects, tetPointScalarFields);
        readFields(tetMesh, objects, tetPointVectorFields);
        readFields(tetMesh, objects, tetPointSphericalTensorFields);
        readFields(tetMesh, objects, tetPointSymmTensorFields);
        readFields(tetMesh, objects, tetPointTensorFields);

        readFields(tetMesh, objects, elementScalarFields);
        readFields(tetMesh, objects, elementVectorFields);
    }


    // Construct the Lagrangian fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fileNameList cloudDirs
    (
        readDir(runTime.timePath()/cloud::prefix, fileName::DIRECTORY)
    );

    // Particles
    PtrList<Cloud<indexedParticle> > lagrangianPositions(cloudDirs.size());
    // Particles per cell
    PtrList< List<SLList<indexedParticle*>*> > cellParticles(cloudDirs.size());

    PtrList<PtrList<labelIOField> > lagrangianLabelFields(cloudDirs.size());
    PtrList<PtrList<scalarIOField> > lagrangianScalarFields(cloudDirs.size());
    PtrList<PtrList<vectorIOField> > lagrangianVectorFields(cloudDirs.size());
    PtrList<PtrList<sphericalTensorIOField> > lagrangianSphericalTensorFields
    (
        cloudDirs.size()
    );
    PtrList<PtrList<symmTensorIOField> > lagrangianSymmTensorFields
    (
        cloudDirs.size()
    );
    PtrList<PtrList<tensorIOField> > lagrangianTensorFields
    (
        cloudDirs.size()
    );

    label cloudI = 0;

    forAll(cloudDirs, i)
    {
        IOobjectList sprayObjs
        (
            mesh,
            runTime.timeName(),
            cloud::prefix/cloudDirs[i]
        );

        IOobject* positionsPtr = sprayObjs.lookup("positions");

        if (positionsPtr)
        {
            // Read lagrangian particles
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            Info<< "Identified lagrangian data set: " << cloudDirs[i] << endl;

            lagrangianPositions.set
            (
                cloudI,
                new Cloud<indexedParticle>
                (
                    mesh,
                    cloudDirs[i],
                    false
                )
            );


            // Sort particles per cell
            // ~~~~~~~~~~~~~~~~~~~~~~~

            cellParticles.set
            (
                cloudI,
                new List<SLList<indexedParticle*>*>
                (
                    mesh.nCells(),
                    static_cast<SLList<indexedParticle*>*>(NULL)
                )
            );

            label i = 0;

            forAllIter
            (
                Cloud<indexedParticle>,
                lagrangianPositions[cloudI],
                iter
            )
            {
                iter().index() = i++;

                label celli = iter().cell();

                // Check
                if (celli < 0 || celli >= mesh.nCells())
                {
                    FatalErrorIn(args.executable())
                        << "Illegal cell number " << celli
                        << " for particle with index " << iter().index()
                        << " at position " << iter().position() << nl
                        << "Cell number should be between 0 and "
                        << mesh.nCells()-1 << nl
                        << "On this mesh the particle should be in cell "
                        << mesh.findCell(iter().position())
                        << exit(FatalError);
                }

                if (!cellParticles[cloudI][celli])
                {
                    cellParticles[cloudI][celli] =
                        new SLList<indexedParticle*>();
                }

                cellParticles[cloudI][celli]->append(&iter());
            }

            // Read fields
            // ~~~~~~~~~~~

            IOobjectList lagrangianObjects
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudDirs[cloudI]
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianLabelFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianScalarFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianVectorFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianSphericalTensorFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianSymmTensorFields
            );

            lagrangianFieldDecomposer::readFields
            (
                cloudI,
                lagrangianObjects,
                lagrangianTensorFields
            );

            cloudI++;
        }
    }

    lagrangianPositions.setSize(cloudI);
    cellParticles.setSize(cloudI);
    lagrangianLabelFields.setSize(cloudI);
    lagrangianScalarFields.setSize(cloudI);
    lagrangianVectorFields.setSize(cloudI);
    lagrangianSphericalTensorFields.setSize(cloudI);
    lagrangianSymmTensorFields.setSize(cloudI);
    lagrangianTensorFields.setSize(cloudI);


    // Any uniform data to copy/link?
    fileName uniformDir("uniform");

    if (isDir(runTime.timePath()/uniformDir))
    {
        Info<< "Detected additional non-decomposed files in "
            << runTime.timePath()/uniformDir
            << endl;
    }
    else
    {
        uniformDir.clear();
    }

    Info<< endl;

    // Split the fields over processors
    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        Info<< "Processor " << procI << ": field transfer" << endl;

        // open the database
        Time processorDb
        (
            Time::controlDictName,
            args.rootPath(),
            args.caseName()/fileName(word("processor") + name(procI))
        );

        processorDb.setTime(runTime);

        // Remove files remnants that can cause horrible problems
        // - mut and nut are used to mark the new turbulence models,
        //   their existence prevents old models from being upgraded
        // 1.6.x merge.  HJ, 25/Aug/2010
        {
            fileName timeDir(processorDb.path()/processorDb.timeName());

            rm(timeDir/"mut");
            rm(timeDir/"nut");
            rm(timeDir/"mut.gz");
            rm(timeDir/"nut.gz");
        }

        // Read the mesh
        fvMesh procMesh
        (
            IOobject
            (
                regionName,
                processorDb.timeName(),
                processorDb
            )
        );

        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // FV fields
        if
        (
            volScalarFields.size()
         || volVectorFields.size()
         || volSphericalTensorFields.size()
         || volSymmTensorFields.size()
         || volTensorFields.size()
         || surfaceScalarFields.size()
         || surfaceVectorFields.size()
         || surfaceSphericalTensorFields.size()
         || surfaceSymmTensorFields.size()
         || surfaceTensorFields.size()
        )
        {
            fvFieldDecomposer fieldDecomposer
            (
                mesh,
                procMesh,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            fieldDecomposer.decomposeFields(volScalarFields);
            fieldDecomposer.decomposeFields(volVectorFields);
            fieldDecomposer.decomposeFields(volSphericalTensorFields);
            fieldDecomposer.decomposeFields(volSymmTensorFields);
            fieldDecomposer.decomposeFields(volTensorFields);

            fieldDecomposer.decomposeFields(surfaceScalarFields);
            fieldDecomposer.decomposeFields(surfaceVectorFields);
            fieldDecomposer.decomposeFields(surfaceSphericalTensorFields);
            fieldDecomposer.decomposeFields(surfaceSymmTensorFields);
            fieldDecomposer.decomposeFields(surfaceTensorFields);
        }


        // Point fields
        if
        (
            pointScalarFields.size()
         || pointVectorFields.size()
         || pointSphericalTensorFields.size()
         || pointSymmTensorFields.size()
         || pointTensorFields.size()
        )
        {
            // Fix ZT/PC, 29-Aug-18
            const pointMesh& procPMesh = pointMesh::New(procMesh, true);

            pointFieldDecomposer fieldDecomposer
            (
                pMesh,
                procPMesh,
                pointProcAddressing,
                boundaryProcAddressing
            );

            fieldDecomposer.decomposeFields(pointScalarFields);
            fieldDecomposer.decomposeFields(pointVectorFields);
            fieldDecomposer.decomposeFields(pointSphericalTensorFields);
            fieldDecomposer.decomposeFields(pointSymmTensorFields);
            fieldDecomposer.decomposeFields(pointTensorFields);
        }


        // tetPoint fields
        if (tetMeshPtr)
        {
            const tetPolyMesh& tetMesh = *tetMeshPtr;
            tetPolyMesh procTetMesh(procMesh);

            tetPointFieldDecomposer fieldDecomposer
            (
                tetMesh,
                procTetMesh,
                pointProcAddressing,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            fieldDecomposer.decomposeFields(tetPointScalarFields);
            fieldDecomposer.decomposeFields(tetPointVectorFields);
            fieldDecomposer.decomposeFields(tetPointSphericalTensorFields);
            fieldDecomposer.decomposeFields(tetPointSymmTensorFields);
            fieldDecomposer.decomposeFields(tetPointTensorFields);

            fieldDecomposer.decomposeFields(elementScalarFields);
            fieldDecomposer.decomposeFields(elementVectorFields);
        }


        // If there is lagrangian data write it out
        forAll(lagrangianPositions, cloudI)
        {
            if (lagrangianPositions[cloudI].size())
            {
                lagrangianFieldDecomposer fieldDecomposer
                (
                    mesh,
                    procMesh,
                    cellProcAddressing,
                    cloudDirs[cloudI],
                    lagrangianPositions[cloudI],
                    cellParticles[cloudI]
                );

                // Lagrangian fields
                if
                (
                    lagrangianLabelFields[cloudI].size()
                 || lagrangianScalarFields[cloudI].size()
                 || lagrangianVectorFields[cloudI].size()
                 || lagrangianSphericalTensorFields[cloudI].size()
                 || lagrangianSymmTensorFields[cloudI].size()
                 || lagrangianTensorFields[cloudI].size()
                )
                {
                    fieldDecomposer.decomposeFields
                    (
                        cloudDirs[cloudI],
                        lagrangianLabelFields[cloudI]
                    );
                    fieldDecomposer.decomposeFields
                    (
                        cloudDirs[cloudI],
                        lagrangianScalarFields[cloudI]
                    );
                    fieldDecomposer.decomposeFields
                    (
                        cloudDirs[cloudI],
                        lagrangianVectorFields[cloudI]
                    );
                    fieldDecomposer.decomposeFields
                    (
                        cloudDirs[cloudI],
                        lagrangianSphericalTensorFields[cloudI]
                    );
                    fieldDecomposer.decomposeFields
                    (
                        cloudDirs[cloudI],
                        lagrangianSymmTensorFields[cloudI]
                    );
                    fieldDecomposer.decomposeFields
                    (
                        cloudDirs[cloudI],
                        lagrangianTensorFields[cloudI]
                    );
                }
            }
        }


        // Any non-decomposed data to copy?
        if (uniformDir.size())
        {
            const fileName timePath = processorDb.timePath();

            if (copyUniform || meshDecomp.distributed())
            {
                cp
                (
                    runTime.timePath()/uniformDir,
                    timePath/uniformDir
                );
            }
            else
            {
                // Link with relative paths
                const string parentPath = string("..")/"..";

                fileName currentDir(cwd());
                chDir(timePath);
                if (!exists(uniformDir))
                {
                    ln
                    (
                        parentPath/runTime.timeName()/uniformDir,
                        uniformDir
                    );
                    chDir(currentDir);
                }
            }
        }
    }


    if (tetMeshPtr)
    {
        delete tetMeshPtr;
        tetMeshPtr = NULL;
    }


    // Finite area mesh and field decomposition

    IOobject faMeshBoundaryIOobj
    (
        "faBoundary",
        mesh.time().findInstance
        (
            mesh.meshDir(),
            "boundary"
        ),
        faMesh::meshSubDir,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );


    if(faMeshBoundaryIOobj.headerOk())
    {
        Info << "\nFinite area mesh decomposition" << endl;

        faMeshDecomposition aMesh(mesh);

        aMesh.decomposeMesh(filterPatches);

        aMesh.writeDecomposition();


        // Construct the area fields
        // ~~~~~~~~~~~~~~~~~~~~~~~~
        PtrList<areaScalarField> areaScalarFields;
        readFields(aMesh, objects, areaScalarFields);

        PtrList<areaVectorField> areaVectorFields;
        readFields(aMesh, objects, areaVectorFields);

        PtrList<areaSphericalTensorField> areaSphericalTensorFields;
        readFields(aMesh, objects, areaSphericalTensorFields);

        PtrList<areaSymmTensorField> areaSymmTensorFields;
        readFields(aMesh, objects, areaSymmTensorFields);

        PtrList<areaTensorField> areaTensorFields;
        readFields(aMesh, objects, areaTensorFields);


        // Construct the edge fields
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PtrList<edgeScalarField> edgeScalarFields;
        readFields(aMesh, objects, edgeScalarFields);

        Info << endl;

        // Split the fields over processors
        for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
        {
            Info<< "Processor " << procI
                << ": finite area field transfer" << endl;

            // open the database
            Time processorDb
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            );

            processorDb.setTime(runTime);

            // Read the mesh
            fvMesh procFvMesh
            (
                IOobject
                (
                    regionName,
                    processorDb.timeName(),
                    processorDb
                )
            );

            faMesh procMesh(procFvMesh);

            word faceLabelsInstance =
                procMesh.time().findInstance
                (
                    procMesh.meshDir(),
                    "faceLabels"
                );


            labelIOList faceProcAddressing
            (
                IOobject
                (
                    "faceProcAddressing",
                    faceLabelsInstance,
                    procMesh.meshSubDir,
                    procFvMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            labelIOList boundaryProcAddressing
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    faceLabelsInstance,
                    procMesh.meshSubDir,
                    procFvMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            // FA fields
            if
            (
                areaScalarFields.size()
             || areaVectorFields.size()
             || areaSphericalTensorFields.size()
             || areaSymmTensorFields.size()
             || areaTensorFields.size()
             || edgeScalarFields.size()
            )
            {
                labelIOList edgeProcAddressing
                (
                    IOobject
                    (
                        "edgeProcAddressing",
                        faceLabelsInstance,
                        procMesh.meshSubDir,
                        procFvMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                faFieldDecomposer fieldDecomposer
                (
                    aMesh,
                    procMesh,
                    edgeProcAddressing,
                    faceProcAddressing,
                    boundaryProcAddressing
                );

                fieldDecomposer.decomposeFields(areaScalarFields);
                fieldDecomposer.decomposeFields(areaVectorFields);
                fieldDecomposer.decomposeFields(areaSphericalTensorFields);
                fieldDecomposer.decomposeFields(areaSymmTensorFields);
                fieldDecomposer.decomposeFields(areaTensorFields);

                fieldDecomposer.decomposeFields(edgeScalarFields);
            }
        }
    }


    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
