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

Application
    simpleSmoothMesh

Description
    Smoothing mesh based on smoothMeshDict in system directory.
    Fields are also mapped.

    Example simpleSmoothMeshDict:

    // Smoothing method
    type                  implicitVolume;

    // Smoothing settings
    implicitVolumeCoeffs
    {
        // Optional: cell zone(s) to smooth
        // If this entry is not found then the entire mesh is smoothed
        cellZones             (billet);

        // Number of outer smoothing iterations
        nCorrectors           10;

        // Overwrite option
        // no: writes to the next time-step
        // yes: overwrites the mesh
        overwrite             no;

        // Optional: if the angle between two faces is greater than this
        // feature angle then this sharp edge is preserved
        featureAngle          45;

        // Optional: points on fixed patches are not moved during smoothing
        fixedPatches ("top" left");
    }


Author
    Philip Cardiff UCD
    Peter De Jaeger
    Zeljko Tukovic

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "meshSmoother.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * * //

template<class Type>
void ReadVolFields
(
    PtrList< GeometricField<Type, fvPatchField, volMesh> >& curFields,
    const IOobjectList& objects,
    const fvMesh& mesh
)
{
    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    // Lookup fields of type fieldClassName
    IOobjectList curFieldNames = objects.lookupClass(fieldClassName);

    label curFieldI = 0;

    curFields.setSize(curFieldNames.size());

    for
    (
        IOobjectList::iterator fieldIter = curFieldNames.begin();
        fieldIter != curFieldNames.end();
        ++fieldIter
    )
    {
        IOobject fieldIOobjectI
        (
            fieldIter()->name(),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldIOobjectI.headerOk())
        {
            Info<< "Reading field " << fieldIter()->name() << endl;

            // Read field and register to the mesh
            curFields.set
            (
                curFieldI,
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    *fieldIter(),
                    mesh
                )
            );

            curFieldI++;
        }
    }

    curFields.resize(curFieldI);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("latestTime", "");
    argList::validOptions.insert("noFields", "");
    argList::noParallel();

#   include "setRootCase.H"
#   include "createTime.H"

    runTime.functionObjects().off();

    if (args.optionFound("latestTime"))
    {
        instantList Times = runTime.times();
        runTime.setTime(Times[Times.size() - 1], Times.size() - 1);
        Info<< "Reading mesh from time: " << runTime.value() << endl;
    }

    const bool noFields = bool(args.optionFound("noFields"));

#   include "createMesh.H"

    // Read the smoothMeshDict
    IOdictionary dict
    (
        IOobject
        (
            "simpleSmoothMeshDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create the smoother
    autoPtr<meshSmoother> smoother = meshSmoother::New(mesh, dict);

    // Clear out previous addressing
    smoother().clearOut();

    // Store the old points before smoothing
    smoother().storeOldPoints(mesh.points());

    // Map the fields
    // Lists to hold the fields
    PtrList<volScalarField> scaFields;
    PtrList<volVectorField> vecFields;
    PtrList<volTensorField> tenFields;
    PtrList<volSymmTensorField> symmTenFields;
    PtrList<volDiagTensorField> diaTenFields;
    if (!noFields)
    {
        // Search for fields
        Info<< nl << "Searching for fields in time = " << mesh.time().timeName()
            << endl;
        IOobjectList objects(mesh, mesh.time().timeName(), mesh.local());

        // Populate the lists
        ReadVolFields<scalar>(scaFields, objects, mesh);
        ReadVolFields<vector>(vecFields, objects, mesh);
        ReadVolFields<tensor>(tenFields, objects, mesh);
        ReadVolFields<symmTensor>(symmTenFields, objects, mesh);
        ReadVolFields<diagTensor>(diaTenFields, objects, mesh);

        Info<< nl << "Note: surfaceFields and pointFields are currently not "
            << "mapped!" << endl;
    }

    // Perform the smoothing
    smoother().smooth();

    // Map the fields if they have not been advected
    if (!smoother().fieldsAdvected() && !noFields)
    {
        smoother().mapFields();
    }

    // Write the mesh to disk
    runTime++;
    smoother().writeMesh();

    // Force fields to write
    if (!noFields)
    {
        // Write the fields
        Info<< nl << "Writing fields to time = " << mesh.time().timeName()
            << endl;
        runTime.write();

        forAll(scaFields, fI)
        {
            scaFields[fI].write();
        }

        forAll(vecFields, fI)
        {
            vecFields[fI].write();
        }

        forAll(tenFields, fI)
        {
            tenFields[fI].write();
        }

        forAll(symmTenFields, fI)
        {
            symmTenFields[fI].write();
        }

        forAll(diaTenFields, fI)
        {
            diaTenFields[fI].write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
