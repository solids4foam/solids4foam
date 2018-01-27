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
    contactFields

Description
    Read the results from a solid mechanics case where there are solidContact
    boundar conditions, and writes out contact fields, such as contactPressure,
    contactArea and contactDistance, for visualisation.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solidContactFvPatchVectorField.H"
#include "pointMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("displacementField", "name");
    argList::validOptions.insert("updatedLagrangian", "");
    argList::validOptions.insert("totalLagrangian", "");

#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
    instantList Times = runTime.times();
#   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);
#   include "createMesh.H"

    // We will assume that the displacement field is called D, but a
    // command-line argument can change this to U (for backwards compatibility)
    word DName = "D";
    if (args.optionFound("displacementField"))
    {
        DName = word(args.optionLookup("displacementField")());
    }

    // Check if the case is nonlinear geometry (updated or total Lagrangian)
    bool totalLagrangian = false;
    bool updatedLagrangian = false;
    if
    (
        args.optionFound("updatedLagrangian")
     && args.optionFound("totalLagrangian")
    )
    {
        FatalError
            << "Invalid option combination: both updatedLagrangian and "
            << "totalLagrangian cannot be simultaneously specified!"
            << abort(FatalError);
    }
    else if (args.optionFound("updatedLagrangian"))
    {
        Info<< nl << "Nonlinear geometry: updated Lagrangian" << endl;
        updatedLagrangian = true;
    }
    else if (args.optionFound("totalLagrangian"))
    {
        Info<< nl << "Nonlinear geometry: total Lagrangian" << endl;
        totalLagrangian = true;
    }
    else
    {
        Info<< "Linear geometry" << endl;
    }

    // Loop thrugh all the time-steps
    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< nl << "Time = " << runTime.timeName() << endl;

        // Read the mesh if it has changed
        mesh.readUpdate();

        // Check if stress field is found
        IOobject sigmaHeader
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check if displacement field is found
        IOobject DHeader
        (
            DName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // If there is no stress or displacement fields then we skip to the
        // next time-step
        if (!sigmaHeader.headerOk() || !DHeader.headerOk())
        {
            Info<< "    skipping time-step: both the sigma field and D "
                << "field were not found"
                << endl;
            continue;
        }

        // Read the stress field
        const volSymmTensorField sigma(sigmaHeader, mesh);

        // Read the displacement field
        volVectorField D(DHeader, mesh);

        // Read the total deformation gradient field, if found
        autoPtr<const volTensorField> FPtr;
        if (totalLagrangian)
        {
            FPtr.set
            (
                new volTensorField
                (
                    IOobject
                    (
                        "F",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }

        // Read the relative deformation gradient field, if found
        autoPtr<const volTensorField> relFPtr;
        if (updatedLagrangian)
        {
            relFPtr.set
            (
                new volTensorField
                (
                    IOobject
                    (
                        "relF",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }

        // Create the contact pressure field
        volScalarField contactPressure
        (
            IOobject
            (
                "contactPressure",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimForce/dimArea, 0.0)
        );

        // Create the contact gap (distance to contact) point field
        pointMesh pMesh(mesh);
        pointScalarField contactGap
        (
            IOobject
            (
                "contactGap",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedScalar("zero", dimLength, 0.0)
        );
        scalarField& contactGapI = contactGap.internalField();

        // Calculate the contact pressure field for solidContact patches
        int nContactPatches = 0;
        forAll(D.boundaryField(), patchI)
        {
            if
            (
                D.boundaryField()[patchI].type()
             == solidContactFvPatchVectorField::typeName
            )
            {
                Info<< "    Calculating contactPressure on patch: "
                    << mesh.boundaryMesh()[patchI].name() << endl;

                nContactPatches++;

                // Patch unit normals
                const vectorField n = mesh.boundary()[patchI].nf();

                // Patch stress
                const symmTensorField& sigmaP = sigma.boundaryField()[patchI];

                // Contact patch
                solidContactFvPatchVectorField& contactPatch =
                    refCast<solidContactFvPatchVectorField>
                    (
                        D.boundaryField()[patchI]
                    );

                // Patch mesh points
                const labelList& meshPoints =
                    mesh.boundaryMesh()[patchI].meshPoints();

                if (totalLagrangian)
                {
                    // Total deformation gradient
                    const tensorField& FP = FPtr().boundaryField()[patchI];

                    // Inverse of the total deformation gradient
                    const tensorField FPinv = inv(FP);

                    // Jacobian of the total deformation gradient
                    const scalarField JP = det(FP);

                    // Patch unit normals in the deformed configuration
                    vectorField curN = JP*FPinv.T() & n;
                    curN /= mag(curN);

                    // The contact pressure is the normal component of stress of
                    // traction on the boundary, where positive pressure is
                    // assumed to act inwards on the surface (hence the negative
                    // sign)
                    contactPressure.boundaryField()[patchI] =
                        -curN & (curN & sigmaP);
                }
                else if (updatedLagrangian)
                {
                    // Relative deformation gradient
                    const tensorField& relFP =
                        relFPtr().boundaryField()[patchI];

                    // Inverse of the relative deformation gradient
                    const tensorField relFPinv = inv(relFP);

                    // Jacobian of the relative deformation gradient
                    const scalarField relJP = det(relFP);

                    // Patch unit normals in the deformed configuration
                    vectorField curN = relJP*relFPinv.T() & n;
                    curN /= mag(curN);

                    // The contact pressure is the normal component of stress of
                    // traction on the boundary, where positive pressure is
                    // assumed to act inwards on the surface (hence the negative
                    // sign)
                    contactPressure.boundaryField()[patchI] =
                        -curN & (curN & sigmaP);
                }
                else
                {
                    // The contact pressure is the normal component of stress of
                    // traction on the boundary, where positive pressure is
                    // assumed to act inwards on the surface (hence the negative
                    // sign)
                    contactPressure.boundaryField()[patchI] = -n & (n & sigmaP);
                }

                // Before calculating the contact gap, we must move the contact
                // zones to the deformed position
                contactPatch.moveZonesToDeformedConfiguration();
                contactPatch.zoneToZone().movePoints
                (
                    tensorField(0), tensorField(0), vectorField(0)
                );

                // Set the contact point gap
                if (contactPatch.master())
                {
                    const scalarField& patchGap =
                        contactPatch.zoneToZone()
                        .masterPointDistanceToIntersection();

                    forAll(patchGap, pI)
                    {
                        const label pointID = meshPoints[pI];
                        contactGapI[pointID] = patchGap[pI];
                    }
                }
                else
                {
                    const scalarField& patchGap =
                        contactPatch.zoneToZone()
                       .slavePointDistanceToIntersection();

                    forAll(patchGap, pI)
                    {
                        const label pointID = meshPoints[pI];
                        contactGapI[pointID] = patchGap[pI];
                    }
                }
            }

            // Check if a least one contact patch was found
            if (nContactPatches == 0)
            {
                FatalError
                    << "The displacement field " << DName
                    << " has no solidContact patches!" << abort(FatalError);
            }
        }

        Info<< "    Writing fields to time = "
            << runTime.timeName() << endl;

        Info<< "    contactPressure volScalarField" << endl;
        contactPressure.write();

        Info<< "    contactGap pointScalarField" << endl;
        contactGap.write();
    }

    Info<< nl << "End" << endl;

    return(0);
}


// ************************************************************************* //
