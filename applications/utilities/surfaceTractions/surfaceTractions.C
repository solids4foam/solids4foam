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
    surfaceTractions

Description
    Calculates and writes the surface tractions as a volVectorField, using
    the sigma volSymmTensorField.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        // Check if sigma field is found
        IOobject sigmaHeader
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check if the deformation gradient field is found
        // This lets us know if it is a nonlinear geometry case
        IOobject FHeader
        (
            "F",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (sigmaHeader.headerOk() && !FHeader.headerOk())
        {
            Info<< "    Detected a linear geometry case" << endl;

            const volSymmTensorField sigma(sigmaHeader, mesh);

            volVectorField totalTraction
            (
                IOobject
                (
                    "totalTraction",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimForce/dimArea, vector::zero)
            );

            volScalarField normalTraction
            (
                IOobject
                (
                    "normalTraction",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimForce/dimArea, 0.0)
            );

            volVectorField shearTraction
            (
                IOobject
                (
                    "shearTraction",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimForce/dimArea, vector::zero)
            );

            forAll(totalTraction.boundaryField(), patchI)
            {
                const vectorField n = mesh.boundary()[patchI].nf();
                const symmTensorField& sigmaP =
                    sigma.boundaryField()[patchI];

                totalTraction.boundaryField()[patchI] = n & sigmaP;

                normalTraction.boundaryField()[patchI] =
                    n & totalTraction.boundaryField()[patchI];

                shearTraction.boundaryField()[patchI] =
                    (I - sqr(n)) & totalTraction.boundaryField()[patchI];
            }

            Info<< "    Writing totalTraction field" << endl;
            totalTraction.write();

            Info<< "    Writing normalTraction field" << endl;
            normalTraction.write();

            Info<< "    Writing shearTraction field" << endl;
            shearTraction.write();
        }
        else if (sigmaHeader.headerOk() && FHeader.headerOk())
        {
            // Check if the increment of displacement field is found
            IOobject DDHeader
            (
                "DD",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (DDHeader.headerOk())
            {
                Info<< "    Detected an updated Lagrangian case" << endl;

                const volSymmTensorField sigma(sigmaHeader, mesh);

                const volTensorField relF
                (
                    IOobject
                    (
                        "relF",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                volVectorField totalTraction
                (
                    IOobject
                    (
                        "totalTraction",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedVector("zero", dimForce/dimArea, vector::zero)
                );

                volScalarField normalTraction
                (
                    IOobject
                    (
                        "normalTraction",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zero", dimForce/dimArea, 0.0)
                );

                volVectorField shearTraction
                (
                    IOobject
                    (
                        "shearTraction",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedVector("zero", dimForce/dimArea, vector::zero)
                );

                forAll(totalTraction.boundaryField(), patchI)
                {
                    const vectorField n = mesh.boundary()[patchI].nf();
                    const symmTensorField& sigmaP =
                        sigma.boundaryField()[patchI];
                    const tensorField& relFP = relF.boundaryField()[patchI];
                    const tensorField relFPinv = inv(relFP);
                    const scalarField relJP = det(relFP);

                    vectorField deformedN = relJP*relFPinv.T() & n;
                    deformedN /= mag(deformedN);

                    totalTraction.boundaryField()[patchI] = deformedN & sigmaP;

                    normalTraction.boundaryField()[patchI] =
                        deformedN & totalTraction.boundaryField()[patchI];

                    shearTraction.boundaryField()[patchI] =
                        (I - sqr(deformedN))
                      & totalTraction.boundaryField()[patchI];
                }

                Info<< "    Writing totalTraction field" << endl;
                totalTraction.write();

                Info<< "    Writing normalTraction field" << endl;
                normalTraction.write();

                Info<< "    Writing shearTraction field" << endl;
                shearTraction.write();
            }
            else
            {
                Info<< "    Detected a total Lagrangian case" << endl;

                const volSymmTensorField sigma(sigmaHeader, mesh);

                const volTensorField F
                (
                    IOobject
                    (
                        "F",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                volVectorField totalTraction
                (
                    IOobject
                    (
                        "totalTraction",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedVector("zero", dimForce/dimArea, vector::zero)
                );

                volScalarField normalTraction
                (
                    IOobject
                    (
                        "normalTraction",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zero", dimForce/dimArea, 0.0)
                );

                volVectorField shearTraction
                (
                    IOobject
                    (
                        "shearTraction",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedVector("zero", dimForce/dimArea, vector::zero)
                );

                forAll(totalTraction.boundaryField(), patchI)
                {
                    const vectorField n = mesh.boundary()[patchI].nf();
                    const symmTensorField& sigmaP =
                        sigma.boundaryField()[patchI];
                    const tensorField& FP = F.boundaryField()[patchI];
                    const tensorField FPinv = inv(FP);
                    const scalarField JP = det(FP);

                    vectorField deformedN = JP*FPinv.T() & n;
                    deformedN /= mag(deformedN);

                    totalTraction.boundaryField()[patchI] = deformedN & sigmaP;

                    normalTraction.boundaryField()[patchI] =
                        deformedN & totalTraction.boundaryField()[patchI];

                    shearTraction.boundaryField()[patchI] =
                        (I - sqr(deformedN))
                      & totalTraction.boundaryField()[patchI];
                }

                Info<< "    Writing totalTraction field" << endl;
                totalTraction.write();

                Info<< "    Writing normalTraction field" << endl;
                normalTraction.write();

                Info<< "    Writing shearTraction field" << endl;
                shearTraction.write();
            }
        }
        else
        {
            Info<< "    No sigma field" << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
