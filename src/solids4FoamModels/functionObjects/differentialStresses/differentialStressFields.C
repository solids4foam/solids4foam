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

\*----------------------------------------------------------------------------*/

#include "differentialStressFields.H"
#include "principalStressFields.H"

// * * * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * //

void Foam::calculateDifferentialStress
(
    const symmTensor& sigma,
    const bool& compressionPositive,
    scalar& sigmaDiff
)
{
    scalar sigmaMax = 0;
    scalar sigmaMid = 0;
    scalar sigmaMin = 0;
    vector sigmaMaxDir = vector::zero;
    vector sigmaMidDir = vector::zero;
    vector sigmaMinDir = vector::zero;

    calculateEigenValues
    (
        sigma,
        compressionPositive,
        sigmaMax,
        sigmaMid,
        sigmaMin,
        sigmaMaxDir,
        sigmaMidDir,
        sigmaMinDir
    );

    sigmaDiff = sigmaMax - sigmaMin;
}


void Foam::writedifferentialStressFields(const volSymmTensorField& sigma, const bool& compressionPositive)
{
    const fvMesh& mesh = sigma.mesh();
    const Time& runTime = mesh.time();

    // Differential stress
    volScalarField sigmaDiff
    (
        IOobject
        (
            "sigmaDiff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sigmaDiff", dimPressure, 0.0)
    );

    // References to internalFields for efficiency
    const symmTensorField& sigmaI = sigma.internalField();
#ifdef OPENFOAMESIORFOUNDATION
    scalarField& sigmaDiffI = sigmaDiff.primitiveFieldRef();
#else
    scalarField& sigmaDiffI = sigmaDiff.internalField();
#endif

    forAll (sigmaI, cellI)
    {
        calculateDifferentialStress
        (
            sigmaI[cellI],
            compressionPositive,
            sigmaDiffI[cellI]
        );
    }

    forAll(sigma.boundaryField(), patchI)
    {
        if
        (
            !sigma.boundaryField()[patchI].coupled()
          && mesh.boundaryMesh()[patchI].type() != "empty"
        )
        {
            const symmTensorField& pSigma = sigma.boundaryField()[patchI];
#ifdef OPENFOAMESIORFOUNDATION
            scalarField& pSigmaDiff = sigmaDiff.boundaryFieldRef()[patchI];
#else
            scalarField& pSigmaDiff = sigmaDiff.boundaryField()[patchI];
#endif

            forAll(pSigmaDiff, faceI)
            {
                calculateDifferentialStress
                (
                    pSigma[faceI],
                    compressionPositive,
                    pSigmaDiff[faceI]
                );
            }
        }
    }

    sigmaDiff.correctBoundaryConditions();

    // Write fields
    // Info<< "    Writing sigmaDiff" << endl;

    sigmaDiff.write();

    Info<< "differential stresses: max = " << gMax(mag(sigmaDiff)()) << endl;
}


// ************************************************************************* //
