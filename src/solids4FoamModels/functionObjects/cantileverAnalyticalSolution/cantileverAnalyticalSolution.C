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

#include "cantileverAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "cantileverStressDisplacement.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cantileverAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cantileverAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::cantileverAnalyticalSolution::writeData()
{
    // Lookup the solid mesh
    const fvMesh* meshPtr = NULL;
    if (time_.foundObject<fvMesh>("solid"))
    {
        meshPtr = &(time_.lookupObject<fvMesh>("solid"));
    }
    else
    {
        meshPtr = &(time_.lookupObject<fvMesh>("region0"));
    }
    const fvMesh& mesh = *meshPtr;

    // Lookup the point mesh
    const pointMesh& pMesh = mesh.lookupObject<pointMesh>("pointMesh");

    // Cell-centres coordinates
    const volVectorField& C = mesh.C();
    const vectorField& CI = C.internalField();

    // Point coordinates
    const pointField& points = mesh.points();

    // Cell analytical fields
    {
        // Analytical stress field
        volSymmTensorField analyticalStress
        (
            IOobject
            (
                "analyticalCellStress",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
            "calculated"
        );

        // Analytical displacement field
        volVectorField analyticalD
        (
            IOobject
            (
                "analyticalD",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimLength, vector::zero),
            "calculated"
        );

        symmTensorField& sI = analyticalStress.internalField();
        vectorField& aDI = analyticalD.internalField();

        forAll(sI, cellI)
        {
            sI[cellI] = cantileverStress(CI[cellI], P_, E_, nu_, L_, D_, I_);
            aDI[cellI] =
                cantileverDisplacement(CI[cellI], P_, E_, nu_, L_, D_, I_);
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
                symmTensorField& sP = analyticalStress.boundaryField()[patchI];
                vectorField& aDP = analyticalD.boundaryField()[patchI];
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                    sP[faceI] =
                        cantileverStress(CP[faceI], P_, E_, nu_, L_, D_, I_);
                    aDP[faceI] =
                        cantileverDisplacement(CP[faceI], P_, E_, nu_, L_, D_, I_);
                }
            }
        }

        // Write out the cell analytical field
        Info<< "Writing analytical stress field (analyticalCellStress)"
            << endl;
        analyticalStress.write();
        Info<< "Writing analytical stress field (analyticalCellStress)"
            << endl;
        analyticalD.write();


        if (mesh.foundObject<volSymmTensorField>("sigma"))
        {
            const volSymmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>("sigma");

            const volSymmTensorField diff
            (
                "cellStressDifference", analyticalStress - sigma
            );
            Info<< "Writing cellStressDifference field" << endl;
            diff.write();

            for (int cmpt = 0; cmpt < pTraits<symmTensor>::nComponents; cmpt++)
            {
                // Only calculate for XX and XY
                if (cmpt != 0 && cmpt != 1)
                {
                    continue;
                }

                const scalarField diffI = diff.internalField().component(cmpt);

                Info<< "    Component: " << cmpt << endl;
                Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                    << "    " << gAverage(mag(diffI))
                    << " " << Foam::sqrt(gSum(magSqr(diffI)))/diffI.size()
                    << " " << gMax(mag(diffI))
                    << endl;
            }
        }

        if (mesh.foundObject<volVectorField>("D"))
        {
            const volVectorField& D =
                mesh.lookupObject<volVectorField>("D");

            const volVectorField diff
            (
                "DDifference", analyticalD - D
            );
            Info<< "Writing DDifference field" << endl;
            diff.write();

            const vectorField& diffI = diff.internalField();
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gSum(magSqr(diffI)))/diffI.size()
                << " " << gMax(mag(diffI))
                << endl;
        }
    }

    // Point analytical fields
    {
        pointSymmTensorField analyticalStress
        (
            IOobject
            (
                "analyticalPointStress",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        );

        pointVectorField analyticalD
        (
            IOobject
            (
                "analyticalPointD",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimLength, vector::zero)
        );

        symmTensorField& sI = analyticalStress.internalField();
        vectorField& aDI = analyticalD.internalField();

        forAll(sI, pointI)
        {
            sI[pointI] =
                cantileverStress(points[pointI], P_, E_, nu_, L_, D_, I_);
            aDI[pointI] =
                cantileverDisplacement(points[pointI], P_, E_, nu_, L_, D_, I_);
        }

        // Write point analytical fields
        Info<< "Writing analyticalPointStress"
            << endl;
        analyticalStress.write();
        Info<< "Writing analyticalPointDisplacement"
            << endl;
        analyticalD.write();

        if (mesh.foundObject<pointVectorField>("pointD"))
        {
            const pointVectorField& pointD =
                mesh.lookupObject<pointVectorField>("pointD");

            const pointVectorField diff
            (
                "pointDDifference", analyticalD - pointD
            );
            Info<< "Writing pointDDifference field" << endl;
            diff.write();

            const vectorField& diffI = diff.internalField();
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gSum(magSqr(diffI)))/diffI.size()
                << " " << gMax(mag(diffI))
                << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cantileverAnalyticalSolution::cantileverAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    P_(readScalar(dict.lookup("P"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    L_(readScalar(dict.lookup("L"))),
    D_(readScalar(dict.lookup("D"))),
    I_(Foam::pow(D_, 3.0)/12.0)
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (L_ < SMALL || D_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "L and D should both be greater than 0!"
            << abort(FatalError);
    }

    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "E and nu should be positive!"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cantileverAnalyticalSolution::start()
{
    return true;
}


#ifdef FOAMEXTEND
bool Foam::cantileverAnalyticalSolution::execute(const bool forceWrite)
#else
bool Foam::cantileverAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::cantileverAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::cantileverAnalyticalSolution::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
