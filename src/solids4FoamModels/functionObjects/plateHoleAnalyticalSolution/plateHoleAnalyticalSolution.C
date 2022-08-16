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

#include "plateHoleAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plateHoleAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        plateHoleAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor Foam::plateHoleAnalyticalSolution::plateHoleStress
(
    const vector& C
)
{
    tensor sigma = tensor::zero;

    // Calculate radial coordinate
    const scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));

    // Calculate circumferential coordinate
    const scalar theta = Foam::atan2(C.y(), C.x());

    const coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    sigma.xx() =
        T_*(1 - sqr(holeR_)/sqr(r))/2
      + T_
       *(1 + 3*pow(holeR_,4)/pow(r,4) - 4*sqr(holeR_)/sqr(r))*::cos(2*theta)/2;

    sigma.xy() =
      - T_
       *(1 - 3*pow(holeR_,4)/pow(r,4) + 2*sqr(holeR_)/sqr(r))*::sin(2*theta)/2;

    sigma.yx() = sigma.xy();

    sigma.yy() =
        T_*(1 + sqr(holeR_)/sqr(r))/2
      - T_*(1 + 3*pow(holeR_,4)/pow(r,4))*::cos(2*theta)/2;


    // Transformation to Cartesian coordinate system
#ifdef OPENFOAMFOUNDATION
    sigma = ((cs.R().R() & sigma) & cs.R().R().T());
#else
    sigma = ((cs.R() & sigma) & cs.R().T());
#endif

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}


Foam::vector Foam::plateHoleAnalyticalSolution::plateHoleDisplacement
(
    const vector& C, const symmTensor& sigma
)
{
    // Shear modulus
    const scalar mu = E_/(2*(1 + nu_));

    // Kappa parameter
    const scalar kappa = 3 - 4*nu_;

    // Polar coordinates
    const scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    const scalar theta = atan2(C.y(), C.x());

    return vector
    (
        (holeR_*T_/(8*mu))
        *(
            (r/holeR_)*(kappa + 1)*cos(theta)
          + (2*holeR_/r)*((1 + kappa)*cos(theta) + cos(3*theta))
          - (2*pow(holeR_, 3)/pow(r,3))*cos(3*theta)
        ),
        (holeR_*T_/(8*mu))
        *(
            (r/holeR_)*(kappa - 3)*sin(theta)
          + (2*holeR_/r)*((1 - kappa)*sin(theta) + sin(3*theta))
          - (2*pow(holeR_, 3)/pow(r,3))*sin(3*theta)
        ),
        0.0
    );
}

bool Foam::plateHoleAnalyticalSolution::writeData()
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
    const vectorField& CI = C;

    // Point coordinates
    const pointField& points = mesh.points();

    if (gMin(mag(points)) < SMALL)
    {
        FatalErrorIn("bool Foam::plateHoleAnalyticalSolution::writeData()")
            << "The hole should be centred on the origin!"
            << abort(FatalError);
    }

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

        symmTensorField& sI = analyticalStress;
        vectorField& aDI = analyticalD;

        forAll(sI, cellI)
        {
            sI[cellI] = plateHoleStress(CI[cellI]);
            aDI[cellI] = plateHoleDisplacement(CI[cellI], sI[cellI]);
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAMESIORFOUNDATION
                symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
                vectorField& aDP = analyticalD.boundaryFieldRef()[patchI];
#else
                symmTensorField& sP = analyticalStress.boundaryField()[patchI];
                vectorField& aDP = analyticalD.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                    sP[faceI] = plateHoleStress(CP[faceI]);
                    aDP[faceI] = plateHoleDisplacement(CP[faceI], sP[faceI]);
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
                // Only calculate for XX, XY and ZZ
                if (cmpt != 0 && cmpt != 1 && cmpt != 3)
                {
                    continue;
                }

                const symmTensorField& diffI = diff;
                const scalarField diffIcmptI(diffI.component(cmpt));

                Info<< "    Component: " << cmpt << endl;
                Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                    << "    " << gAverage(mag(diffIcmptI))
                    << " " << Foam::sqrt(gAverage(magSqr(diffIcmptI)))
                    << " " << gMax(mag(diffIcmptI))
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

            const vectorField& diffI = diff;
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
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

        symmTensorField& sI = analyticalStress;
        vectorField& aDI = analyticalD;

        forAll(sI, pointI)
        {
            sI[pointI] = plateHoleStress(points[pointI]);
            aDI[pointI] = plateHoleDisplacement(points[pointI], sI[pointI]);
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

            const vectorField& diffI = diff;
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plateHoleAnalyticalSolution::plateHoleAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu")))
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (holeR_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "holeRadius should be greater than 0!"
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

bool Foam::plateHoleAnalyticalSolution::start()
{
    return true;
}


#if FOAMEXTEND
    bool Foam::plateHoleAnalyticalSolution::execute(const bool forceWrite)
#else
    bool Foam::plateHoleAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::plateHoleAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAMESIORFOUNDATION
bool Foam::plateHoleAnalyticalSolution::write()
{
    return writeData();
}
#endif

// ************************************************************************* //
