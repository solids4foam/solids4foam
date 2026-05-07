/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    solids4foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with solids4foam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "pressurisedCylinderAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

namespace
{

Foam::tensor cylindricalToCartesianRotation(const Foam::point& pt)
{
    const Foam::scalar r = Foam::mag(Foam::vector(pt.x(), pt.y(), 0));

    if (r < Foam::SMALL)
    {
        return Foam::tensor
        (
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        );
    }

    const Foam::scalar c = pt.x()/r;
    const Foam::scalar s = pt.y()/r;

    return Foam::tensor
    (
        c, -s, 0,
        s,  c, 0,
        0,  0, 1
    );
}

Foam::symmTensor transformStressToCylindrical
(
    const Foam::symmTensor& sigma,
    const Foam::point& pt
)
{
    const Foam::tensor R = cylindricalToCartesianRotation(pt);

    return Foam::symm(R.T() & (sigma & R));
}

Foam::vector transformDisplacementToCylindrical
(
    const Foam::vector& D,
    const Foam::point& pt
)
{
    const Foam::tensor R = cylindricalToCartesianRotation(pt);

    return R.T() & D;
}

}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressurisedCylinderAnalyticalSolution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        pressurisedCylinderAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::pressurisedCylinderAnalyticalSolution::writeData()
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
    if (cellDisplacement_ || cellStress_)
    {
        // Analytical stress field
        volSymmTensorField analyticalStress
        (
            IOobject
            (
                "cylAnalyticalCellStress",
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
                "cylAnalyticalD",
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
            if (cellStress_)
            {
                sI[cellI] = analyticalSol_.cylindricalStress(CI[cellI]);
            }

            if (cellDisplacement_)
            {
                aDI[cellI] = analyticalSol_.cylindricalDisplacement(CI[cellI]);
            }
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
                vectorField& aDP = analyticalD.boundaryFieldRef()[patchI];
#else
                symmTensorField& sP = analyticalStress.boundaryField()[patchI];
                vectorField& aDP = analyticalD.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                    if (cellStress_)
                    {
                        sP[faceI] = analyticalSol_.cylindricalStress(CP[faceI]);
                    }

                    if (cellDisplacement_)
                    {
                        aDP[faceI] =
                            analyticalSol_.cylindricalDisplacement(CP[faceI]);
                    }
                }
            }
        }

        // Write out the cell analytical fields
        if (cellStress_)
        {
            Info<< "Writing cylAnalyticalCellStress field"
                << nl << endl;
            analyticalStress.write();
        }

        if (cellDisplacement_)
        {
            Info<< "Writing cylAnalyticalD field"
                << nl << endl;
            analyticalD.write();
        }

        if (cellStress_ && mesh.foundObject<volSymmTensorField>("sigma"))
        {
            const volSymmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>("sigma");

            volSymmTensorField diff
            (
                IOobject
                (
                    "cylCellStressDifference",
                    time_.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
                "calculated"
            );

#ifdef OPENFOAM_NOT_EXTEND
            symmTensorField& diffI = diff.primitiveFieldRef();
#else
            symmTensorField& diffI = diff.internalField();
#endif

            forAll(diffI, cellI)
            {
                diffI[cellI] =
                    analyticalStress[cellI]
                  - transformStressToCylindrical(sigma[cellI], CI[cellI]);
            }

            forAll(diff.boundaryField(), patchI)
            {
                if (mesh.boundary()[patchI].type() != "empty")
                {
#ifdef OPENFOAM_NOT_EXTEND
                    symmTensorField& diffP = diff.boundaryFieldRef()[patchI];
#else
                    symmTensorField& diffP = diff.boundaryField()[patchI];
#endif
                    const symmTensorField& sigmaP = sigma.boundaryField()[patchI];
                    const vectorField& CP = C.boundaryField()[patchI];

                    forAll(diffP, faceI)
                    {
                        diffP[faceI] =
                            analyticalStress.boundaryField()[patchI][faceI]
                          - transformStressToCylindrical
                            (
                                sigmaP[faceI],
                                CP[faceI]
                            );
                    }
                }
            }

            Info<< "Writing cylCellStressDifference field" << endl;
            diff.write();

            for (int cmpt = 0; cmpt < pTraits<symmTensor>::nComponents; cmpt++)
            {
                // Only calculate for sigmaR, tauRTheta and sigmaTheta
                if (cmpt != 0 && cmpt != 1 && cmpt != 3)
                {
                    continue;
                }

                const scalarField diffCmpt =
                    diff.internalField().component(cmpt);

                Info<< "    Component: " << cmpt << endl;
                Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                    << "    " << gAverage(mag(diffCmpt))
                    << " " << Foam::sqrt(gAverage(magSqr(diffCmpt)))
                    << " " << gMax(mag(diffCmpt))
                    << nl << endl;
            }
        }

        if (cellDisplacement_ && mesh.foundObject<volVectorField>("D"))
        {
            const volVectorField& D =
                mesh.lookupObject<volVectorField>("D");

            volVectorField diff
            (
                IOobject
                (
                    "cylDDifference",
                    time_.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimLength, vector::zero),
                "calculated"
            );

#ifdef OPENFOAM_NOT_EXTEND
            vectorField& diffI = diff.primitiveFieldRef();
#else
            vectorField& diffI = diff.internalField();
#endif

            forAll(diffI, cellI)
            {
                diffI[cellI] =
                    analyticalD[cellI]
                  - transformDisplacementToCylindrical(D[cellI], CI[cellI]);
            }

            forAll(diff.boundaryField(), patchI)
            {
                if (mesh.boundary()[patchI].type() != "empty")
                {
#ifdef OPENFOAM_NOT_EXTEND
                    vectorField& diffP = diff.boundaryFieldRef()[patchI];
#else
                    vectorField& diffP = diff.boundaryField()[patchI];
#endif
                    const vectorField& Dp = D.boundaryField()[patchI];
                    const vectorField& CP = C.boundaryField()[patchI];

                    forAll(diffP, faceI)
                    {
                        diffP[faceI] =
                            analyticalD.boundaryField()[patchI][faceI]
                          - transformDisplacementToCylindrical
                            (
                                Dp[faceI],
                                CP[faceI]
                            );
                    }
                }
            }

            Info<< "Writing cylDDifference field" << endl;
            diff.write();

            const vectorField& diffField = diff.internalField();
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffField))
                << " " << Foam::sqrt(gAverage(magSqr(diffField)))
                << " " << gMax(mag(diffField))
                << nl << endl;
        }
    }

    // Point analytical fields
    if (pointDisplacement_ || pointStress_)
    {
        pointSymmTensorField analyticalStress
        (
            IOobject
            (
                "cylAnalyticalPointStress",
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
                "cylAnalyticalPointD",
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
            if (pointStress_)
            {
                sI[pointI] = analyticalSol_.cylindricalStress(points[pointI]);
            }

            if (pointDisplacement_)
            {
                aDI[pointI] =
                    analyticalSol_.cylindricalDisplacement(points[pointI]);
            }
        }

        // Write point analytical fields
        if (pointStress_)
        {
            Info<< "Writing cylAnalyticalPointStress field"
                << nl << endl;
            analyticalStress.write();
        }

        if (pointDisplacement_)
        {
            Info<< "Writing cylAnalyticalPointDisplacement field"
                << nl << endl;
            analyticalD.write();
        }

        if
        (
            pointDisplacement_
         && mesh.foundObject<pointVectorField>("pointD")
        )
        {
            const pointVectorField& pointD =
                mesh.lookupObject<pointVectorField>("pointD");

            pointVectorField diff
            (
                IOobject
                (
                    "cylPointDDifference",
                    time_.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                pMesh,
                dimensionedVector("zero", dimLength, vector::zero)
            );

            vectorField& diffI = diff;

            forAll(diffI, pointI)
            {
                diffI[pointI] =
                    analyticalD[pointI]
                  - transformDisplacementToCylindrical
                    (
                        pointD[pointI],
                        points[pointI]
                    );
            }

            Info<< "Writing cylPointDDifference field" << endl;
            diff.write();

            const vectorField& diffField = diff.internalField();
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffField))
                << " " << Foam::sqrt(gAverage(magSqr(diffField)))
                << " " << gMax(mag(diffField))
                << nl << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressurisedCylinderAnalyticalSolution::
pressurisedCylinderAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    analyticalSol_(dict),
    cellDisplacement_(dict.lookupOrDefault<Switch>("cellDisplacement", true)),
    pointDisplacement_(dict.lookupOrDefault<Switch>("pointDisplacement", true)),
    cellStress_(dict.lookupOrDefault<Switch>("cellStress", true)),
    pointStress_(dict.lookupOrDefault<Switch>("pointStress", true))
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pressurisedCylinderAnalyticalSolution::start()
{
    return true;
}

#ifdef FOAMEXTEND
bool Foam::pressurisedCylinderAnalyticalSolution::execute(const bool forceWrite)
#else
bool Foam::pressurisedCylinderAnalyticalSolution::execute()
#endif
{
    return writeData();
}

bool Foam::pressurisedCylinderAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}

#ifdef OPENFOAM_NOT_EXTEND
bool Foam::pressurisedCylinderAnalyticalSolution::write()
{
    return false;
}
#endif
// ************************************************************************* //
