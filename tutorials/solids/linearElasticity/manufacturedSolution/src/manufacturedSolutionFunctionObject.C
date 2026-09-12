/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*----------------------------------------------------------------------------*/

#include "manufacturedSolutionFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manufacturedSolutionFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        manufacturedSolutionFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::manufacturedSolutionFunctionObject::writeData()
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

    // Create the MMS object, if needed
    if (!mmsPtr_.valid())
    {
        mmsPtr_.reset(new manufacturedSolution(mesh, dict_));
    }

    // Lookup the point mesh
    const pointMesh& pMesh = mesh.lookupObject<pointMesh>("pointMesh");

    // Point coordinates
    const pointField& points = mesh.points();

    //Cell-centre coordinates
    const volVectorField& C = mesh.C();
    const vectorField& CI = C;

    // Point analytical fields
    if (pointDisplacement_ || pointStress_)
    {
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

        pointSymmTensorField analyticalPointEpsilon
        (
            IOobject
            (
                "analyticalPointEpsilon",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
            "calculated"
        );

        pointVectorField analyticalPointD
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
            dimensionedVector("zero", dimLength, vector::zero)
        );

        symmTensorField& sI = analyticalStress;
        symmTensorField& pEI = analyticalPointEpsilon;
        vectorField& aPDI = analyticalPointD;
        vectorField& aDI = analyticalD;

        forAll(sI, cellI)
        {
            sI[cellI] = mmsPtr_->calculateStress(CI[cellI]);
            aDI[cellI] = mmsPtr_->calculateDisplacement(CI[cellI]);
        }

        forAll(pEI, pointI)
        {
            pEI[pointI] = mmsPtr_->calculateStrain(points[pointI]);
            aPDI[pointI] = mmsPtr_->calculateDisplacement(points[pointI]);
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                symmTensorField& sP = analyticalStress.boundaryFieldRef()[patchI];
                vectorField& DP = analyticalD.boundaryFieldRef()[patchI];
#else
                symmTensorField& sP = analyticalStress.boundaryField()[patchI];
                vectorField& DP = analyticalD.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                    sP[faceI] = mmsPtr_->calculateStress(CP[faceI]);
                    DP[faceI] = mmsPtr_->calculateDisplacement(CP[faceI]);
                }
            }
        }

        analyticalStress.correctBoundaryConditions();
        analyticalD.correctBoundaryConditions();

        // Write point analytical fields
        if (pointStress_)
        {
            Info<< "Writing analyticalPointStress"
                << nl << endl;
            analyticalStress.write();

            Info<< "Writing analyticalPointEpsilon"
                << nl << endl;
            analyticalPointEpsilon.write();
        }

        if (pointDisplacement_)
        {
            Info<< "Writing analyticalPointDisplacement and analyticalDisplacement"
                << nl << endl;
            analyticalPointD.write();
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

            const pointVectorField diff
            (
                "pointDDifference", analyticalPointD - pointD
            );
            Info<< "Writing pointDDifference field" << endl;
            diff.write();

            const vectorField& diffI = diff;
            Info<< "    Displacement error norms: mean L1, mean L2, LInf: " << nl
                << "    Magnitude: " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << endl;
            for (int cmptI = 0; cmptI < 3; cmptI++)
            {
                Info<< "    " << cmptI << " "
                    << gAverage(mag(diffI.component(cmptI)))
                    << " " << Foam::sqrt(gAverage(magSqr(diffI.component(cmptI))))
                    << " " << gMax(mag(diffI.component(cmptI)))
                    << endl;
            }
        }

        if
        (
            cellDisplacement_
         && mesh.foundObject<volVectorField>("D")
        )
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
            Info<< "    Displacement error norms: mean L1, mean L2, LInf: " << nl
                << "    Magnitude: " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << endl;
            for (int cmptI = 0; cmptI < 3; cmptI++)
            {
                Info<< "    " << cmptI << " "
                    << gAverage(mag(diffI.component(cmptI)))
                    << " " << Foam::sqrt(gAverage(magSqr(diffI.component(cmptI))))
                    << " " << gMax(mag(diffI.component(cmptI)))
                    << endl;
            }
        }

        if
        (
            cellStress_
         && mesh.foundObject<volSymmTensorField>("sigma")
        )
        {
            const volSymmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>("sigma");

            const volSymmTensorField diff
            (
                "sigmaDifference", analyticalStress - sigma
            );
            Info<< "Writing sigmaDifference field" << endl;
            diff.write();

            const symmTensorField& diffI = diff;
            Info<< "    Stress error norms: mean L1, mean L2, LInf: " << nl
                << "    Magnitude: " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << endl;
            for (int cmptI = 0; cmptI < 6; cmptI++)
            {
                Info<< "    " << cmptI << " "
                    << gAverage(mag(diffI.component(cmptI)))
                    << " " << Foam::sqrt(gAverage(magSqr(diffI.component(cmptI))))
                    << " " << gMax(mag(diffI.component(cmptI)))
                    << endl;
            }
        }

        if
        (
            pointStress_
         && mesh.foundObject<pointSymmTensorField>("pEpsilon")
        )
        {
            const pointSymmTensorField& pointEpsilon =
                mesh.lookupObject<pointSymmTensorField>("pEpsilon");

            const pointSymmTensorField diffEpsilon
            (
                "pointEpsilonDifference", analyticalPointEpsilon - pointEpsilon
            );
            Info<< "Writing pointPointEpsilonDifference field" << endl;
            diffEpsilon.write();

            const symmTensorField& diffEpsilonI = diffEpsilon;
            Info<< "    pEpsilon error norms: mean L1, mean L2, LInf: " << nl
                << "    Component XX: " << gAverage(mag(diffEpsilonI.component(0)))
                << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(0))))
                << " " << gMax(mag(diffEpsilonI.component(0)))
                << nl
                << "    Component XY: " << gAverage(mag(diffEpsilonI.component(1)))
                << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(1))))
                << " " << gMax(mag(diffEpsilonI.component(1)))
                << nl
                << "    Component XZ: " << gAverage(mag(diffEpsilonI.component(2)))
                << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(2))))
                << " " << gMax(mag(diffEpsilonI.component(2)))
                << nl
                << "    Component YY: " << gAverage(mag(diffEpsilonI.component(3)))
                << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(3))))
                << " " << gMax(mag(diffEpsilonI.component(3)))
                << nl
                << "    Component YZ: " << gAverage(mag(diffEpsilonI.component(4)))
                << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(4))))
                << " " << gMax(mag(diffEpsilonI.component(4)))
                << nl
                << "    Component ZZ: " << gAverage(mag(diffEpsilonI.component(5)))
                << " " << Foam::sqrt(gAverage(magSqr(diffEpsilonI.component(5))))
                << " " << gMax(mag(diffEpsilonI.component(5)))
                << nl << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manufacturedSolutionFunctionObject::manufacturedSolutionFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    mmsPtr_(),
    dict_(dict),
    cellDisplacement_
    (
        dict.lookupOrDefault<Switch>("cellDisplacement", true)
    ),
    pointDisplacement_
    (
        dict.lookupOrDefault<Switch>("pointDisplacement", true)
    ),
    cellStress_
    (
        dict.lookupOrDefault<Switch>("cellStress", true)
    ),
    pointStress_
    (
        dict.lookupOrDefault<Switch>("pointStress", true)
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::manufacturedSolutionFunctionObject::start()
{
    return true;
}


#if FOAMEXTEND
    bool Foam::manufacturedSolutionFunctionObject::execute(const bool forceWrite)
#else
    bool Foam::manufacturedSolutionFunctionObject::execute()
#endif
{
    return writeData();
}


bool Foam::manufacturedSolutionFunctionObject::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::manufacturedSolutionFunctionObject::write()
{
    return true;
}
#endif

// ************************************************************************* //
