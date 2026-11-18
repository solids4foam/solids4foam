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
#include "lookupSolidModel.H"
#include "kExactLeastSquares.H"
#include "compatibilityFunctions.H"

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

void Foam::manufacturedSolutionFunctionObject::
calculateAnalyticalCellDisplacement
(
    const fvMesh& mesh,
    vectorField& analyticalD
)
{
    const solidModel& solMod = lookupSolidModel(mesh);
    bool useCellAverage = false;

    if (solMod.highOrderResidual())
    {
        useCellAverage =
            isA<kExactLeastSquares>(solMod.displacementLeastSquares());
    }

    if (useCellAverage)
    {
        const fvMeshQuadrature& quadrature =
            solMod.displacementLeastSquares().quadrature();
        auto& cellQuadPoints =
            compactListListCRef(quadrature.cellQuadPoints());
        auto& cellQuadWeights =
            compactListListCRef(quadrature.cellQuadWeights());
        const scalarField& cellVolumes = mesh.V();

        Info<< "Using cell-average analytical displacement" << endl;

        forAll(analyticalD, cellI)
        {
            analyticalD[cellI] = vector::zero;

            forAll(cellQuadPoints[cellI], pointI)
            {
                analyticalD[cellI] +=
                    cellQuadWeights[cellI][pointI]
                   *mmsPtr_->calculateDisplacement
                    (
                        cellQuadPoints[cellI][pointI]
                    );
            }

            analyticalD[cellI] /= cellVolumes[cellI];
        }
    }
    else
    {
        Info<< "Using point-valued analytical displacement" << endl;

        const vectorField& cellCentres = mesh.C().internalField();
        forAll(analyticalD, cellI)
        {
            analyticalD[cellI] =
                mmsPtr_->calculateDisplacement(cellCentres[cellI]);
        }
    }
}


void Foam::manufacturedSolutionFunctionObject::writePressureError
(
    const fvMesh& mesh
)
{
    // Read solvePressure as the solid model itself does (solidModel.C);
    // solidModel::solvePressure() is not public
    const bool solvePressure =
        lookupSolidModel(mesh).solidModelDict().lookupOrDefault<Switch>
        (
            "solvePressure", false
        );

    const bool reportPressure =
        cellPressureSpecified_ ? bool(cellPressure_) : solvePressure;

    if (!reportPressure)
    {
        return;
    }

    if (!solvePressure || !mesh.foundObject<volScalarField>("p"))
    {
        FatalErrorInFunction
            << "cellPressure is enabled but the solid model does not solve "
            << "for a pressure field p: enable solvePressure or set "
            << "cellPressure no" << exit(FatalError);
    }

    // The mixed formulation uses the compression-positive pressure,
    //     sigma = sigmaTilde - p*I,  with  -p/K - tr(grad(D)) = 0,
    // so for this linear elastic solution the analytical pressure is
    //     p = -tr(sigma)/3 = -(lambda + 2*mu/3)*div(u) = -K*div(u),
    // evaluated from the analytical stress at the cell centres
    volScalarField analyticalP
    (
        IOobject
        (
            "analyticalP",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0),
        "calculated"
    );

    const vectorField& CI = mesh.C();
    scalarField& aPI = analyticalP;
    forAll(aPI, cellI)
    {
        aPI[cellI] = -tr(mmsPtr_->calculateStress(CI[cellI]))/3.0;
    }

    forAll(analyticalP.boundaryField(), patchI)
    {
        if (mesh.boundary()[patchI].type() != "empty")
        {
#ifdef OPENFOAM_NOT_EXTEND
            scalarField& pP = analyticalP.boundaryFieldRef()[patchI];
#else
            scalarField& pP = analyticalP.boundaryField()[patchI];
#endif
            const vectorField& CP = mesh.C().boundaryField()[patchI];

            forAll(pP, faceI)
            {
                pP[faceI] = -tr(mmsPtr_->calculateStress(CP[faceI]))/3.0;
            }
        }
    }

    Info<< "Writing analyticalP" << nl << endl;
    analyticalP.write();

    const volScalarField& p = mesh.lookupObject<volScalarField>("p");

    const volScalarField diff("pDifference", analyticalP - p);
    Info<< "Writing pDifference field" << endl;
    diff.write();

    // Cell-wise norms, as for the displacement and stress above: the mean of
    // |e| and the root mean square over cells (unweighted), and the maximum.
    // The volume-weighted RMS, sqrt(sum(V*e^2)/sum(V)), is the continuous L2
    // norm of the cell-wise constant error divided by sqrt(domain volume); it
    // equals the cell RMS on a uniform mesh only
    const scalarField& diffI = diff;
    const scalarField& V = mesh.V();
    Info<< "    Pressure error norms: mean L1, cell RMS, LInf, "
        << "volume-weighted RMS: " << nl
        << "    Magnitude: " << gAverage(mag(diffI))
        << " " << Foam::sqrt(gAverage(sqr(diffI)))
        << " " << gMax(mag(diffI))
        << " " << Foam::sqrt(gSum(V*sqr(diffI))/gSum(V))
        << endl;
}


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
        }

        calculateAnalyticalCellDisplacement(mesh, aDI);

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


    // Pressure error of the mixed formulation (optional)
    writePressureError(mesh);

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
    ),
    cellPressureSpecified_(dict.found("cellPressure")),
    cellPressure_
    (
        dict.lookupOrDefault<Switch>("cellPressure", false)
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
