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

\*---------------------------------------------------------------------------*/

#include "neoHookeanManufacturedSolutionFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "lookupSolidModel.H"
#include "kExactLeastSquares.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanManufacturedSolutionFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        neoHookeanManufacturedSolutionFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::neoHookeanManufacturedSolutionFunctionObject::writeNorms
(
    const word& title,
    const Field<Type>& diff
) const
{
    Info<< "    " << title << " error norms: mean L1, mean L2, LInf: " << nl
        << "    Magnitude: " << gAverage(mag(diff))
        << " " << Foam::sqrt(gAverage(magSqr(diff)))
        << " " << gMax(mag(diff))
        << endl;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Info<< "    " << cmpt << " "
            << gAverage(mag(diff.component(cmpt)))
            << " " << Foam::sqrt(gAverage(magSqr(diff.component(cmpt))))
            << " " << gMax(mag(diff.component(cmpt)))
            << endl;
    }
}


bool Foam::neoHookeanManufacturedSolutionFunctionObject::writeData()
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
        mmsPtr_.reset(new neoHookeanManufacturedSolution(mesh));
    }
    const neoHookeanManufacturedSolution& mms = mmsPtr_();

    const scalar t = time_.value();
    const pointMesh& pMesh = pointMesh::New(mesh);
    const pointField& points = mesh.points();
    const volVectorField& C = mesh.C();

#ifdef FOAMEXTEND
    const bool writeTime = time_.outputTime();
#else
    const bool writeTime = time_.writeTime();
#endif

    // Analytical fields

    volVectorField analyticalD
    (
        IOobject
        (
            "analyticalD",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero),
        "calculated"
    );

    volSymmTensorField analyticalSigma
    (
        IOobject
        (
            "analyticalSigma",
            time_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
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
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    {
        vectorField& aDI = analyticalD;
        symmTensorField& aSI = analyticalSigma;
        const vectorField& CI = C;

        forAll(CI, cellI)
        {
            aDI[cellI] = mms.displacement(CI[cellI], t);
            aSI[cellI] = mms.cauchyStress(CI[cellI], t);
        }

        // kExactLeastSquares stores cell averages, so its displacement is
        // compared with the cell average of the exact solution
        const solidModel& solMod = lookupSolidModel(mesh);
        if
        (
            solMod.highOrderResidual()
         && isA<kExactLeastSquares>(solMod.displacementLeastSquares())
        )
        {
            const fvMeshQuadrature& quadrature =
                solMod.displacementLeastSquares().quadrature();
            auto& cellQuadPoints =
                compactListListCRef(quadrature.cellQuadPoints());
            auto& cellQuadWeights =
                compactListListCRef(quadrature.cellQuadWeights());

            Info<< "Using cell-average analytical displacement" << endl;

            forAll(aDI, cellI)
            {
                aDI[cellI] = vector::zero;

                forAll(cellQuadPoints[cellI], pointI)
                {
                    aDI[cellI] +=
                        cellQuadWeights[cellI][pointI]
                       *mms.displacement(cellQuadPoints[cellI][pointI], t);
                }

                aDI[cellI] /= mesh.V()[cellI];
            }
        }

        forAll(mesh.boundary(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                vectorField& aDP = analyticalD.boundaryFieldRef()[patchI];
                symmTensorField& aSP =
                    analyticalSigma.boundaryFieldRef()[patchI];
#else
                vectorField& aDP = analyticalD.boundaryField()[patchI];
                symmTensorField& aSP =
                    analyticalSigma.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(CP, faceI)
                {
                    aDP[faceI] = mms.displacement(CP[faceI], t);
                    aSP[faceI] = mms.cauchyStress(CP[faceI], t);
                }
            }
        }

        vectorField& aPDI = analyticalPointD;
        forAll(points, pointI)
        {
            aPDI[pointI] = mms.displacement(points[pointI], t);
        }
    }

    // Errors

    if (mesh.foundObject<volVectorField>("D"))
    {
        const volVectorField& D = mesh.lookupObject<volVectorField>("D");

        const volVectorField diff("DDifference", analyticalD - D);
        Info<< "Writing DDifference field" << endl;
        writeNorms("Displacement", vectorField(diff));

        if (writeTime)
        {
            analyticalD.write();
            diff.write();
        }
    }

    if (mesh.foundObject<pointVectorField>("pointD"))
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        const pointVectorField diff
        (
            "pointDDifference", analyticalPointD - pointD
        );
        Info<< "Writing pointDDifference field" << endl;
        writeNorms("Point displacement", vectorField(diff));

        if (writeTime)
        {
            analyticalPointD.write();
            diff.write();
        }
    }

    if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        const volSymmTensorField diff
        (
            "sigmaDifference", analyticalSigma - sigma
        );
        Info<< "Writing sigmaDifference field" << endl;
        writeNorms("Stress", symmTensorField(diff));

        if (writeTime)
        {
            analyticalSigma.write();
            diff.write();
        }
    }

    Info<< endl;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neoHookeanManufacturedSolutionFunctionObject::
neoHookeanManufacturedSolutionFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    mmsPtr_()
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::neoHookeanManufacturedSolutionFunctionObject::start()
{
    return true;
}


#if FOAMEXTEND
bool Foam::neoHookeanManufacturedSolutionFunctionObject::execute
(
    const bool forceWrite
)
#else
bool Foam::neoHookeanManufacturedSolutionFunctionObject::execute()
#endif
{
    return writeData();
}


bool Foam::neoHookeanManufacturedSolutionFunctionObject::read
(
    const dictionary& dict
)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::neoHookeanManufacturedSolutionFunctionObject::write()
{
    return true;
}
#endif

// ************************************************************************* //
