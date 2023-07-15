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
    if (cellDisplacement_ || cellStress_)
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
            if (cellStress_)
            {
                sI[cellI] =
                    cantileverStress(CI[cellI], P_, E_, nu_, L_, D_, I_);
            }

            if (cellDisplacement_)
            {
                aDI[cellI] =
                    cantileverDisplacement(CI[cellI], P_, E_, nu_, L_, D_, I_);
            }
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
                    if (cellStress_)
                    {
                        sP[faceI] =
                            cantileverStress
                            (
                                CP[faceI], P_, E_, nu_, L_, D_, I_
                            );
                    }

                    if (cellDisplacement_)
                    {
                        aDP[faceI] =
                            cantileverDisplacement
                            (
                                CP[faceI], P_, E_, nu_, L_, D_, I_
                            );
                    }
                }
            }
        }

        // Write out the cell analytical fields
        if (cellStress_)
        {
            Info<< "Writing analyticalCellStress field"
                << nl << endl;
            analyticalStress.write();
        }

        if (cellDisplacement_)
        {
            Info<< "Writing analyticalD field"
                << nl << endl;
            analyticalD.write();
        }


        if (cellStress_ && mesh.foundObject<volSymmTensorField>("sigma"))
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
                // Only calculate for XX, XY and YY
                if (cmpt != 0 && cmpt != 1 && cmpt != 3)
                {
                    continue;
                }

                const scalarField diffI = diff.internalField().component(cmpt);

                Info<< "    Component: " << cmpt << endl;
                Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                    << "    " << gAverage(mag(diffI))
                    << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                    << " " << gMax(mag(diffI))
                    << nl << endl;
            }
        }

        if (cellDisplacement_ && mesh.foundObject<volVectorField>("D"))
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
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
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
            if (pointStress_)
            {
                sI[pointI] =
                    cantileverStress(points[pointI], P_, E_, nu_, L_, D_, I_);
            }

            if (pointDisplacement_)
            {
                aDI[pointI] =
                    cantileverDisplacement
                    (
                        points[pointI], P_, E_, nu_, L_, D_, I_
                    );
            }
        }

        // Write point analytical fields
        if (pointStress_)
        {
            Info<< "Writing analyticalPointStress field"
                << nl << endl;
            analyticalStress.write();
        }

        if (pointDisplacement_)
        {
            Info<< "Writing analyticalPointDisplacement field"
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

            const pointVectorField diff
            (
                "pointDDifference", analyticalD - pointD
            );
            Info<< "Writing pointDDifference field" << endl;
            diff.write();

            const vectorField& diffI = diff.internalField();
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << nl << endl;
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
    I_(Foam::pow(D_, 3.0)/12.0),
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
    return false;
}
#endif

// ************************************************************************* //
