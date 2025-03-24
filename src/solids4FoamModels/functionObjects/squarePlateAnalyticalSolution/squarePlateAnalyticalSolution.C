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

#ifndef OPENFOAM_ORG

#include "squarePlateAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "squarePlateStressDisplacement.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(squarePlateAnalyticalSolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        squarePlateAnalyticalSolution,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::squarePlateAnalyticalSolution::writeData()
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

    // Pointer to analytical plate displacement function
    scalar (*plateDisplacement)
    (
        const vector&,
        const scalar,
        const scalar,
        const scalar,
        const scalar,
        const scalar,
        const scalar
    ) = nullptr;

    if (boundaryConditions_ == "allEdgesClamped")
    {
        plateDisplacement = &allEdgesClampedMaxDisplacement;

        // For all edges clamped getting analytical deflection field is quite
        // complicated, so we will print only central deflection
        vector zero = vector::zero;
        const scalar maxD = plateDisplacement(zero, p_, E_, nu_, a_, b_, h_);

        Info<<"Analytical solution for maximal plate deflection at centre "
            << "point is: " << maxD << endl;

        return true;
    }
    else if (boundaryConditions_ == "allEdgesSupported")
    {
        plateDisplacement = &allEdgesSupportedDisplacement;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown boundaryConditions: " << boundaryConditions_ << nl
            << "Valid options are: allEdgesClamped or allEdgesSupported"
            << exit(FatalError);
    }

    // Cell analytical fields
    if (cellDisplacement_)
    {
        // Analytical displacement field
        volScalarField analyticalD
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
            dimensionedScalar("zero", dimLength, 0.0),
            "calculated"
        );

        scalarField& aDI = analyticalD;

        forAll(aDI, cellI)
        {
            aDI[cellI] = plateDisplacement(CI[cellI], p_, E_, nu_, a_, b_, h_);
        }

        forAll(analyticalD.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
#ifdef OPENFOAM_NOT_EXTEND
                scalarField& aDP = analyticalD.boundaryFieldRef()[patchI];
#else
                scalarField& aDP = analyticalD.boundaryField()[patchI];
#endif
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(aDP, faceI)
                {
                    aDP[faceI] =
                        plateDisplacement
                        (
                            CP[faceI], p_, E_, nu_, a_, b_, h_
                        );
                }
            }
        }

        Info<< "Writing analyticalD field" << nl << endl;
        analyticalD.write();

        if (cellDisplacement_ && mesh.foundObject<volVectorField>("D"))
        {
            const volScalarField& D =
                mesh.lookupObject<volScalarField>("wVf");

            const volScalarField diff
            (
                "DDifference", analyticalD - D
            );
            Info<< "Writing DDifference field" << endl;
            diff.write();

            const scalarField& diffI = diff.internalField();
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << nl << endl;
        }
    }

    // Point analytical fields
    if (pointDisplacement_)
    {
        pointScalarField analyticalD
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
            dimensionedScalar("zero", dimLength, 0.0)
        );

        scalarField& aDI = analyticalD;

        forAll(aDI, pointI)
        {
            aDI[pointI] =
                plateDisplacement
                (
                    points[pointI], p_, E_, nu_, a_, b_, h_
                );
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

            const pointScalarField diff
            (
                "pointDDifference", analyticalD - (pointD & vector(0,0,1))
            );
            Info<< "Writing pointDDifference field" << endl;
            diff.write();

            const scalarField& diffI = diff.internalField();
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

Foam::squarePlateAnalyticalSolution::squarePlateAnalyticalSolution
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    p_(readScalar(dict.lookup("p"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    a_(readScalar(dict.lookup("a"))),
    b_(readScalar(dict.lookup("b"))),
    h_(readScalar(dict.lookup("h"))),
    boundaryConditions_
    (
        dict.get<word>("boundaryConditions")
    ),
    cellDisplacement_
    (
        dict.lookupOrDefault<Switch>("cellDisplacement", true)
    ),
    pointDisplacement_
    (
        dict.lookupOrDefault<Switch>("pointDisplacement", false)
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (a_ < SMALL || b_ < SMALL || h_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "a, b and h should be greater than 0!"
            << abort(FatalError);
    }

    if (E_ < SMALL || nu_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "E and nu should be positive!"
            << abort(FatalError);
    }

     if ((sqrt(a_*b_)/h_) < 10)
    {
        InfoInFunction << this->name() << " function object constructor"
            << "Analytical solution is derived using Kirchhoff–Love shell"
            << "theory which is applicable only in cases L/h > 10"
            << abort(FatalError);
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::squarePlateAnalyticalSolution::start()
{
    return true;
}


#ifdef FOAMEXTEND
bool Foam::squarePlateAnalyticalSolution::execute(const bool forceWrite)
#else
bool Foam::squarePlateAnalyticalSolution::execute()
#endif
{
    return writeData();
}


bool Foam::squarePlateAnalyticalSolution::read(const dictionary& dict)
{
    return true;
}


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::squarePlateAnalyticalSolution::write()
{
    return false;
}
#endif

#endif // #ifndef OPENFOAM_ORG

// ************************************************************************* //
