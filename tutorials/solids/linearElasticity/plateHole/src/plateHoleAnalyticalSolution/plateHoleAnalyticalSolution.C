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

#include "plateHoleAnalyticalSolution.H"
#include "addToRunTimeSelectionTable.H"
#include "lookupSolidModel.H"
#include "pointFields.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "plateHoleAnalyticalFields.H"
#include "compatibilityFunctions.H"

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

const Foam::fvMesh& Foam::plateHoleAnalyticalSolution::mesh() const
{
    if (time_.foundObject<fvMesh>(regionName_))
    {
        return time_.lookupObject<fvMesh>(regionName_);
    }
    else if (time_.foundObject<fvMesh>("solid"))
    {
        return time_.lookupObject<fvMesh>("solid");
    }

    return time_.lookupObject<fvMesh>("region0");
}


bool Foam::plateHoleAnalyticalSolution::planeStress(const fvMesh& mesh) const
{
    if (mesh.foundObject<IOdictionary>("mechanicalProperties"))
    {
        return Switch
        (
            mesh.lookupObject<IOdictionary>
            (
                "mechanicalProperties"
            ).lookup("planeStress")
        );
    }
    else if
    (
        mesh.objectRegistry::parent().foundObject<objectRegistry>("region0")
    )
    {
        return Switch
        (
            mesh.objectRegistry::parent().subRegistry
            (
                "region0"
            ).lookupObject<IOdictionary>
            (
                "mechanicalProperties"
            ).lookup("planeStress")
        );
    }

    return Switch
    (
        mesh.objectRegistry::parent().subRegistry
        (
            "solid"
        ).lookupObject<IOdictionary>
        (
            "mechanicalProperties"
        ).lookup("planeStress")
    );
}


Foam::symmTensor Foam::plateHoleAnalyticalSolution::plateHoleStress
(
    const vector& C
) const
{
    return plateHoleAnalyticalFields::stress(C, T_, holeR_);
}


Foam::vector Foam::plateHoleAnalyticalSolution::plateHoleDisplacement
(
    const vector& C,
    const fvMesh& mesh
) const
{
    if (!deriveMaterialProperties_)
    {
        return plateHoleAnalyticalFields::displacement
        (
            C,
            T_,
            holeR_,
            E_,
            nu_,
            planeStress(mesh)
        );
    }

    const solidModel& solMod = lookupSolidModel(mesh);

    const scalar mu = solMod.mechanical().shearModulus()()[0];
    const scalar K = solMod.mechanical().bulkModulus()()[0];

    scalar nu = 0.5;
    if (K + SMALL < GREAT)
    {
        nu = (3*K - 2*mu)/(2*(3*K + mu));

        if (planeStress(mesh))
        {
            nu = (K - mu)/(K + mu);
        }
    }

    scalar kappa = 3 - 4*nu;
    if (planeStress(mesh))
    {
        kappa = (3.0 - nu)/(1.0 + nu);
    }

    return plateHoleAnalyticalFields::displacement
    (
        C,
        T_,
        holeR_,
        mu,
        kappa
    );
}


Foam::scalar Foam::plateHoleAnalyticalSolution::plateHoleHydPressure
(
    const vector& C,
    const fvMesh& mesh
) const
{
    scalar nu = nu_;

    if (deriveMaterialProperties_)
    {
        const solidModel& solMod = lookupSolidModel(mesh);

        const scalar mu = solMod.mechanical().shearModulus()()[0];
        const scalar K = solMod.mechanical().bulkModulus()()[0];

        nu = 0.5;
        if (K + SMALL < GREAT)
        {
            nu = (3*K - 2*mu)/(2*(3*K + mu));

            if (planeStress(mesh))
            {
                nu = (K - mu)/(K + mu);
            }
        }
    }

    return plateHoleAnalyticalFields::hydPressure(C, T_, holeR_, nu);
}


void Foam::plateHoleAnalyticalSolution::writePressureDisplacementData
(
    const fvMesh& mesh,
    const pointMesh& pMesh
) const
{
    if (time_.outputTime())
    {
        volVectorField Derror
        (
            IOobject
            (
                "Derror",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimLength, vector::zero)
        );

        const volVectorField& D =
            mesh.lookupObject<volVectorField>("D");

        const vectorField& DI = D.internalField();
        const vectorField& C = mesh.C().internalField();
        vectorField& DError = Foam::primitiveFieldRef(Derror);

        forAll(DError, cellI)
        {
            const vector curR(C[cellI].x(), C[cellI].y(), 0);
            const vector curDa = plateHoleDisplacement(curR, mesh);

            DError[cellI] = DI[cellI] - curDa;
        }

        Info<< "DError, max : " << gMax(mag(DError)) << endl;
        Derror.write();
    }

    if (mesh.foundObject<pointVectorField>("pointD"))
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        const vectorField& pointDI = pointD.internalField();
        const vectorField& points = mesh.points();
        scalarField pointDError(pointDI.size(), 0);

        forAll(pointDError, pointI)
        {
            const vector curR(points[pointI].x(), points[pointI].y(), 0);
            const vector curDa = plateHoleDisplacement(curR, mesh);

            pointDError[pointI] = mag(pointDI[pointI] - curDa);
        }

        Info<< "pointDError, max : " << gMax(pointDError) << endl;

        if (time_.outputTime())
        {
            pointVectorField pointDerror
            (
                IOobject
                (
                    "pointDerror",
                    time_.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                pMesh,
                dimensionedVector("zero", dimLength, vector::zero)
            );

            forAll(pointDError, pointI)
            {
                const vector curR(points[pointI].x(), points[pointI].y(), 0);
                const vector curDa = plateHoleDisplacement(curR, mesh);

                Foam::primitiveFieldRef(pointDerror)[pointI] =
                    pointDI[pointI] - curDa;
            }

            pointDerror.write();
        }
    }

    if
    (
        time_.outputTime()
     && mesh.foundObject<volSymmTensorField>("sigma")
     && mesh.foundObject<volScalarField>("p")
    )
    {
        volScalarField sigmaXXErr
        (
            IOobject
            (
                "sigmaXXErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0)
        );

        volScalarField sigmaXYErr
        (
            IOobject
            (
                "sigmaXYErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0)
        );

        volScalarField sigmaYYErr
        (
            IOobject
            (
                "sigmaYYErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0)
        );

        volScalarField pErr
        (
            IOobject
            (
                "pErr",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0)
        );

        scalarField& sigmaXXErrI = Foam::primitiveFieldRef(sigmaXXErr);
        scalarField& sigmaXYErrI = Foam::primitiveFieldRef(sigmaXYErr);
        scalarField& sigmaYYErrI = Foam::primitiveFieldRef(sigmaYYErr);
        scalarField& pErrI = Foam::primitiveFieldRef(pErr);

        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");
        const volScalarField& p =
            mesh.lookupObject<volScalarField>("p");

        const symmTensorField& sigmaI = sigma.internalField();
        const scalarField& pI = p.internalField();
        const vectorField& C = mesh.C().internalField();

        forAll(sigmaI, cellI)
        {
            const vector curR(C[cellI].x(), C[cellI].y(), 0);
            const symmTensor curSigmaA = plateHoleStress(curR);
            const scalar curHydPressureA = plateHoleHydPressure(curR, mesh);

            sigmaXXErrI[cellI] = mag(sigmaI[cellI].xx() - curSigmaA.xx());
            sigmaXYErrI[cellI] = mag(sigmaI[cellI].xy() - curSigmaA.xy());
            sigmaYYErrI[cellI] = mag(sigmaI[cellI].yy() - curSigmaA.yy());
            pErrI[cellI] = mag(pI[cellI] - curHydPressureA);
        }

        Info<< "sigmaXXErr, max : " << gMax(sigmaXXErr) << endl;
        Info<< "sigmaXYErr, max : " << gMax(sigmaXYErr) << endl;
        Info<< "sigmaYYErr, max : " << gMax(sigmaYYErr) << endl;
        Info<< "pErr, max : " << gMax(pErr) << endl;

        sigmaXXErr.write();
        sigmaXYErr.write();
        sigmaYYErr.write();
        pErr.write();
    }
}

bool Foam::plateHoleAnalyticalSolution::writeData()
{
    const fvMesh& mesh = this->mesh();

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
                sI[cellI] = plateHoleStress(CI[cellI]);
            }

            if (cellDisplacement_)
            {
                aDI[cellI] = plateHoleDisplacement(CI[cellI], mesh);
            }
        }

        forAll(analyticalStress.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].type() != "empty")
            {
                symmTensorField& sP = boundaryFieldRef(analyticalStress)[patchI];
                vectorField& aDP = boundaryFieldRef(analyticalD)[patchI];
                const vectorField& CP = C.boundaryField()[patchI];

                forAll(sP, faceI)
                {
                    if (cellStress_)
                    {
                        sP[faceI] = plateHoleStress(CP[faceI]);
                    }

                    if (cellDisplacement_)
                    {
                        aDP[faceI] =
                            plateHoleDisplacement(CP[faceI], mesh);
                    }
                }
            }
        }

        // Write out the cell analytical field
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

            const vectorField& diffI = diff;
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
                sI[pointI] = plateHoleStress(points[pointI]);
            }

            if (pointDisplacement_)
            {
                aDI[pointI] = plateHoleDisplacement(points[pointI], mesh);
            }
        }

        // Write point analytical fields
        if (pointStress_)
        {
            Info<< "Writing analyticalPointStress"
                << nl << endl;
            analyticalStress.write();
        }

        if (pointDisplacement_)
        {
            Info<< "Writing analyticalPointDisplacement"
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

            const vectorField& diffI = diff;
            Info<< "    Norms: mean L1, mean L2, LInf: " << nl
                << "    " << gAverage(mag(diffI))
                << " " << Foam::sqrt(gAverage(magSqr(diffI)))
                << " " << gMax(mag(diffI))
                << nl << endl;
        }
    }

    if (pressureDisplacement_)
    {
        writePressureDisplacementData(mesh, pMesh);
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
    regionName_(dict.lookupOrDefault<word>("region", polyMesh::defaultRegion)),
    T_(readScalar(dict.lookup("farFieldTractionX"))),
    holeR_(readScalar(dict.lookup("holeRadius"))),
    E_(dict.found("E") ? readScalar(dict.lookup("E")) : -1.0),
    nu_(dict.found("nu") ? readScalar(dict.lookup("nu")) : -1.0),
    deriveMaterialProperties_
    (
        dict.lookupOrDefault<Switch>("deriveMaterialProperties", false)
     || !dict.found("E")
     || !dict.found("nu")
    ),
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
    pressureDisplacement_
    (
        dict.lookupOrDefault<Switch>("pressureDisplacement", false)
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;

    if (holeR_ < SMALL)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "holeRadius should be greater than 0!"
            << abort(FatalError);
    }

    if
    (
        !deriveMaterialProperties_
     && (E_ < SMALL || nu_ < SMALL)
    )
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


#ifdef OPENFOAM_NOT_EXTEND
bool Foam::plateHoleAnalyticalSolution::write()
{
    return false;
}
#endif

// ************************************************************************* //
