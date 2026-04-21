/*---------------------------------------------------------------------------* \
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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "setPlateHoleBC.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "surfaceFields.H"

#include "ValuePointPatchField.H"
#include "tractionPressureDisplacementFvPatchVectorField.H"
#include "solidTractionFvPatchVectorField.H"
#include "linearElastic.H"
#include "lookupSolidModel.H"
#include "OFstream.H"

#include "linearElastic.H"
#include "pdNeoHookeanElastic.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setPlateHoleBC, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setPlateHoleBC,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setPlateHoleBC::planeStress() const
{
    // Reference to fvMesh
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

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
    else
    {
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
}

Foam::symmTensor Foam::setPlateHoleBC::plateHoleStress(const vector& C) const
{
    tensor sigma = tensor::zero;

    // scalar T = 1e6;
    // scalar a = 1;

    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    scalar theta = atan2(C.y(), C.x());

    sigma.xx() =
        T_
       *(
           1.0
         - (sqr(a_)/sqr(r))*(3*cos(2*theta)/2 + cos(4*theta))
         + (3*pow(a_,4)/(2*pow(r,4)))*cos(4*theta)
        );

    sigma.yy() =
        T_
       *(
         - (sqr(a_)/sqr(r))*(cos(2*theta)/2 - cos(4*theta))
         - (3*pow(a_,4)/(2*pow(r,4)))*cos(4*theta)
        );

    sigma.xy() =
        T_
       *(
         - (sqr(a_)/sqr(r))*(sin(2*theta)/2 + sin(4*theta))
         + (3*pow(a_,4)/(2*pow(r,4)))*sin(4*theta)
        );

    sigma.yx() = sigma.xy();

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}


Foam::vector Foam::setPlateHoleBC::plateHoleDisplacement(const vector& C) const
{
    vector displacement = vector::zero;

    // scalar T = 1e6;
    // scalar a = 1;

    // Reference to fvMesh
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(mesh);

    const PtrList<mechanicalLaw>& laws = solMod.mechanical();
    scalar mu = 0;
    if (isA<linearElastic>(laws[0]))
    {
        mu = refCast<const linearElastic>(laws[0]).mu().value();
    }
    else if (isA<pdNeoHookeanElastic>(laws[0]))
    {
        mu =
            refCast<const pdNeoHookeanElastic>
            (
                laws[0]
            ).mu().value();
    }
    else
    {
        FatalErrorIn("setPlateHoleBC::plateHoleDisplacement(const vector)")
            << "Mehanical law: " << laws[0].type()
            << " is not supported."
            << abort(FatalError);
    }

    scalar K = solMod.mechanical().bulkModulus()()[0];

    scalar nu = 0.5;
    if (K+SMALL < GREAT)
    {
        nu = (3*K - 2*mu)/(2*(3*K + mu));

        if (planeStress())
        {
            nu = (K - mu)/(K + mu);
        }
    }

    // Kolosov constnt
    scalar kappa = 3-4*nu;
    if (planeStress())
    {
        kappa = (3.0-nu)/(1.0+nu);
    }

    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    scalar theta = atan2(C.y(), C.x());

    displacement.x() =
        (a_*T_/(8*mu))
       *(
           (r/a_)*(kappa+1)*cos(theta)
         + (2*a_/r)*((1+kappa)*cos(theta) + cos(3*theta))
         - (2*pow(a_,3)/pow(r,3))*cos(3*theta)
        );

    displacement.y() =
        (a_*T_/(8*mu))
       *(
           (r/a_)*(kappa-3)*sin(theta)
         + (2*a_/r)*((1-kappa)*sin(theta) + sin(3*theta))
         - (2*pow(a_,3)/pow(r,3))*sin(3*theta)
        );

    return displacement;
}

Foam::scalar Foam::setPlateHoleBC::plateHoleHydPressure
(
    const vector& C
) const
{
    scalar hydPressure = 0;

    // Reference to fvMesh
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(mesh);

    const PtrList<mechanicalLaw>& laws = solMod.mechanical();
    scalar mu = 0;
    if (isA<linearElastic>(laws[0]))
    {
        mu = refCast<const linearElastic>(laws[0]).mu().value();
    }
    else if (isA<pdNeoHookeanElastic>(laws[0]))
    {
        mu =
            refCast<const pdNeoHookeanElastic>
            (
                laws[0]
            ).mu().value();
    }
    else
    {
        FatalErrorIn("setPlateHoleBC::plateHoleDisplacement(const vector)")
            << "Mehanical law: " << laws[0].type()
            << " is not supported."
            << abort(FatalError);
    }

    scalar K = solMod.mechanical().bulkModulus()()[0];

    scalar nu = 0.5;
    if (K+SMALL < GREAT)
    {
        nu = (3*K - 2*mu)/(2*(3*K + mu));

        if (planeStress())
        {
            nu = (K - mu)/(K + mu);
        }
    }

    symmTensor sigma = plateHoleStress(C);

    scalar sigmaZZ = nu*(sigma.xx() + sigma.yy());

    hydPressure = -(sigma.xx() + sigma.yy() + sigmaZZ)/3;

    return hydPressure;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setPlateHoleBC::setPlateHoleBC
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    a_(readScalar(dict.lookup("radius"))),
    T_(readScalar(dict.lookup("traction")))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setPlateHoleBC::start()
{
    setBC();

    return true;
}


// #if FOAMEXTEND > 40
// bool Foam::setPlateHoleBC::execute(const bool forceWrite)
// #else
// bool Foam::setPlateHoleBC::execute()
// #endif

bool Foam::setPlateHoleBC::execute(const bool forceWrite)
{
    calcError();

    return false;
}

void Foam::setPlateHoleBC::setBC()
{
    Info << "Update boundary conditions" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(mesh);

    word DName("D");

    if (solMod.incremental())
    {
        DName = word("DD");
    }

    // Info << "DName: " << DName << endl;

    volVectorField& D =
        const_cast<volVectorField&>(mesh.lookupObject<volVectorField>(DName));

    // scalar a = 1.0;

    forAll(D.boundaryField(), patchI)
    {
        if
        (
            D.boundaryField()[patchI].type()
         == tractionPressureDisplacementFvPatchVectorField::typeName
        )
        {
            if (mesh.boundary()[patchI].name() != "hole")
            {
                tractionPressureDisplacementFvPatchVectorField& patchD =
                    refCast<tractionPressureDisplacementFvPatchVectorField>
                    (
                        D.boundaryField()[patchI]
                    );

                vectorField nf = patchD.patch().nf();
                scalarField magSf = patchD.patch().magSf();

                const vectorField& Cf = mesh.Cf().boundaryField()[patchI];

                forAll(patchD.traction(), faceI)
                {
                    vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
                    vector curN = nf[faceI];

                    // if (mesh.boundary()[patchI].name() == "hole")
                    // {
                    //     curC /= mag(curC);
                    //     curC *= a_;

                    //     curN = -curC/mag(curC);
                    // }

                    patchD.traction()[faceI] = (curN & plateHoleStress(curC));
                }

                Info << mesh.boundary()[patchI].name() << ": "
                     << sum(patchD.traction()*magSf) << endl;
            }
        }
        else if
        (
            D.boundaryField()[patchI].type()
         == solidTractionFvPatchVectorField::typeName
        )
        {
            solidTractionFvPatchVectorField& patchD =
                refCast<solidTractionFvPatchVectorField>
                (
                    D.boundaryField()[patchI]
                );

            vectorField nf = patchD.patch().nf();

            const vectorField& Cf = mesh.Cf().boundaryField()[patchI];

            forAll(patchD.traction(), faceI)
            {
                vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
                vector curN = nf[faceI];

                if (mesh.boundary()[patchI].name() == "hole")
                {
                    curC /= mag(curC);
                    curC *= a_;

                    curN = -curC/mag(curC);

                    // patchD.traction()[faceI] = vector::zero;
                    // patchD.traction()[faceI] =
                    //     (curN & plateHoleStress(curC));
                }
                // else
                // {
                //     patchD.traction()[faceI] = (curN & plateHoleStress(curC));
                // }

                patchD.traction()[faceI] = (nf[faceI] & plateHoleStress(curC));
            }

            // Info << mesh.boundary()[patchI].name() << ", "
            //      << average(mag(patchD.traction())) << endl;
        }
    }
}


void Foam::setPlateHoleBC::calcError() const
{
    Info << "Calculating errors" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Cell displacement
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
            dimensionedVector("0", dimLength, vector::zero)
        );

        const volVectorField& D =
            mesh.lookupObject<volVectorField>("D");

        const vectorField& DI = D.internalField();

        const vectorField& C = mesh.C().internalField();

        vectorField& DError = Derror.internalField();

        forAll(DError, cellI)
        {
            vector curR = vector(C[cellI].x(), C[cellI].y(), 0);
            vector curDa = plateHoleDisplacement(curR);

            DError[cellI] = (DI[cellI] - curDa);
        }

        Info << "DError, max : " << gMax(mag(DError)) << endl;
        Info << "DError, avg : " << gAverage(mag(DError)) << endl;
        Info << "DError, L2 : " << sqrt(sum(magSqr(DError))/C.size()) << endl;

        Derror.write();
    }

    // Point displacement
    if (true)
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>("pointD");

        const vectorField& pointDI = pointD.internalField();

        const vectorField& C = mesh.points();

        scalarField pointDError(pointDI.size(), 0);

        forAll(pointDError, pointI)
        {
            vector curR = vector(C[pointI].x(), C[pointI].y(), 0);
            vector curDa = plateHoleDisplacement(curR);

            pointDError[pointI] = mag(pointDI[pointI] - curDa);
        }

        Info << "pointDError, max : " << gMax(pointDError) << endl;
    }

    // Face stress
    if (true)
    {
        const surfaceSymmTensorField& sigmaf =
            mesh.lookupObject<surfaceSymmTensorField>("sigmaf");

        const surfaceVectorField& Cf = mesh.Cf();

        surfaceSymmTensorField sigmafA
        (
            IOobject
            (
                "sigmafA",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
        );
        symmTensorField& sigmafAI = sigmafA.internalField();


        forAll(sigmafAI, faceI)
        {
            vector curR = vector(Cf[faceI].x(), Cf[faceI].y(), 0);
            sigmafAI[faceI] = plateHoleStress(curR);
        }

        forAll(sigmafA.boundaryField(), patchI)
        {
            symmTensorField& pSigmafA =
                sigmafA.boundaryField()[patchI];

            const vectorField& pCf = Cf.boundaryField()[patchI];

            forAll(pSigmafA, faceI)
            {
                vector curR = vector(pCf[faceI].x(), pCf[faceI].y(), 0);
                pSigmafA[faceI] = plateHoleStress(curR);
            }
        }

        surfaceSymmTensorField sigmafErr = sigmaf - sigmafA;

//         scalar maxSigmafXxErr = max(mag(sigmafErr.component(symmTensor::XX))).value();
//         scalar maxSigmafXyErr = max(mag(sigmafErr.component(symmTensor::XY))).value();
//         scalar maxSigmafYyErr = max(mag(sigmafErr.component(symmTensor::YY))).value();
//         Info << "sigmafXxErr, max : " << maxSigmafXxErr << endl;
//         Info << "sigmafXyErr, max : " << maxSigmafXyErr << endl;
//         Info << "sigmafYyErr, max : " << maxSigmafYyErr << endl;

        scalar maxSigmafIXxErr = max(mag(sigmafErr.internalField().component(symmTensor::XX)));
        scalar maxSigmafIXyErr = max(mag(sigmafErr.internalField().component(symmTensor::XY)));
        scalar maxSigmafIYyErr = max(mag(sigmafErr.internalField().component(symmTensor::YY)));

        Info << "sigmafIXxErr, max : " << maxSigmafIXxErr << endl;
        Info << "sigmafIXyErr, max : " << maxSigmafIXyErr << endl;
        Info << "sigmafIYyErr, max : " << maxSigmafIYyErr << endl;

        forAll(sigmafErr.boundaryField(), patchI)
        {
            if (mesh.boundary()[patchI].name() == "hole")
            {
                scalar maxSigmafXxErrHole =
                    max(mag(sigmafErr.boundaryField()[patchI].component(symmTensor::XX)));
                scalar maxSigmafXyErrHole =
                    max(mag(sigmafErr.boundaryField()[patchI].component(symmTensor::XY)));
                scalar maxSigmafYyErrHole =
                    max(mag(sigmafErr.boundaryField()[patchI].component(symmTensor::YY)));

                Info << "sigmafXxErrHole, max : " << maxSigmafXxErrHole << endl;
                Info << "sigmafXyErrHole, max : " << maxSigmafXyErrHole << endl;
                Info << "sigmafYyErrHole, max : " << maxSigmafYyErrHole << endl;
            }
        }
    }

    if (time_.outputTime())
    {
        // Stress error field for plate with hole case
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
            dimensionedScalar("zero", dimless, 0)
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
            dimensionedScalar("zero", dimless, 0)
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
            dimensionedScalar("zero", dimless, 0)
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
            dimensionedScalar("zero", dimless, 0)
        );

        scalarField& sigmaXXErrI = sigmaXXErr.internalField();
        scalarField& sigmaXYErrI = sigmaXYErr.internalField();
        scalarField& sigmaYYErrI = sigmaYYErr.internalField();
        scalarField& pErrI = pErr.internalField();

        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        const volScalarField& p =
            mesh.lookupObject<volScalarField>("p");

        const scalarField& pI = p.internalField();
        const symmTensorField& sigmaI = sigma.internalField();
        const vectorField& C = mesh.C().internalField();

//         scalar maxSigmaXX = max(mag(sigma.component(symmTensor::XX))).value();
//         scalar maxSigmaXY = max(mag(sigma.component(symmTensor::XY))).value();
//         scalar maxSigmaYY = max(mag(sigma.component(symmTensor::YY))).value();

        forAll(sigmaI, cellI)
        {
            vector curR = vector(C[cellI].x(), C[cellI].y(), 0);
            symmTensor curSigmaA = plateHoleStress(curR);
            scalar curHydPressureA = plateHoleHydPressure(curR);

            sigmaXXErrI[cellI] =
                mag(sigmaI[cellI].xx() - curSigmaA.xx());
//                 mag(sigmaI[cellI].xx() - curSigmaA.xx())/maxSigmaXX;

            sigmaXYErrI[cellI] =
                mag(sigmaI[cellI].xy() - curSigmaA.xy());
//                 mag(sigmaI[cellI].xy() - curSigmaA.xy())/maxSigmaXY;

            sigmaYYErrI[cellI] =
                mag(sigmaI[cellI].yy() - curSigmaA.yy());
//                 mag(sigmaI[cellI].yy() - curSigmaA.yy())/maxSigmaYY;

            pErrI[cellI] =
                mag(pI[cellI] - curHydPressureA);
        }

        Info << "sigmaXXErr, max : " << gMax(sigmaXXErr) << endl;
        Info << "sigmaXYErr, max : " << gMax(sigmaXYErr) << endl;
        Info << "sigmaYYErr, max : " << gMax(sigmaYYErr) << endl;
        Info << "pErr, max : " << gMax(pErr) << endl;

        Info << "sigmaXXErr, avg : " << gAverage(sigmaXXErr) << endl;
        Info << "sigmaXYErr, avg : " << gAverage(sigmaXYErr) << endl;
        Info << "sigmaYYErr, avg : " << gAverage(sigmaYYErr) << endl;

        sigmaXXErr.write();
        sigmaXYErr.write();
        sigmaYYErr.write();
        pErr.write();

        // Point displacement
        if (true)
        {
            pointMesh pMesh(mesh);

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
                dimensionedVector("0", dimLength, vector::zero)
            );

            const pointVectorField& pointD =
                mesh.lookupObject<pointVectorField>("pointD");

            const vectorField& pointDI = pointD.internalField();

            const vectorField& C = mesh.points();

            scalarField pointDError(pointDI.size(), 0);

            forAll(pointDError, pointI)
            {
                vector curR = vector(C[pointI].x(), C[pointI].y(), 0);
                vector curDa = plateHoleDisplacement(curR);

                pointDerror.internalField()[pointI] = (pointDI[pointI] - curDa);
            }

            pointDerror.write();
        }

        // Write stress and displacement to file
        {
            OFstream fileLeft("left.dat");
            forAll(mesh.boundary(), patchI)
            {
                if (mesh.boundary()[patchI].name() == word("left"))
                {
                    const symmTensorField& sigmaLeft =
                        sigma.boundaryField()[patchI];
                    const vectorField& CLeft =
                        mesh.C().boundaryField()[patchI];
                    const volVectorField& D =
                        mesh.lookupObject<volVectorField>("D");
                    const vectorField& DLeft =
                        D.boundaryField()[patchI];

                    //header
                    fileLeft << "#x" << " " << "y" << " "
                             << "sigma_xx" << " "
                             << "sigma_xx analitic" << " "
                             << "sigma_yy" << " "
                             << "sigma_yy analitic" << " "
                             << "sigma_xy" << " "
                             << "sigma_xy analitic" << " "
                             << "D_x"<< " "
                             << "D_x_analitic"<< " "
                             << "D_y" << " "
                             << "D_y_analitic" <<endl;

                    forAll(sigmaLeft, faceI)
                    {
                        vector curR =
                            vector(CLeft[faceI].x(), CLeft[faceI].y(), 0);
                        symmTensor curSigmaA = plateHoleStress(curR);
                        vector curD = plateHoleDisplacement(curR);

                        fileLeft << CLeft[faceI].x() << " "
                                 << CLeft[faceI].y() << " "
                                 << sigmaLeft[faceI].xx() << " "
                                 << curSigmaA.xx() << " "
                                 << sigmaLeft[faceI].yy() << " "
                                 << curSigmaA.yy() << " "
                                 << sigmaLeft[faceI].xy() << " "
                                 << curSigmaA.xy() << " "
                                 << DLeft[faceI].x() << " "
                                 << curD.x() << " "
                                 << DLeft[faceI].y() << " "
                                 << curD.y() << endl;
                    }
                }
            }
        }

        {
            OFstream fileDown("down.dat");
            forAll(mesh.boundary(), patchI)
            {
                if (mesh.boundary()[patchI].name() == word("down"))
                {
                    const symmTensorField& sigmaLeft =
                        sigma.boundaryField()[patchI];
                    const vectorField& CLeft =
                        mesh.C().boundaryField()[patchI];
                    const volVectorField& D =
                        mesh.lookupObject<volVectorField>("D");
                    const vectorField& DLeft =
                        D.boundaryField()[patchI];

                    //header
                    fileDown << "#x" << " " << "y" << " "
                             << "sigma_xx" << " "
                             << "sigma_xx analitic" << " "
                             << "sigma_yy" << " "
                             << "sigma_yy analitic" << " "
                             << "sigma_xy" << " "
                             << "sigma_xy analitic" << " "
                             << "D_x"<< " "
                             << "D_x_analitic"<< " "
                             << "D_y" << " "
                             << "D_y_analitic" <<endl;

                    forAll(sigmaLeft, faceI)
                    {
                        vector curR =
                            vector(CLeft[faceI].x(), CLeft[faceI].y(), 0);
                        symmTensor curSigmaA = plateHoleStress(curR);
                        vector curD = plateHoleDisplacement(curR);

                        fileDown << CLeft[faceI].x() << " "
                                 << CLeft[faceI].y() << " "
                                 << sigmaLeft[faceI].xx() << " "
                                 << curSigmaA.xx() << " "
                                 << sigmaLeft[faceI].yy() << " "
                                 << curSigmaA.yy() << " "
                                 << sigmaLeft[faceI].xy() << " "
                                 << curSigmaA.xy() << " "
                                 << DLeft[faceI].x() << " "
                                 << curD.x() << " "
                                 << DLeft[faceI].y() << " "
                                 << curD.y() << endl;
                    }
                }
            }
        }

        {
            OFstream fileNorm("norm.dat");

            const fvMesh& mesh =
                time_.lookupObject<fvMesh>(regionName_);

            const volVectorField& D =
                mesh.lookupObject<volVectorField>("D");

            const vectorField& DI = D.internalField();
            const vectorField& C = mesh.C().internalField();
            // int meshSize = mesh.cells().size();
            label avgCellSize =
                gAverage(1.0/mesh.deltaCoeffs().internalField());
            vectorField DError = D.internalField();

            forAll(DError, cellI)
            {
                vector curR = vector(C[cellI].x(), C[cellI].y(), 0);
                vector curDa = plateHoleDisplacement(curR);
                DError[cellI] = (DI[cellI] - curDa);
            }

            scalarField sigmaXXErrI = sigmaXXErr.internalField();
            scalarField sigmaXYErrI = sigmaXYErr.internalField();
            scalarField sigmaYYErrI = sigmaYYErr.internalField();


            const volSymmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>("sigma");

            const symmTensorField& sigmaI = sigma.internalField();

            forAll(sigmaI, cellI)
            {
                vector curR = vector(C[cellI].x(), C[cellI].y(), 0);
                symmTensor curSigmaA = plateHoleStress(curR);

                sigmaXXErrI[cellI] =
                    mag(sigmaI[cellI].xx() - curSigmaA.xx());

                sigmaXYErrI[cellI] =
                    mag(sigmaI[cellI].xy() - curSigmaA.xy());

                sigmaYYErrI[cellI] =
                    mag(sigmaI[cellI].yy() - curSigmaA.yy());

            }

            fileNorm << "# avg cell size" << "\t"
                     << "D infinity norm" << "\t"
                     << "D L1 norm" << "\t"
                     << "D L2 norm" << "\t"
                     << "sigma_xx infinity norm" << "\t"
                     << "sigma_xx L1 norm" << "\t"
                     << "sigma_xx L2 norm" <<endl;
            fileNorm << avgCellSize << " "
                     << gMax(mag(DError)) << " "
                     << average(mag(DError)) << " "
                     << sqrt(average(magSqr(DError))) << " "
                     << gMax(mag(sigmaXXErrI)) << " "
                     << sum(mag(sigmaXXErrI)) << " "
                     << sqrt(sum(magSqr(sigmaXXErrI))) <<  endl;
            //}
         }
    }
}


bool Foam::setPlateHoleBC::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
