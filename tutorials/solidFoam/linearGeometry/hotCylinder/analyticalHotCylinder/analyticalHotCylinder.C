/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Generate analytical solution for a thick-walled cylinder with a
    temperature gradient.
    Temperature field T and stress field sigma and generated.
    Based on solution outlined in Timoshenko, Theory of Elasticity.

Author
    philip.cardiff@ucd.ie

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    runTime++;

    Info<< "Writing analytical solution for a plain strain cylinder "
        << "with concentric hole,\nwhere"
        << "\n\tinner radius = 0.5"
        << "\n\touter radius = 0.7"
        << "\n\tinner temperature = 100"
        << "\n\touter temperature = 0"
        << "\n\tinner pressure = 0"
        << "\n\touter pressure = 0"
        << "\n\tE = 200e9"
        << "\n\tu = 0.3"
        << "\n\talpha = 1e-5"
        << nl << endl;

    //- inner and outer radii and temperatures
    scalar a = 0.5;
    scalar b = 0.7;
    scalar Ti = 100;
    scalar To = 0;

    //- mechanical and thermal properties
    scalar E = 200e9;
    scalar nu = 0.3;
    scalar alpha = 1e-5;

    const volVectorField& C = mesh.C();

    //- radial coordinate
    volScalarField radii
    (
        sqrt
        (
            sqr(C.component(vector::X))
          + sqr(C.component(vector::Y))
        )/dimensionedScalar("one", dimLength, 1)
    );

    const scalarField& rIn = radii.internalField();

    Info << "Writing analytical termpature field" << endl;
    //- create T field
    volScalarField T
    (
        IOobject
        (
            "analyticalT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((Ti - To)/Foam::log(b/a))*Foam::log(b/radii)
    );
    T.write();

    //- create sigma field
    Info << "\nWriting analytical sigmaR field" << endl;
    volScalarField sigmaR
    (
        IOobject
        (
            "sigmaR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((alpha*E*(Ti - To))/(2*(1 - nu)*Foam::log(b/a)))*
        (
            -Foam::log(b/radii)
          - (sqr(a)/(sqr(b) - sqr(a)))*(1 - sqr(b)/sqr(radii))*Foam::log(b/a)
        )
    );
    sigmaR.write();


    Info << "\nWriting analytical sigmaTheta field" << endl;
    volScalarField sigmaTheta
    (
        IOobject
        (
            "sigmaTheta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((alpha*E*(Ti - To))/(2*(1 - nu)*Foam::log(b/a)))*
        (
            1 - Foam::log(b/radii)
          - (sqr(a)/(sqr(b) - sqr(a)))*(1 + sqr(b)/sqr(radii))*Foam::log(b/a)
        )
    );
    sigmaTheta.write();

    Info << "\nWriting analytical sigmaZ field" << endl;
    volScalarField sigmaZ
    (
        IOobject
        (
            "sigmaZ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        // Timoshenko says this but I am not sure I am not sure the BCs in
        // the z direction
        // ((alpha*E*(Ti - To))/(2*(1 - nu)*Foam::log(b/a)))*
        // (1 - 2*Foam::log(b/radii) - ( 2*sqr(a)/(sqr(b) - sqr(a)))*Foam::log(b/a));
        0.3*(sigmaR + sigmaTheta) - E*alpha*(T)
    );
    sigmaZ.write();

    //- create theta field
    volScalarField yOverX
    (
        "yOverX",
        Foam::max
        (
            scalar(-1),
            Foam::min
            (
                scalar(1),
                mesh.C().component(vector::Y)/
                stabilise
                (
                    mesh.C().component(vector::X),
                    dimensionedScalar("small", dimLength, SMALL)
                )
            )
        )
    );

    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Foam::atan(yOverX)
    );

    //- rotation matrix to convert polar stresses to cartesian
    volTensorField rotMat
    (
        IOobject
        (
            "rotMat",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    tensorField& rotMatIn = rotMat.internalField();
    const scalarField tIn = theta.internalField();

    forAll (rotMatIn, celli)
    {
        const scalar& t = tIn[celli];

        rotMatIn[celli] =
            tensor
            (
                Foam::cos(t),  Foam::sin(t), 0,
               -Foam::sin(t),  Foam::cos(t), 0,
                0, 0, 1
            );
    }


    forAll (rotMat.boundaryField(), patchi)
    {
        forAll (rotMat.boundaryField()[patchi], facei)
        {
            const scalar& t = theta.boundaryField()[patchi][facei];

            rotMat.boundaryField()[patchi][facei] =
                tensor
                (
                    Foam::cos(t),  Foam::sin(t), 0,
                   -Foam::sin(t),  Foam::cos(t), 0,
                    0, 0, 1
                );
        }
    }

    volSymmTensorField sigma
    (
        IOobject
        (
            "analyticalSigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    {
        symmTensorField& sigmaIn = sigma.internalField();

        const scalarField& rIn = sigmaR.internalField();
        const scalarField& tIn = sigmaTheta.internalField();
        const scalarField& zIn = sigmaZ.internalField();

        forAll (sigmaIn, celli)
        {
            symmTensor sigmaCart
            (
                rIn[celli], 0, 0,
                tIn[celli], 0,
                zIn[celli]
            );

            const tensor& rot = rotMatIn[celli];

            sigmaIn[celli] = symm(rot.T() & sigmaCart & rot);

            // for general 2-D plain strain problems, the axial stress is:
            // (which is not equal to the solution by Timoshenko... hmmmnn)
            //       sigmaIn[celli][symmTensor::ZZ] =
            //           0.3*(sigmaIn[celli][symmTensor::XX]
            //         + sigmaIn[celli][symmTensor::YY])
            //         - E*alpha*(T.internalField()[celli]);
        }
    }

    forAll (sigma.boundaryField(), patchi)
    {
        symmTensorField& pSigma = sigma.boundaryField()[patchi];
        const scalarField& pR = sigmaR.boundaryField()[patchi];
        const scalarField& pT = sigmaTheta.boundaryField()[patchi];
        const scalarField& pZ = sigmaZ.boundaryField()[patchi];

        const tensorField pRot = rotMat.boundaryField()[patchi];

        forAll (pSigma, facei)
        {
            const tensor& rot = pRot[facei];

            symmTensor sigmaCart
            (
                pR[facei], 0, 0,
                pT[facei], 0,
                pZ[facei]
            );

            pSigma[facei] = symm(rot.T() & sigmaCart & rot);
        }
    }

    Info << "\nWriting analytical sigma tensor" << endl;
    sigma.write();

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
