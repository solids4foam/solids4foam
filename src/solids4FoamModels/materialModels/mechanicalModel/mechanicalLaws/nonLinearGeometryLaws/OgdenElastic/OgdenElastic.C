/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "OgdenElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "eig3Field.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OgdenElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, OgdenElastic, nonLinGeomMechLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::OgdenElastic::OgdenElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu1_(dict.lookup("mu1")),
    mu2_(dict.lookup("mu2")),
    mu3_(dict.lookup("mu3")),
    alpha1_(dict.lookup("alpha1")),
    alpha2_(dict.lookup("alpha2")),
    alpha3_(dict.lookup("alpha3")),
    K_(dict.lookup("K"))
{
    Info<< "Material properties " << nl
        << "    mu1 = " << mu1_.value() << nl
        << "    mu2 = " << mu2_.value() << nl
        << "    mu3 = " << mu3_.value() << nl
        << "    alpha1 = " << alpha1_.value() << nl
        << "    alpha2 = " << alpha2_.value() << nl
        << "    alpha3 = " << alpha3_.value() << nl
        << "    K = "   << K_.value() << endl;

    // Compute lineaarised shear modulus
    mu() = mu1_ + mu2_ + mu3_;;

    // Update surface fields
    muf() = fvc::interpolate(mu());
    Kf()  = fvc::interpolate(K());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField> Foam::OgdenElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField((4.0/3.0)*mu() + K())
    );
}


void Foam::OgdenElastic::correct
(
    volSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    if (updateF(sigma, mu(), K()))
    {
        return;
    }

    // Calculate the Jacobian of the deformation gradient
    const volScalarField J(det(F()));

    // Update the hydrostatic stress
    updateSigmaHyd
    (
        0.5*K()*(pow(J, 2.0) - 1.0),
        (4.0/3.0)*mu() + K()
    );

    // Calculate the right Cauchy Green tensor
    const volSymmTensorField C(symm(F().T() & F()));

    // Eigen value field of C
    // We will store the eigen values in a vector instead of a diagTensor
    // because the tranform function is not definite for diagTensors on a wedge
    // boundary
    volVectorField lambda
    (
        IOobject
        (
            "eigenVal(" + C.name() + ")",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedVector("zero", C.dimensions(), vector::zero)
    );

    // Eigen vectors will be store in the rows i.e. the first eigen vector
    // is (eigenVec.xx() eigenVec.xy() eigenVec.xz())
    volTensorField eigVec
    (
        IOobject
        (
            "eigenVec(" + C.name() + ")",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    // Calculate eigen values and eigen vectors of C
    eig3Field(C, eigVec, lambda);

    // Initialize the deviatoric stress tensor
    volSymmTensorField s
    (
        IOobject
        (
            "s",
            C.time().timeName(),
            C.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C.mesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
    );

    // Calculate the deviatoric stress
    const vectorField& lambdaI = lambda.internalField();
    const tensorField& eigVecI = eigVec.internalField();
    symmTensorField& sI = s.primitiveFieldRef();
    forAll(sI, cellI)
    {
        // Calculate principal stresses individually for each cell
        const scalar p1
        (
            mu1_.value()
           *pow(max(sqrt(lambdaI[cellI].x()), VSMALL), alpha1_.value())
         +  mu2_.value()
           *pow(max(sqrt(lambdaI[cellI].x()), VSMALL), alpha2_.value())
         +  mu3_.value()
           *pow(max(sqrt(lambdaI[cellI].x()), VSMALL), alpha3_.value())
        );
        const scalar p2
        (
            mu1_.value()
           *pow(max(sqrt(lambdaI[cellI].y()), VSMALL), alpha1_.value())
         +  mu2_.value()
           *pow(max(sqrt(lambdaI[cellI].y()), VSMALL), alpha2_.value())
         +  mu3_.value()
           *pow(max(sqrt(lambdaI[cellI].y()), VSMALL), alpha3_.value())
        );
        const scalar p3
        (
            mu1_.value()
           *pow(max(sqrt(lambdaI[cellI].z()), VSMALL), alpha1_.value())
         +  mu2_.value()
           *pow(max(sqrt(lambdaI[cellI].z()), VSMALL), alpha2_.value())
         +  mu3_.value()
            *pow(max(sqrt(lambdaI[cellI].z()), VSMALL), alpha3_.value())
        );

        // Set up the pricipal stress tensor
        const symmTensor prinStress
        (
            p1, 0, 0,
                p2, 0,
                    p3
        );

        // Rotate back to Cauchy stress for each cell
        sI[cellI] = transform(eigVecI[cellI].T(), prinStress);
    }

    forAll(s.boundaryField(), patchI)
    {
        const vectorField& lambdaP = lambda.boundaryField()[patchI];
        const tensorField& eigVecP = eigVec.boundaryField()[patchI];

#ifdef FOAMEXTEND
        symmTensorField& sP = s.boundaryField()[patchI];
#else
        symmTensorField& sP = s.boundaryFieldRef()[patchI];
#endif

        forAll(sP, cellI)
        {
            // Calculate principal stresses individually for each cell
            const scalar p1
            (
                mu1_.value()
               *pow(max(sqrt(lambdaP[cellI].x()), VSMALL), alpha1_.value())
             +  mu2_.value()
               *pow(max(sqrt(lambdaP[cellI].x()), VSMALL), alpha2_.value())
             +  mu3_.value()
               *pow(max(sqrt(lambdaP[cellI].x()), VSMALL), alpha3_.value())
            );
            const scalar p2
            (
                mu1_.value()
               *pow(max(sqrt(lambdaP[cellI].y()), VSMALL), alpha1_.value())
             +  mu2_.value()
               *pow(max(sqrt(lambdaP[cellI].y()), VSMALL), alpha2_.value())
             +  mu3_.value()
               *pow(max(sqrt(lambdaP[cellI].y()), VSMALL), alpha3_.value())
            );
            const scalar p3
            (
                mu1_.value()
               *pow(max(sqrt(lambdaP[cellI].z()), VSMALL), alpha1_.value())
             +  mu2_.value()
               *pow(max(sqrt(lambdaP[cellI].z()), VSMALL), alpha2_.value())
             +  mu3_.value()
               *pow(max(sqrt(lambdaP[cellI].z()), VSMALL), alpha3_.value())
            );

            // Set up the pricipal stress tensor
            const symmTensor prinStress
            (
                p1, 0, 0,
                    p2, 0,
                        p3
            );

            // Rotate back to Cauchy stress for each cell
            sP[cellI] = transform(eigVecP[cellI].T(), prinStress);
        }
    }

    s.correctBoundaryConditions();

    // Calculate the Cauchy stress
    sigma =
        (1.0/J)
       *(
           dev(s - (mu1_ + mu2_ + mu3_)*I) + sigmaHyd()*I + symm(F() & sigma0() & F().T())
        );
}


void Foam::OgdenElastic::correct(surfaceSymmTensorField& sigma)
{
    notImplemented("OgdenElastic::correct(surfaceSymmTensorField& sigma)");
}


// ************************************************************************* //
