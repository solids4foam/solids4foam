/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "neoHookeanElasticMisesPlasticRubin.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticMisesPlasticRubin, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElasticMisesPlasticRubin, dictionary
    );
} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::neoHookeanElasticMisesPlasticRubin::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn
        (
            "void Foam::neoHookeanElasticMisesPlasticRubin::makeRelF()"
        )   << "pointer already set" << abort(FatalError);
    }

    relFPtr_ =
        new volTensorField
        (
            IOobject
            (
                "lawRelF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::neoHookeanElasticMisesPlasticRubin::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}


void Foam::neoHookeanElasticMisesPlasticRubin::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticMisesPlasticRubin::makeJ()")
            << "pointer already set" << abort(FatalError);
    }

    JPtr_ =
        new volScalarField
        (
            IOobject
            (
                "lawJ",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );

    // Store the old-time
    JPtr_->oldTime();
}


Foam::volScalarField& Foam::neoHookeanElasticMisesPlasticRubin::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}


void Foam::neoHookeanElasticMisesPlasticRubin::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticMisesPlasticRubin::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ =
        new volTensorField
        (
            IOobject
            (
                "lawF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );

    // Store the old-time
    FPtr_->oldTime();
}


Foam::volTensorField& Foam::neoHookeanElasticMisesPlasticRubin::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElasticMisesPlasticRubin::Ibar
(
    const volSymmTensorField& bEbarTrial
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar o denforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<volScalarField> tIbar
    (
        new volScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& Ibar = tIbar();

    // Calculate deviatoric component of bEbarTrial
    const volSymmTensorField devBEbarTrial = dev(bEbarTrial);

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarTrialI = devBEbarTrial.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarTrialI[cellI]);
        const scalar dotprod = devBEbarTrialI[cellI] && devBEbarTrialI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryField()[patchI];
            const symmTensorField& devBEbarTrialP =
                devBEbarTrial.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarTrialP[faceI]);
                const scalar dotprod =
                    devBEbarTrialP[faceI] && devBEbarTrialP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElasticMisesPlasticRubin::neoHookeanElasticMisesPlasticRubin
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    relFPtr_(NULL),
    JPtr_(NULL),
    FPtr_(NULL),
    stressPlasticStrainSeries_(dict),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    Je_
    (
        IOobject
        (
            "Je",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    DEpsilonPEq_
    (
        IOobject
        (
            "DEpsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    )
    // plasticN_
    // (
    //     IOobject
    //     (
    //         "plasticN",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    // ),
{
    // Force storage of old time for adjustable time-step calculations
    bEbar_.oldTime();
    Je_.oldTime();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E/(2.0*(1.0 + nu));

        // Set the bulk modulus
        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn
        (
            "neoHookeanElasticMisesPlasticRubin::"
            "neoHookeanElasticMisesPlasticRubin::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElasticMisesPlasticRubin::~neoHookeanElasticMisesPlasticRubin()
{
    deleteDemandDrivenData(relFPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(FPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlasticRubin::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "lawRho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlasticRubin::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void Foam::neoHookeanElasticMisesPlasticRubin::correct
(
    volSymmTensorField& sigma
)
{
    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the relative deformation gradient
        relF() = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF() & F().oldTime();
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Update the total deformation gradient
        F() = I + gradD.T();

        // Update the relative deformation gradient
        relF() = F() & inv(F().oldTime());
    }

    // Update the Jacobian of the total deformation gradient
    J() = det(F());

    // Calculate the relative Jacobian
    const volScalarField relJ = J()/J().oldTime();

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    const volTensorField relFbar = pow(relJ, -1.0/3.0)*relF();

    // Update bE trial
    bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Calculate trial deviatoric stress
    const volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

    // Calculate sigmaStar
    // Check: is this an alternative definition to sigmaEq?
    const volSymmTensorField devBEbarTrial = dev(bEbarTrial_);
    const volScalarField sigmaStar =
        mu_*sqrt((3.0/2.0)*tr(devBEbarTrial & devBEbarTrial));

    // Calculate yield function
    const volScalarField fTrial = sigmaStar/sigmaY_;

    // Calculate plastic multiplier
    lambda_ = 1.0/max(1.0, fTrial);

    // Calculate bEBar
    bEbar_.storePrevIter();
    bEbar_ = lambda_*dev(bEbarTrial_) + Ibar(bEbarTrial_)*I;

    // Calcualte Je
    Je_ =
        relJ*Je_.oldTime()
       *exp
        (
           (3.0/2.0)*((1.0 - lambda_)/lambda_)
          *(
              1.0 - (9.0/(tr(bEbar_)*tr(inv(bEbar_))))
          )
        );

    // Calculate the deviatoric stress
    const volSymmTensorField s = mu_*dev(bEbar_);

    // Update hydrostatic stress (negative of pressure)
    sigmaHyd_ = 0.5*K_*(pow(J(), 2.0) - 1.0);

    // Update the Cauchy stress
    sigma = (1.0/J())*(sigmaHyd_*I + s);
}


void Foam::neoHookeanElasticMisesPlasticRubin::correct
(
    surfaceSymmTensorField& sigma
)
{
    FatalErrorIn
    (
        "void Foam::neoHookeanElasticMisesPlasticRubin::correct\n"
        "(\n"
        "    surfaceSymmTensorField& sigma\n"
        ")"
    )   << "Not implemented" << abort(FatalError);
}


Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin::residual()
{
    // Calculate residual based on change in plastic strain increment
    {
        return
            gMax
            (
                mag
                (
                    bEbar_.internalField()
                  - bEbar_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(bEbar_.prevIter().internalField()));
    }
}


void Foam::neoHookeanElasticMisesPlasticRubin::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;

    //Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    // epsilonPEq_ += DEpsilonPEq_;
    // epsilonPEqf_ += DEpsilonPEqf_;
    // epsilonP_ += DEpsilonP_;
    // epsilonPf_ += DEpsilonPf_;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.internalField(), cellI)
    {
        if (lambda_.internalField()[cellI] < (1.0 - SMALL))
        {
            activeYield_.internalField()[cellI] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.internalField()[cellI] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchI)
    {
        if (!activeYield_.boundaryField()[patchI].coupled())
        {
            forAll(activeYield_.boundaryField()[patchI], faceI)
            {
                if (lambda_.boundaryField()[patchI][faceI] < (1.0 - SMALL))
                {
                    activeYield_.boundaryField()[patchI][faceI] = 1.0;
                }
                else
                {
                    activeYield_.boundaryField()[patchI][faceI] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;

    // if (mesh().time().outputTime())
    // {
    //     Info<< "Writing det(bEbar)" << endl;

    //     volScalarField detBEbarMinus1
    //     (
    //         "detBEbarMinus1",
    //         det(bEbar_) - 1.0
    //     );

    //     detBEbarMinus1.write();

    //     volScalarField trBEbarOver3Minus1
    //     (
    //         "trBEbarOver3Minus1",
    //         (tr(bEbar_)/3.0) - 1.0
    //     );

    //     trBEbarOver3Minus1.write();
    // }
}


Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin::newDeltaT()
{
    notImplemented(type() + "::newDeltaT()");

    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformation gradient
    // if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    // {
    //     F() = fvc::average(relFf()) & F().oldTime();
    // }
    //else
    // {
    //     F() = relF() & F().oldTime();
    // }

    // // Calculate the total true (Hencky) strain
    // const volSymmTensorField epsilon = 0.5*log(symm(F().T() & F()));

    // // Calculate equivalent strain, for normalisation of the error
    // const volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));

    // // Take reference to internal fields
    // const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    // const symmTensorField& plasticNI = plasticN_.internalField();
    // const symmTensorField& plasticNIold =
    //    plasticN_.oldTime().internalField();
    // const scalarField& epsilonEqI = epsilonEq.internalField();

    // // Calculate error field
    // const symmTensorField DEpsilonPErrorI =
    //     Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
    //    /(epsilonEqI + SMALL);

    // // Max error
    // const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    // if (maxMagDEpsilonPErr > SMALL)
    // {
    //     Info<< "    " << name() << ": max time integration error = "
    //         << maxMagDEpsilonPErr
    //         << endl;

    //     if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
    //     {
    //         WarningIn
    //         (
    //             "Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin"
    //             "::newDeltaT()"
    //             " const"
    //        )   << "The error in the plastic strain is lover 50 times larger "
    //             << "than the desired value!\n    Consider starting the "
    //             << "simulation with a smaller initial time-step" << endl;
    //     }

    //     // Calculate the time-step scaling factor, where maxDeltaErr_ is the
    //     // maximum allowed error
    //     const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

    //     // Return the new time-step size
    //     return scaleFac*mesh().time().deltaTValue();
    // }

    return mesh().time().endTime().value();
}


// ************************************************************************* //
