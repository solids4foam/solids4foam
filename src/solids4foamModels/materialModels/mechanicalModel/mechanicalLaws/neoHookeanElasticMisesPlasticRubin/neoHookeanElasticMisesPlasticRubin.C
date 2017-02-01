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

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Store sqrt(2/3) as we use it often
    scalar neoHookeanElasticMisesPlasticRubin::sqrtTwoOverThree_ =
        ::sqrt(2.0/3.0);

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


void Foam::neoHookeanElasticMisesPlasticRubin::makeRelFf()
{
    if (relFfPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticMisesPlasticRubin::makeRelFf()")
            << "pointer already set" << abort(FatalError);
    }

    relFfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "lawRelFf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::surfaceTensorField& Foam::neoHookeanElasticMisesPlasticRubin::relFf()
{
    if (!relFfPtr_)
    {
        makeRelFf();
    }

    return *relFfPtr_;
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


void Foam::neoHookeanElasticMisesPlasticRubin::makeJf()
{
    if (JfPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticMisesPlasticRubin::makeJf()")
            << "pointer already set" << abort(FatalError);
    }

    JfPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "lawJf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );

    // Store the old-time
    JfPtr_->oldTime();
}


Foam::surfaceScalarField& Foam::neoHookeanElasticMisesPlasticRubin::Jf()
{
    if (!JfPtr_)
    {
        makeJf();
    }

    return *JfPtr_;
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


void Foam::neoHookeanElasticMisesPlasticRubin::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanElasticMisesPlasticRubin::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "lawFf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );

    // Store the old-time
    FfPtr_->oldTime();
}


Foam::surfaceTensorField& Foam::neoHookeanElasticMisesPlasticRubin::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElasticMisesPlasticRubin::Ibar
(
    const volSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
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

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
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
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
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
    relFfPtr_(NULL),
    JPtr_(NULL),
    JfPtr_(NULL),
    FPtr_(NULL),
    FfPtr_(NULL),
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
        dimensionedScalar("1", dimless, 1.0)
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
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
        dimensionedScalar("one", dimless, 1.0)
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
    ),
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    plasticCompaction_(dict.lookup("plasticCompaction")),
    rubin_(dict.lookupOrDefault<Switch>("rubin", true))
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
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
    deleteDemandDrivenData(relFfPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(JfPtr_);
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::
neoHookeanElasticMisesPlasticRubin::rho() const
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

    if (rubin_)
    {
        // Calculate sigmaStar
        // Note: tr(bEbar & bEbar) == magSqr(bEbar)

        const volSymmTensorField devBEbarTrial = dev(bEbarTrial_);

        const volScalarField sigmaStar =
            mu_*sqrt((3.0/2.0)*tr(devBEbarTrial & devBEbarTrial));

        // Calculate yield function
        const volScalarField fTrial = sigmaStar/sigmaY_;

        // Or we can enforce Cauchy yield stress instead of Kirchhoff yield
        // stress
        //const volScalarField fTrial = sigmaStar/(J()*sigmaY_);

        // Store previous iteration for calculation of the residual
        lambda_.storePrevIter();

        // Calculate plastic multiplier
        lambda_ = 1.0/max(1.0, fTrial);

        // Deviatoric component of bEbar
        const volSymmTensorField devBEbar = lambda_*dev(bEbarTrial_);

        // Calculate bEBar
        bEbar_ = devBEbar + Ibar(devBEbar)*I;

        // Calcualte Je
        if (plasticCompaction_)
        {
            Je_ =
                relJ*Je_.oldTime()
               *exp
                (
                    (3.0/2.0)*((1.0 - lambda_)/lambda_)
                    *(
                        1.0 - (9.0/(tr(bEbar_)*tr(inv(bEbar_))))
                    )
                );
        }
        else
        {
            Je_ = 1.0*J();
        }

        // Reciprocal of J
        const volScalarField rJ = 1.0/J();

        // Calculate the deviatoric Cauchy stress
        const volSymmTensorField Tprime = rJ*mu_*dev(bEbar_);

        // Update hydrostatic pressure
        // Note: P == -sigmaHyd
        sigmaHyd_ = -rJ*0.5*K_*(1.0 - pow(Je_, 2.0));

        // Update the Cauchy stress
        sigma = Tprime + sigmaHyd_*I;
    }
    else
    {
        // Calculate trial deviatoric stress
        const volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

        const volScalarField Ibar = tr(bEbarTrial_)/3.0;
        const volScalarField muBar = Ibar*mu_;

        // Check for plastic loading
        // and calculate increment of plastic equivalent strain
        // i.e. the plastic multiplier

        // Trial yield function
        // sigmaY is the Cauchy yield stress so we scale it by J
        const volScalarField fTrial =
            mag(sTrial) - sqrtTwoOverThree_*J()*sigmaY_;

        // Return direction
        plasticN_ =
            sTrial
           /(mag(sTrial) + dimensionedScalar("SMALL", dimPressure, SMALL));

        // Store previous iteration for calculation of the residual
        DLambda_.storePrevIter();

        // Plastic multiplier
        DLambda_ = pos(fTrial)*(fTrial/(2.0*muBar));

        // Calculate deviatoric stress
        const volSymmTensorField s = sTrial - 2.0*mu_*Ibar*DLambda_*plasticN_;

        // Calculate the deviatoric component of bEbar
        const volSymmTensorField devBEbar = s/mu_;

        // Update bEbar
        bEbar_ = devBEbar + this->Ibar(devBEbar)*I;

        // Update hydrostatic stress (negative of pressure)
        sigmaHyd_ = 0.5*K_*(pow(J(), 2.0) - 1.0);

        // Update the Cauchy stress
        sigma = (1.0/J())*(sigmaHyd_*I + s);
    }
}


void Foam::neoHookeanElasticMisesPlasticRubin::correct
(
    surfaceSymmTensorField& sigma
)
{
    notImplemented(type() + " for surface fields");
}


Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin::residual()
{
    // Note: we remove mag(I) when we normalise the residual so that the
    // residual is like a strain residual
    if (rubin_)
    {
        return
        gMax
        (
            mag
            (
                lambda_.internalField()
              - lambda_.prevIter().internalField()
            )
        )/0.1;
    }
    else
    {
        return
        gMax
        (
            mag
            (
                DLambda_.internalField()
              - DLambda_.prevIter().internalField()
            )
        )/(gMax(mag(DLambda_.prevIter().internalField())) + SMALL);
    }
}


void Foam::neoHookeanElasticMisesPlasticRubin::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;

    // Count cells actively yielding
    int numCellsYielding = 0;

    scalarField& activeYieldI = activeYield_.internalField();
    const scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& lambdaI = lambda_.internalField();

    forAll(activeYieldI, cellI)
    {
        if (DLambdaI[cellI] > SMALL || lambdaI[cellI] < (1.0 - SMALL))
        {
            activeYieldI[cellI] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYieldI[cellI] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchI)
    {
        if (!activeYield_.boundaryField()[patchI].coupled())
        {
            scalarField& activeYieldP = activeYield_.boundaryField()[patchI];
            const scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
            const scalarField& lambdaP = lambda_.boundaryField()[patchI];

            forAll(activeYieldP, faceI)
            {
                if (DLambdaP[faceI] > SMALL || lambdaP[faceI] < (1.0 - SMALL))
                {
                    activeYieldP[faceI] = 1.0;
                }
                else
                {
                    activeYieldP[faceI] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;

    if (mesh().time().outputTime())
    {
        Info<< "Writing tauEq" << endl;

        // Calculate the deviatoric stress
        const volSymmTensorField s = mu_*dev(bEbar_);

        // Calculate the Kirchhoff stress
        const volSymmTensorField tau("tau", s + sigmaHyd_*I);

        // Calculate the equivalent Kirchhoff stress
        const volScalarField tauEq("tauEq", sqrt((3.0/2.0)*magSqr(dev(tau))));

        // Write the field
        tauEq.write();

        // Plastic component of total Jacobian
        Info<< "Writing Jp" << endl;

        const volScalarField Jp("Jp", J()/Je_);

        // Write the field
        Jp.write();
    }
}


Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin::newDeltaT()
{
    return mesh().time().endTime().value();
}


// ************************************************************************* //
