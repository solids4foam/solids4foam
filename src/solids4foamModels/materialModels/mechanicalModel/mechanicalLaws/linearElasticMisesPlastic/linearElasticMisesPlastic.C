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

#include "linearElasticMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticMisesPlastic, linGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar linearElasticMisesPlastic::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label linearElasticMisesPlastic::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar linearElasticMisesPlastic::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar linearElasticMisesPlastic::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::scalar Foam::linearElasticMisesPlastic::curYieldStress
(
    const scalar curEpsilonPEq    // Current equivalent plastic strain
) const
{
    return stressPlasticStrainSeries_(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::linearElasticMisesPlastic::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar            // Scaled shear modulus
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
           *curYieldStress
            (
                epsilonPEqOld + sqrtTwoOverThree_*DLambda
            );
}


void Foam::linearElasticMisesPlastic::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar);

        if (i == MaxNewtonIter_)
        {
            WarningIn("linearElasticMisesPlastic::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Note: we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticMisesPlastic::linearElasticMisesPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    stressPlasticStrainSeries_(dict),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    sigmaHydf_
    (
        IOobject
        (
            "sigmaHydf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
    sigmaYf_
    (
        IOobject
        (
            "sigmaYf",
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
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    DSigmaYf_
    (
        IOobject
        (
            "DSigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonf_
    (
        IOobject
        (
            "epsilonf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    epsilonPf_
    (
        IOobject
        (
            "epsilonPf",
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
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    DEpsilonPEqf_
    (
        IOobject
        (
            "DEpsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
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
    DLambdaf_
    (
        IOobject
        (
            "DLambdaf",
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
    epsilonPEqf_
    (
        IOobject
        (
            "epsilonPEqf",
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
    plasticNf_
    (
        IOobject
        (
            "plasticNf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    epsilon_.oldTime();
    plasticN_.oldTime();

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
        K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;

        // This should no longer be performed because we directly update
        // epsilon_ZZ
        // if (planeStress())
        // {
        //     K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        // }
        // else
        // {
        //     K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        // }
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
            "linearElasticMisesPlastic::linearElasticMisesPlastic::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Check if plasticity is a nonlinear function of plastic strain
    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
    }
    else
    {
        if (stressPlasticStrainSeries_.size() == 1)
        {
            Info<< "    Perfect Plasticity" << endl;
        }
        else
        {
            Info<< "    Plasticity is linear" << endl;

            // Define linear plastic modulus
            Hp_ =
                (
                    stressPlasticStrainSeries_[1].second()
                  - stressPlasticStrainSeries_[0].second()
                )
               /(
                    stressPlasticStrainSeries_[1].first()
                  - stressPlasticStrainSeries_[0].first()
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticMisesPlastic::~linearElasticMisesPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlastic::rho() const
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
Foam::linearElasticMisesPlastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField scaleFactor = 1.0 - (2.0*mu_*DLambda_/magSTrial);

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
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


Foam::tmp<Foam::volDiagTensorField>
Foam::linearElasticMisesPlastic::impKdiagTensor() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField theta = 1.0 - (2.0*mu_*DLambda_/magSTrial);

    // Calculate N squared where N is the plastic return direction
    const volTensorField NsquaredTensor = plasticN_ & plasticN_;
    volDiagTensorField Nsquared
    (
        IOobject
        (
            "Nsquared",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedDiagTensor("zero", dimless, diagTensor::zero)
    );

    Nsquared.internalField() = diag(NsquaredTensor.internalField());

    forAll(Nsquared.boundaryField(), patchI)
    {
        Nsquared.boundaryField()[patchI] =
            diag(NsquaredTensor.boundaryField()[patchI]);
    }

    const diagTensor Idiag = diagTensor(1.0, 1.0, 1.0);

    return tmp<volDiagTensorField>
    (
        new volDiagTensorField
        (
            IOobject
            (
                "impKdiagTensor",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            K_*Idiag + mu_*0.5*theta*(Idiag*4.0/3.0 - 2.0*Nsquared)
            //K_*Idiag + mu_*theta*(Idiag*4.0/3.0)
            //K_*Idiag + mu_*(theta/theta)*(Idiag*4.0/3.0)
        )
    );
}


void Foam::linearElasticMisesPlastic::correct(volSymmTensorField& sigma)
{
    // Calculate total strain
    if (mesh().foundObject<volTensorField>("grad(DD)"))
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon_ = epsilon_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon_ = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        // Poisson's ratio
        const dimensionedScalar nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));

        // Young's modulus
        const dimensionedScalar E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        epsilon_.replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
          - (
                epsilonP_.component(symmTensor::XX)
              + epsilonP_.component(symmTensor::YY)
            )
        );
    }

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Calculate the yield function
    const volScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*sigmaY_;

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon_.internalField())), SMALL);

    // Store previous iteration for under-relaxation
    DEpsilonP_.storePrevIter();

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& sigmaYI = sigmaY_.internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();

    // Calculate DLambda_ and plasticN_
    forAll(fTrialI, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);
        if (magS > SMALL)
        {
            plasticNI[cellI] = sTrialI[cellI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
            plasticNI[cellI] = symmTensor::zero;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain where t is start of time-step
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaI[cellI],
                    curSigmaY,
                    epsilonPEqOldI[cellI],
                    magS,
                    mu_.value(),
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[cellI] = fTrialI[cellI]/(2*mu_.value());

                if (magHp > SMALL)
                {
                    DLambdaI[cellI] /= 1.0 + Hp_/(3*mu_.value());

                    // Update increment of yield stress
                    DSigmaYI[cellI] = DLambdaI[cellI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = sTrialP[faceI]/magS;
            }

            // Calculate DEpsilonPEq
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
                plasticNP[faceI] = symmTensor::zero;
            }
            else
            {
                // yielding
                if (nonLinearPlasticity_)
                {
                    scalar curSigmaY = 0.0; // updated in loop below

                    // Calculate DEpsilonPEq and curSigmaY
                    newtonLoop
                        (
                            DLambdaP[faceI],
                            curSigmaY,
                            epsilonPEqOldP[faceI],
                            magS,
                            mu_.value(),
                            maxMagBE
                        );

                    // Update increment of yield stress
                    DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                }
                else
                {
                    // Plastic modulus is linear
                    DLambdaP[faceI] = fTrialP[faceI]/(2.0*mu_.value());

                    if (magHp > SMALL)
                    {
                        DLambdaP[faceI] /= 1.0 + Hp_/(3.0*mu_.value());

                        // Update increment of yield stress
                        DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                    }
                }
            }
        }
    }

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;
    DEpsilonP_ = DLambda_*plasticN_;
    DEpsilonP_.relax();

    // Calculate deviatoric stress
    const volSymmTensorField s = sTrial - 2*mu_*DEpsilonP_;

    // Calculate the hydrostatic pressure directly from the displacement
    // field
    sigmaHyd_ = K_*tr(epsilon_);

    // Update the stress
    sigma = sigmaHyd_*I + s;
}


void Foam::linearElasticMisesPlastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate total strain
    if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        epsilonf_ = epsilonf_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        epsilonf_ = symm(gradD);
    }

    // Calculate deviatoric strain
    const surfaceSymmTensorField e = dev(epsilonf_);

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial =
        2.0*mu_*(e - dev(epsilonPf_.oldTime()));

    // Calculate the yield function
    const surfaceScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*sigmaYf_;

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilonf_.internalField())), SMALL);

    // Store previous iteration for under-relaxation
    DEpsilonPf_.storePrevIter();

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticNf_.internalField();
    scalarField& DSigmaYI = DSigmaYf_.internalField();
    scalarField& DLambdaI = DLambdaf_.internalField();
    const scalarField& sigmaYI = sigmaYf_.internalField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().internalField();

    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrialI, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);
        if (magS > SMALL)
        {
            plasticNI[cellI] = sTrialI[cellI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
            plasticNI[cellI] = symmTensor::zero;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain where t is start of time-step
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaI[cellI],
                    curSigmaY,
                    epsilonPEqOldI[cellI],
                    magS,
                    mu_.value(),
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[cellI] = fTrialI[cellI]/(2*mu_.value());

                if (magHp > SMALL)
                {
                    DLambdaI[cellI] /= 1.0 + Hp_/(3*mu_.value());

                    // Update increment of yield stress
                    DSigmaYI[cellI] = DLambdaI[cellI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticNf_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaYf_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambdaf_.boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaYf_.boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEqf_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = sTrialP[faceI]/magS;
            }

            // Calculate DEpsilonPEq
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
                plasticNP[faceI] = symmTensor::zero;
            }
            else
            {
                // yielding
                if (nonLinearPlasticity_)
                {
                    scalar curSigmaY = 0.0; // updated in loop below

                    // Calculate DEpsilonPEq and curSigmaY
                    newtonLoop
                    (
                        DLambdaP[faceI],
                        curSigmaY,
                        epsilonPEqOldP[faceI],
                        magS,
                        mu_.value(),
                        maxMagBE
                    );

                    // Update increment of yield stress
                    DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                }
                else
                {
                    // Plastic modulus is linear
                    DLambdaP[faceI] = fTrialP[faceI]/(2.0*mu_.value());

                    if (magHp > SMALL)
                    {
                        DLambdaP[faceI] /= 1.0 + Hp_/(3.0*mu_.value());

                        // Update increment of yield stress
                        DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                    }
                }
            }
        }
    }

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;
    DEpsilonPf_ = DLambdaf_*plasticNf_;
    DEpsilonPf_.relax();

    // Calculate deviatoric stress
    const surfaceSymmTensorField s = sTrial - 2*mu_*DEpsilonPf_;

    // Calculate the hydrostatic pressure directly from the displacement
    // field
    sigmaHydf_ = K_*tr(epsilonf_);

    // Update the stress
    sigma = sigmaHydf_*I + s;
}


Foam::scalar Foam::linearElasticMisesPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().foundObject<surfaceTensorField>("grad(D)f")
     || mesh().foundObject<surfaceTensorField>("grad(DD)f")
    )
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonPf_.internalField()
                  - DEpsilonPf_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().internalField()));
    }
    else
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonP_.internalField()
                  - DEpsilonP_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
    }
}


void Foam::linearElasticMisesPlastic::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;
    sigmaYf_ += DSigmaYf_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonPEqf_ += DEpsilonPEqf_;
    epsilonP_ += DEpsilonP_;
    epsilonPf_ += DEpsilonPf_;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.internalField(), celli)
    {
        if (DEpsilonPEq_.internalField()[celli] > SMALL)
        {
            activeYield_.internalField()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.internalField()[celli] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchi)
    {
        if (!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
                {
                    activeYield_.boundaryField()[patchi][facei] = 1.0;
                }
                else
                {
                    activeYield_.boundaryField()[patchi][facei] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    const int nTotalCells = returnReduce(mesh().nCells(), sumOp<int>());

    Info<< "    " << numCellsYielding << " cells ("
        << 100.0*scalar(numCellsYielding)/scalar(nTotalCells)
        << "% of the cells in this material) are actively yielding"
        << nl << endl;
}


Foam::scalar Foam::linearElasticMisesPlastic::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the Euler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon_)));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();

    // Calculate error field
    const symmTensorField DEpsilonPErrorI =
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL);

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max time integration error = "
            << maxMagDEpsilonPErr
            << endl;

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningIn
            (
                "Foam::scalar Foam::linearElasticMisesPlastic::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is over 50 times larger "
                << "than the desired value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


// ************************************************************************* //
