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

#include "neoHookeanElasticMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElasticMisesPlastic, dictionary
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

        // Tolerance for Newton loop
        scalar neoHookeanElasticMisesPlastic::LoopTol_ = 1e-8;

        // Maximum number of iterations for Newton loop
        label neoHookeanElasticMisesPlastic::MaxNewtonIter_ = 200;

        // finiteDiff is the delta for finite difference differentiation
        scalar neoHookeanElasticMisesPlastic::finiteDiff_ = 0.25e-6;

        // Store sqrt(2/3) as we use it often
        scalar neoHookeanElasticMisesPlastic::sqrtTwoOverThree_ =
            ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::neoHookeanElasticMisesPlastic::calcCurMaterialf() const
{
    if (curMaterialfPtr_)
    {
        FatalErrorIn
        (
            "void Foam::neoHookeanElasticMisesPlastic::calcCurMaterialf() const"
        )   << "pointer already set" << abort(FatalError);
    }

    curMaterialfPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "curMaterialf",
                mesh().time().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        );

    const volScalarField& curMat = curMaterial();
    const scalarField& curMatI = curMat.internalField();
    scalarField& curMatfI = curMaterialfPtr_->internalField();
    const unallocLabelList& own = mesh().owner();

    forAll(curMatfI, faceI)
    {
        // Set to owner: OK except for interfaces
        curMatfI[faceI] = curMatI[own[faceI]];
    }

    forAll(curMat.boundaryField(), patchI)
    {
        curMaterialfPtr_->boundaryField()[patchI] =
            curMat.boundaryField()[patchI];
    }

    curMaterialfPtr_->correctBoundaryConditions();
}


const Foam::surfaceScalarField&
Foam::neoHookeanElasticMisesPlastic::curMaterialf() const
{
    if (!curMaterialfPtr_)
    {
        calcCurMaterialf();
    }

    return *curMaterialfPtr_;
}


Foam::scalar Foam::neoHookeanElasticMisesPlastic::curYieldStress
(
    const scalar curEpsilonPEq,    // Current equivalent plastic strain
    const scalar J                 // Current Jacobian
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*stressPlasticStrainSeries_(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::neoHookeanElasticMisesPlastic::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar,            // Scaled shear modulus
    const scalar J                 // Current Jacobian
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
                epsilonPEqOld + sqrtTwoOverThree_*DLambda,
                J
            );
}


void Foam::neoHookeanElasticMisesPlastic::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar J,                // Current Jacobian
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newtion's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar, J);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluaitons are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar, J
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar,  J);

        if (i == MaxNewtonIter_)
        {
            WarningIn("neoHookeanElasticMisesPlastic::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Bug-fix PC ZT 7-Nov-14
    // Note we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda,  J
        )/J;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElasticMisesPlastic::neoHookeanElasticMisesPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mechanicalLaw(name, mesh, dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    stressPlasticStrainSeries_(dict),
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
    bEbarTrialf_
    (
        IOobject
        (
            "bEbarTrialf",
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
    bEbarf_
    (
        IOobject
        (
            "bEbarf",
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
    curMaterialfPtr_(NULL),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
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

Foam::neoHookeanElasticMisesPlastic::~neoHookeanElasticMisesPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElasticMisesPlastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
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
Foam::neoHookeanElasticMisesPlastic::impK() const
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


void Foam::neoHookeanElasticMisesPlastic::correct(volSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Compute elastic predictor

    // Lookup relative deformation gradient
    const volTensorField& relF = mesh.lookupObject<volTensorField>("relF");

    // Lookup relative Jacobian
    const volScalarField& relJ = mesh.lookupObject<volScalarField>("relJ");

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    volTensorField relFbar = pow(relJ, -1.0/3.0)*relF;

    // Update bE trial
    bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Lookup the Jacobian of the deformation gradient
    const volScalarField& J = mesh.lookupObject<volScalarField>("J");

    // Calculate trial deviatoric stress
    volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

    const volScalarField Ibar = tr(bEbarTrial_)/3.0;
    const volScalarField muBar = Ibar*mu_;

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrial_.internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const volScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*J*sigmaY_;

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& muBarI = muBar.internalField();
    const scalarField& JI = J.internalField();
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
                    muBarI[cellI],
                    JI[cellI],
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[cellI] = fTrialI[cellI]/(2*muBarI[cellI]);

                if (magHp > SMALL)
                {
                    DLambdaI[cellI] /= 1.0 + Hp_/(3*muBarI[cellI]);

                    // Update increment of yield stress
                    DSigmaYI[cellI] = DLambdaI[cellI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        if (!fTrial.boundaryField()[patchI].coupled())
        {
            // Take references to the boundary patch fields for efficiency
            const scalarField& fTrialP = fTrial.boundaryField()[patchI];
            const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
            symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
            scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
            scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
            const scalarField& muBarP = muBar.boundaryField()[patchI];
            const scalarField& JP = J.boundaryField()[patchI];
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
                            muBarP[faceI],
                            JP[faceI],
                            maxMagBE
                        );

                        // Update increment of yield stress
                        DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                    }
                    else
                    {
                        // Plastic modulus is linear
                        DLambdaP[faceI] = fTrialP[faceI]/(2.0*muBarP[faceI]);

                        if (magHp > SMALL)
                        {
                            DLambdaP[faceI] /= 1.0 + Hp_/(3.0*muBarP[faceI]);

                            // Update increment of yield stress
                            DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                        }
                    }
                }
            }
        }
    }

    DSigmaY_.correctBoundaryConditions();
    DLambda_.correctBoundaryConditions();
    plasticN_.correctBoundaryConditions();

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;
    DEpsilonP_ = Ibar*DLambda_*plasticN_;
    DEpsilonP_.relax();

    // Calculate deviatoric stress
    const volSymmTensorField s = sTrial - 2*mu_*DEpsilonP_;

    // Calculate new Kirchhoff stress
    tau = 0.5*K_*(pow(J, 2) - 1)*I + s;

    // Update bEbar
    bEbar_ = (s/mu_) + Ibar*I;


    // Multi-material:
    //volSymmTensorField newTau("newTau", 0.5*K_*(pow(J, 2) - 1)*I + s);

    // Assign Kirchhoff stress
    // For now, to deal with multi-materials, we will multiply by curMaterial
    // index field so only cells in the current material are calculated:
    // we should be able to do this in a better way

    //tau = curMaterial()*newTau + (1.0 - curMaterial())*tau;

    // Update bEbar
    //const volSymmTensorField newBEbar = (s/mu_) + Ibar*I;
    //bEbar_ = curMaterial()*newBEbar + (1.0 - curMaterial())*bEbar_;
}


void Foam::neoHookeanElasticMisesPlastic::correct(surfaceSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Compute elastic predictor

    // Lookup relative deformation gradient
    const surfaceTensorField& relF =
        mesh.lookupObject<surfaceTensorField>("relFf");

    // Lookup relative Jacobian
    const surfaceScalarField& relJ =
        mesh.lookupObject<surfaceScalarField>("relJf");

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    surfaceTensorField relFbar = pow(relJ, -1.0/3.0)*relF;

    // Update bE trial
    bEbarTrialf_ = transform(relFbar, bEbarf_.oldTime());

    // Lookup the Jacobian of the deformation gradient
    const surfaceScalarField& J = mesh.lookupObject<surfaceScalarField>("Jf");

    // Calculate trial deviatoric stress
    surfaceSymmTensorField sTrial = mu_*dev(bEbarTrialf_);

    const surfaceScalarField Ibar = tr(bEbarTrialf_)/3.0;
    const surfaceScalarField muBar = Ibar*mu_;

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonPf_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrialf_.internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const surfaceScalarField fTrial =
        mag(sTrial) - sqrtTwoOverThree_*J*sigmaYf_;

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticNf_.internalField();
    scalarField& DSigmaYI = DSigmaYf_.internalField();
    scalarField& DLambdaI = DLambdaf_.internalField();
    const scalarField& muBarI = muBar.internalField();
    const scalarField& JI = J.internalField();
    const scalarField& sigmaYI = sigmaYf_.internalField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().internalField();


    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrialI, faceI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[faceI]);
        if (magS > SMALL)
        {
            plasticNI[faceI] = sTrialI[faceI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[faceI] < SMALL)
        {
            // elastic
            DSigmaYI[faceI] = 0.0;
            DLambdaI[faceI] = 0.0;
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
                    DLambdaI[faceI],
                    curSigmaY,
                    epsilonPEqOldI[faceI],
                    magS,
                    muBarI[faceI],
                    JI[faceI],
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[faceI] = curSigmaY - sigmaYI[faceI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[faceI] = fTrialI[faceI]/(2*muBarI[faceI]);

                if (magHp > SMALL)
                {
                    DLambdaI[faceI] /= 1.0 + Hp_/(3*muBarI[faceI]);

                    // Update increment of yield stress
                    DSigmaYI[faceI] = DLambdaI[faceI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Calculate on coupled for surface fields
        //if (!fTrial.boundaryField()[patchI].coupled())
        {
            // Take references to the boundary patch fields for efficiency
            const scalarField& fTrialP = fTrial.boundaryField()[patchI];
            const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
            symmTensorField& plasticNP = plasticNf_.boundaryField()[patchI];
            scalarField& DSigmaYP = DSigmaYf_.boundaryField()[patchI];
            scalarField& DLambdaP = DLambdaf_.boundaryField()[patchI];
            const scalarField& muBarP = muBar.boundaryField()[patchI];
            const scalarField& JP = J.boundaryField()[patchI];
            const scalarField& sigmaYP = sigmaYf_.boundaryField()[patchI];
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
                            muBarP[faceI],
                            JP[faceI],
                            maxMagBE
                        );

                        // Update increment of yield stress
                        DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                    }
                    else
                    {
                        // Plastic modulus is linear
                        DLambdaP[faceI] = fTrialP[faceI]/(2.0*muBarP[faceI]);

                        if (magHp > SMALL)
                        {
                            DLambdaP[faceI] /= 1.0 + Hp_/(3.0*muBarP[faceI]);

                            // Update increment of yield stress
                            DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                        }
                    }
                }
            }
        }
    }

    DSigmaYf_.correctBoundaryConditions();
    DLambdaf_.correctBoundaryConditions();
    plasticNf_.correctBoundaryConditions();

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;
    DEpsilonPf_ = Ibar*DLambdaf_*plasticNf_;
    DEpsilonPf_.relax();

    // Calculate deviatoric stress
    const surfaceSymmTensorField s = sTrial - 2*mu_*DEpsilonPf_;

    // Calculate new Kirchhoff stress
    tau = 0.5*K_*(pow(J, 2) - 1)*I + s;

    // Update bEbar
    bEbarf_ = (s/mu_) + Ibar*I;


    // Multi-material:
    // surfaceSymmTensorField newTau("newTau", 0.5*K_*(pow(J, 2) - 1)*I + s);

    // Assign Kirchhoff stress
    // For now, to deal with multi-materials, we will multiply by curMaterial
    // index field so only cells in the current material are calculated:
    // we should be able to do this in a better way

    // tau = curMaterialf()*newTau + (1.0 - curMaterialf())*tau;

    // Update bEbar
    //const surfaceSymmTensorField newBEbar = (s/mu_) + Ibar*I;
    //bEbarf_ = curMaterialf()*newBEbar + (1.0 - curMaterialf())*bEbarf_;
}


Foam::scalar Foam::neoHookeanElasticMisesPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if (mesh().foundObject<surfaceTensorField>("Ff"))
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


void Foam::neoHookeanElasticMisesPlastic::updateTotalFields()
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

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


// Foam::tmp<Foam::volScalarField>
// Foam::neoHookeanElasticMisesPlastic::plasticDissipationRate() const
// {
//     const volSymmTensorField& sigmaCauchy =
//         mesh().lookupObject<volSymmTensorField>("sigmaCauchy");

//     // We assume 90% conversion
//     return tmp<volScalarField>
//     (
//         new volScalarField
//         (
//             "plasticDissipationRate",
//             max
//             (
//                 dimensionedScalar("zero", dimForce/(dimArea*dimTime), 0.0),
//                 0.9*(sigmaCauchy && DEpsilonP_)/mesh().time().deltaT()
//             )
//         )
//     );
// }


void Foam::neoHookeanElasticMisesPlastic::setMaterialIndex(label curMatIndex)
{
    // Set current material index
    curMaterialIndex() = curMatIndex;
}


// ************************************************************************* //
