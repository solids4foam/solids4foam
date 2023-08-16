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


void Foam::linearElasticMisesPlastic::updatePlasticity
(
    symmTensor& plasticN,          // Plastic return direction
    scalar& DLambda,               // Plastic multiplier increment
    scalar& DSigmaY,               // Increment of yield stress
    scalar& sigmaY,                // Yield stress
    const scalar sigmaYOld,        // Yield stress old time
    const scalar fTrial,           // Trial yield function
    const symmTensor& sTrial,      // Trial deviatoric stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagBE          // Max strain increment magnitude
) const
{
    // Calculate DLambda/DEpsilonPEq
    if (fTrial < SMALL)
    {
        // Elasticity
        plasticN = symmTensor(I);
        DLambda = 0.0;
        DSigmaY = 0.0;
        sigmaY = sigmaYOld;
    }
    else
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrial);
        if (magS > SMALL)
        {
            plasticN = sTrial/magS;
        }
        else
        {
            // Deviatoric stress is zero so plasticN value does not matter, but
            // we will set it to the identity
            plasticN = symmTensor(I);
        }

        if (nonLinearPlasticity_)
        {
            // Update plastic multiplier (DLambda) and current yield stress
            // (sigmaY)
            newtonLoop
            (
                DLambda,
                sigmaY,
                epsilonPEqOld,
                magS,
                mu_.value(),
                maxMagBE
            );

            // Update increment of yield stress
            DSigmaY = sigmaY - sigmaYOld;
        }
        else
        {
            // Update DLambda
            DLambda = fTrial/(2*mu_.value());

            // If the isotropic linear modulus is non-zero
            if (mag(Hp_) > SMALL)
            {
                DLambda /= 1.0 + Hp_/(3*mu_.value());

                // Update increment of yield stress
                DSigmaY = sqrtTwoOverThree_*DLambda*Hp_;

                // Update yield stress
                sigmaY = sigmaYOld + DSigmaY;
            }
        }
    }
}


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
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda
        );
}


void Foam::linearElasticMisesPlastic::calculateStress
(
 surfaceSymmTensorField& sigma,
 const surfaceSymmTensorField& epsilon
 ) const
{
    // Calculate deviatoric strain
    const surfaceSymmTensorField e(dev(epsilon));

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial
    (
     2.0*mu_*(e - dev(epsilonPf_.oldTime()))
     );

    // Calculate the yield function
    const surfaceScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*sigmaYf_);

    // Make a copy of history fields that are updated
    surfaceSymmTensorField plasticN("plasticNtmp", 1.0*plasticNf_);
    surfaceScalarField DSigmaY("DSigmaYtmp", 1.0*DSigmaYf_);
    surfaceScalarField DLambda("DLambdatmp", 1.0*DLambdaf_);
    surfaceScalarField sigmaY("sigmaYtmp", 1.0*sigmaYf_);

#ifdef OPENFOAM_NOT_EXTEND
    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon.primitiveField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticN.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaY.primitiveFieldRef();
    scalarField& DLambdaI = DLambda.primitiveFieldRef();
    scalarField& sigmaYI = sigmaY.primitiveFieldRef();
    const scalarField& sigmaYOldI = sigmaYf_.oldTime().primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().primitiveField();
#else
    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon.internalField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN.internalField();
    scalarField& DSigmaYI = DSigmaY.internalField();
    scalarField& DLambdaI = DLambda.internalField();
    scalarField& sigmaYI = sigmaY.internalField();
    const scalarField& sigmaYOldI = sigmaYf_.oldTime().internalField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().internalField();
#endif

    // Calculate DLambdaf_ and plasticNf_
    // int numYield = 0;
    forAll(fTrialI, faceI)
    {
        // Update plasticN, DLambda, DSigmaY and sigmaY for this face
        updatePlasticity
        (
         plasticNI[faceI],
         DLambdaI[faceI],
         DSigmaYI[faceI],
         sigmaYI[faceI],
         sigmaYOldI[faceI],
         fTrialI[faceI],
         sTrialI[faceI],
         epsilonPEqOldI[faceI],
         mu_.value(),
         maxMagBE
         );

        // if (fTrialI[faceI] > 0)
        // {
        //     numYield++;
        // }
    }
    // Info<< "        tang: numYield = " <<  numYield << endl;

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
#ifdef OPENFOAM_NOT_EXTEND
        symmTensorField& plasticNP = plasticN.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaY.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambda.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaY.boundaryFieldRef()[patchI];
#else
        symmTensorField& plasticNP = plasticN.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaY.boundaryField()[patchI];
#endif
        const scalarField& sigmaYOldP =
        sigmaYf_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
        epsilonPEqf_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            updatePlasticity
            (
             plasticNP[faceI],
             DLambdaP[faceI],
             DSigmaYP[faceI],
             sigmaYP[faceI],
             sigmaYOldP[faceI],
             fTrialP[faceI],
             sTrialP[faceI],
             epsilonPEqOldP[faceI],
             mu_.value(),
             maxMagBE
             );
        }
    }

    // Update DEpsilonPEq
    // DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;

    // Store previous iteration for residual calculation
    // DEpsilonPf_.storePrevIter();

    // Update DEpsilonP
    const surfaceSymmTensorField DEpsilonP(DLambda*plasticN);

    // Update total plastic strain
    // epsilonPf_ = epsilonPf_.oldTime() + DEpsilonPf_;

    // Update equivalent total plastic strain
    // epsilonPEqf_ = epsilonPEqf_.oldTime() + DEpsilonPEqf_;

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(sTrial - 2*mu_*DEpsilonP);

    // Calculate the hydrostatic pressure directly from the displacement
    // field
    // const surfaceScalarField trEpsilon(tr(epsilon));
    // calculateHydrostaticStress(sigmaHydf_, trEpsilon);

    // Update the stress
    sigma = K_*tr(epsilon)*I + s;
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
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
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
    if (planeStress())
    {
        FatalErrorIn
        (
            "linearElasticMisesPlastic::linearElasticMisesPlastic::()"
        )   << "Not implemented for planeStress. If needed, you can solve the "
            << "case in 3-D and set the back to a symmetryPlane and the front "
            << "to traction-free" << abort(FatalError);
    }

    // Force storage of old-time fields
    epsilon().oldTime();
    epsilonf().oldTime();
    epsilonP_.oldTime();
    epsilonPf_.oldTime();
    epsilonPEq_.oldTime();
    epsilonPEqf_.oldTime();
    plasticN_.oldTime();
    sigmaY_.oldTime();
    sigmaYf_.oldTime();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_ = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        nu_ = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        // Set the bulk modulus
        K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        // Read shear modulus
        mu_ = dimensionedScalar(dict.lookup("mu"));

        // Read bulk modulus
        K_ = dimensionedScalar(dict.lookup("K"));

        // Calculate Young's modulus
        E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        // Calculate Poisson's ratio
        nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
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

Foam::tmp<Foam::volScalarField>
Foam::linearElasticMisesPlastic::bulkModulus() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "bulkModulus",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            K_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    tresult.ref().correctBoundaryConditions();
#else
    tresult().correctBoundaryConditions();
#endif

    return tresult;
}


Foam::tmp<Foam::volScalarField>
Foam::linearElasticMisesPlastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(epsilon()));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*mu_*DLambda_/magSTrial));

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


#ifndef OPENFOAM_NOT_EXTEND
Foam::tmp<Foam::volDiagTensorField>
Foam::linearElasticMisesPlastic::impKdiagTensor() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(epsilon()));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField theta(1.0 - (2.0*mu_*DLambda_/magSTrial));

    // Calculate N squared where N is the plastic return direction
    const volTensorField NsquaredTensor(plasticN_ & plasticN_);
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
#endif


#ifdef OPENFOAM_NOT_EXTEND
Foam::tmp<Foam::Field<Foam::scalarSquareMatrix>>
Foam::linearElasticMisesPlastic::materialTangentField() const
{
    // Prepare tmp field
    tmp<Field<Foam::scalarSquareMatrix>> tresult
    (
        new Field<Foam::scalarSquareMatrix>
        (
            mesh().nFaces(), Foam::scalarSquareMatrix(6, 0.0)
        )
    );
    Field<Foam::scalarSquareMatrix>& result = tresult.ref();

    // Update total strain
    const_cast<linearElasticMisesPlastic&>(*this).updateEpsilonf();

    // Calculate deviatoric strain
    const surfaceSymmTensorField e(dev(epsilonf()));

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonPf_.oldTime())));

    // Magnitude of the deviatoric trial stress
    const surfaceScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate the yield function
    const surfaceScalarField fTrial
    (
        mag(sTrial) - sqrtTwoOverThree_*sigmaYf_.oldTime()
    );
    const scalarField& fTrialI = fTrial.internalField();

    // Return direction
    const surfaceSymmTensorField plasticN(sTrial/magSTrial);
    const symmTensorField& plasticNI = plasticN.internalField();

    // Calculate tangent field
    const Switch numericalTangent(dict().lookup("numericalTangent"));
    if (numericalTangent)
    {
        // Lookup current stress and store it as the reference
        const surfaceSymmTensorField& sigmaRef =
            mesh().lookupObject<surfaceSymmTensorField>("sigmaf");
        const surfaceSymmTensorField& epsilonRef = epsilonf();

        // Create fields to be used for perturbations
        surfaceSymmTensorField sigmaPerturb("sigmaPerturb", sigmaRef);
        surfaceSymmTensorField epsilonPerturb("epsilonPerturb", epsilonRef);

        // Small number used for perturbations
        const scalar eps(readScalar(dict().lookup("tangentEps")));

        // For each component of epsilon, sequentially apply a perturbation and
        // then calculate the resulting sigma
        for (label cmptI = 0; cmptI < symmTensor::nComponents; cmptI++)
        {
            // Reset epsilonPerturb
            // We multiply by 1.0 to avoid issues with epsilonf being removed
            // from the object registry
            epsilonPerturb = 1.0*epsilonRef;

            // Perturb this component of epsilon
            epsilonPerturb.replace(cmptI, epsilonPerturb.component(cmptI) + eps);

            // Calculate perturbed stress
            calculateStress(sigmaPerturb, epsilonPerturb);

            // Calculate tangent component
            const surfaceSymmTensorField tangCmpt((sigmaPerturb - sigmaRef)/eps);
            const symmTensorField& tangCmptI = tangCmpt.internalField();

            // Insert tangent component
            forAll(tangCmptI, faceI)
            {
                if (cmptI == symmTensor::XX)
                {
                    result[faceI](0,0) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,0) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,0) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,0) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,0) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,0) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == symmTensor::YY)
                {
                    result[faceI](0,1) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,1) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,1) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,1) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,1) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,1) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == symmTensor::ZZ)
                {
                    result[faceI](0,2) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,2) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,2) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,2) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,2) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,2) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == symmTensor::XY)
                {
                    result[faceI](0,3) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,3) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,3) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,3) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,3) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,3) = tangCmptI[faceI][symmTensor::XZ];
                }
                else if (cmptI == symmTensor::YZ)
                {
                    result[faceI](0,4) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,4) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,4) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,4) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,4) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,4) = tangCmptI[faceI][symmTensor::XZ];
                }
                else // if (cmptI == symmTensor::XZ)
                {
                    result[faceI](0,5) = tangCmptI[faceI][symmTensor::XX];
                    result[faceI](1,5) = tangCmptI[faceI][symmTensor::YY];
                    result[faceI](2,5) = tangCmptI[faceI][symmTensor::ZZ];
                    result[faceI](3,5) = tangCmptI[faceI][symmTensor::XY];
                    result[faceI](4,5) = tangCmptI[faceI][symmTensor::YZ];
                    result[faceI](5,5) = tangCmptI[faceI][symmTensor::XZ];
                }
            }

            forAll(tangCmpt.boundaryField(), patchI)
            {
                const symmTensorField& tangCmptP =
                    tangCmpt.boundaryField()[patchI];
                const label start = mesh().boundaryMesh()[patchI].start();

                forAll(tangCmptP, fI)
                {
                    const label faceID = start + fI;

                    if (cmptI == symmTensor::XX)
                    {
                        result[faceID](0,0) = tangCmptI[fI][symmTensor::XX];
                        result[faceID](1,0) = tangCmptI[fI][symmTensor::YY];
                        result[faceID](2,0) = tangCmptI[fI][symmTensor::ZZ];
                        result[faceID](3,0) = tangCmptI[fI][symmTensor::XY];
                        result[faceID](4,0) = tangCmptI[fI][symmTensor::YZ];
                        result[faceID](5,0) = tangCmptI[fI][symmTensor::XZ];
                    }
                    else if (cmptI == symmTensor::YY)
                    {
                        result[faceID](0,1) = tangCmptI[fI][symmTensor::XX];
                        result[faceID](1,1) = tangCmptI[fI][symmTensor::YY];
                        result[faceID](2,1) = tangCmptI[fI][symmTensor::ZZ];
                        result[faceID](3,1) = tangCmptI[fI][symmTensor::XY];
                        result[faceID](4,1) = tangCmptI[fI][symmTensor::YZ];
                        result[faceID](5,1) = tangCmptI[fI][symmTensor::XZ];
                    }
                    else if (cmptI == symmTensor::ZZ)
                    {
                        result[faceID](0,2) = tangCmptI[fI][symmTensor::XX];
                        result[faceID](1,2) = tangCmptI[fI][symmTensor::YY];
                        result[faceID](2,2) = tangCmptI[fI][symmTensor::ZZ];
                        result[faceID](3,2) = tangCmptI[fI][symmTensor::XY];
                        result[faceID](4,2) = tangCmptI[fI][symmTensor::YZ];
                        result[faceID](5,2) = tangCmptI[fI][symmTensor::XZ];
                    }
                    else if (cmptI == symmTensor::XY)
                    {
                        result[faceID](0,3) = tangCmptI[fI][symmTensor::XX];
                        result[faceID](1,3) = tangCmptI[fI][symmTensor::YY];
                        result[faceID](2,3) = tangCmptI[fI][symmTensor::ZZ];
                        result[faceID](3,3) = tangCmptI[fI][symmTensor::XY];
                        result[faceID](4,3) = tangCmptI[fI][symmTensor::YZ];
                        result[faceID](5,3) = tangCmptI[fI][symmTensor::XZ];
                    }
                    else if (cmptI == symmTensor::YZ)
                    {
                        result[faceID](0,4) = tangCmptI[fI][symmTensor::XX];
                        result[faceID](1,4) = tangCmptI[fI][symmTensor::YY];
                        result[faceID](2,4) = tangCmptI[fI][symmTensor::ZZ];
                        result[faceID](3,4) = tangCmptI[fI][symmTensor::XY];
                        result[faceID](4,4) = tangCmptI[fI][symmTensor::YZ];
                        result[faceID](5,4) = tangCmptI[fI][symmTensor::XZ];
                    }
                    else // if (cmptI == symmTensor::XZ)
                    {
                        result[faceID](0,5) = tangCmptI[fI][symmTensor::XX];
                        result[faceID](1,5) = tangCmptI[fI][symmTensor::YY];
                        result[faceID](2,5) = tangCmptI[fI][symmTensor::ZZ];
                        result[faceID](3,5) = tangCmptI[fI][symmTensor::XY];
                        result[faceID](4,5) = tangCmptI[fI][symmTensor::YZ];
                        result[faceID](5,5) = tangCmptI[fI][symmTensor::XZ];
                    }
                }
            }
        }

        // Include 0.5 factor for shear components
        forAll(result, faceI)
        {

            for (int i = 3; i < 6; i++)
            {
                for (int j = 3; j < 6; j++)
                {
                    result[faceI](i,j) *= 0.5;
                }
            }
        }
    }
    else // Analytical tangent
    {
        // Defined in box 3.2 of Simo and Hughes
        notImplemented("Analytical tangent not implemented");
    }

    return tresult;
}
#endif


void Foam::linearElasticMisesPlastic::correct(volSymmTensorField& sigma)
{
    // Calculate total strain
    updateEpsilon();

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(epsilon()));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Calculate the yield function
    const volScalarField fTrial
    (
        mag(sTrial) - sqrtTwoOverThree_*sigmaY_.oldTime()
    );

#ifdef OPENFOAM_NOT_EXTEND
    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon().primitiveField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticN_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaY_.primitiveFieldRef();
    scalarField& DLambdaI = DLambda_.primitiveFieldRef();
    scalarField& sigmaYI = sigmaY_.primitiveFieldRef();
    const scalarField& sigmaYOldI = sigmaY_.oldTime().primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().primitiveField();
#else
    const scalar maxMagBE = max(gMax(mag(epsilon().internalField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    scalarField& sigmaYI = sigmaY_.internalField();
    const scalarField& sigmaYOldI = sigmaY_.oldTime().internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();
#endif

    forAll(fTrialI, cellI)
    {
        // Update plasticN, DLambda, DSigmaY and sigmaY for this cell
        updatePlasticity
        (
            plasticNI[cellI],
            DLambdaI[cellI],
            DSigmaYI[cellI],
            sigmaYI[cellI],
            sigmaYOldI[cellI],
            fTrialI[cellI],
            sTrialI[cellI],
            epsilonPEqOldI[cellI],
            mu_.value(),
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
#ifdef OPENFOAM_NOT_EXTEND
        symmTensorField& plasticNP = plasticN_.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryFieldRef()[patchI];
#else
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
#endif
        const scalarField& sigmaYOldP =
            sigmaY_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            updatePlasticity
            (
                plasticNP[faceI],
                DLambdaP[faceI],
                DSigmaYP[faceI],
                sigmaYP[faceI],
                sigmaYOldP[faceI],
                fTrialP[faceI],
                sTrialP[faceI],
                epsilonPEqOldP[faceI],
                mu_.value(),
                maxMagBE
            );
        }
    }

    // Update DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;

    // Store previous iteration for residual calculation
    DEpsilonP_.storePrevIter();

    // Update DEpsilonP
    DEpsilonP_ = DLambda_*plasticN_;

    // Update total plastic strain
    epsilonP_ = epsilonP_.oldTime() + DEpsilonP_;

    // Update equivalent total plastic strain
    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq_;

    // Calculate deviatoric stress
    const volSymmTensorField s(sTrial - 2*mu_*DEpsilonP_);

    // Update hydrostatic stress
    updateSigmaHyd(K_*tr(epsilon()), (4.0/3.0)*mu_ + K_);

    // Update the stress
    sigma = sigmaHyd()*I + s;
}


void Foam::linearElasticMisesPlastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate total strain
    updateEpsilonf();

    // Calculate deviatoric strain
    const surfaceSymmTensorField e(dev(epsilonf()));

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial
    (
        2.0*mu_*(e - dev(epsilonPf_.oldTime()))
    );

    // Calculate the yield function
    const surfaceScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*sigmaYf_);

#ifdef OPENFOAM_NOT_EXTEND
    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilonf().primitiveField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticNf_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaYf_.primitiveFieldRef();
    scalarField& DLambdaI = DLambdaf_.primitiveFieldRef();
    scalarField& sigmaYI = sigmaYf_.primitiveFieldRef();
    const scalarField& sigmaYOldI = sigmaYf_.oldTime().primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().primitiveField();
#else
    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilonf().internalField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticNf_.internalField();
    scalarField& DSigmaYI = DSigmaYf_.internalField();
    scalarField& DLambdaI = DLambdaf_.internalField();
    scalarField& sigmaYI = sigmaYf_.internalField();
    const scalarField& sigmaYOldI = sigmaYf_.oldTime().internalField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().internalField();
#endif

    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrialI, faceI)
    {
        // Update plasticN, DLambda, DSigmaY and sigmaY for this face
        updatePlasticity
        (
            plasticNI[faceI],
            DLambdaI[faceI],
            DSigmaYI[faceI],
            sigmaYI[faceI],
            sigmaYOldI[faceI],
            fTrialI[faceI],
            sTrialI[faceI],
            epsilonPEqOldI[faceI],
            mu_.value(),
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
#ifdef OPENFOAM_NOT_EXTEND
        symmTensorField& plasticNP = plasticNf_.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaYf_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambdaf_.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaYf_.boundaryFieldRef()[patchI];
#else
        symmTensorField& plasticNP = plasticNf_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaYf_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambdaf_.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaYf_.boundaryField()[patchI];
#endif
        const scalarField& sigmaYOldP =
            sigmaYf_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEqf_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            updatePlasticity
            (
                plasticNP[faceI],
                DLambdaP[faceI],
                DSigmaYP[faceI],
                sigmaYP[faceI],
                sigmaYOldP[faceI],
                fTrialP[faceI],
                sTrialP[faceI],
                epsilonPEqOldP[faceI],
                mu_.value(),
                maxMagBE
            );
        }
    }

    // Update DEpsilonPEq
    DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;

    // Store previous iteration for residual calculation
    DEpsilonPf_.storePrevIter();

    // Update DEpsilonP
    DEpsilonPf_ = DLambdaf_*plasticNf_;

    // Update total plastic strain
    epsilonPf_ = epsilonPf_.oldTime() + DEpsilonPf_;

    // Update equivalent total plastic strain
    epsilonPEqf_ = epsilonPEqf_.oldTime() + DEpsilonPEqf_;

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(sTrial - 2*mu_*DEpsilonPf_);

    if (solvePressureEqn())
    {
        // Solve pressure equation at cells
        updateEpsilon();
        updateSigmaHyd(K_*tr(epsilon()), (4.0/3.0)*mu_ + K_);

        // Interpolate to faces
        const surfaceScalarField sigmaHydf(fvc::interpolate(sigmaHyd()));

        // Update the stress
        sigma = sigmaHydf*I + s;
    }
    else
    {
        // Calculate hydrostatic stress at the faces
        const surfaceScalarField sigmaHydf(K_*tr(epsilonf()));

        // Update the stress
        sigma = sigmaHydf*I + s;
    }
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
#ifdef OPENFOAM_NOT_EXTEND
            gMax
            (
                mag
                (
                    DEpsilonPf_.primitiveField()
                  - DEpsilonPf_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().primitiveField()));
#else
            gMax
            (
                mag
                (
                    DEpsilonPf_.internalField()
                  - DEpsilonPf_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().internalField()));
#endif
    }
    else
    {
        return
#ifdef OPENFOAM_NOT_EXTEND
            gMax
            (
                mag
                (
                    DEpsilonP_.primitiveField()
                  - DEpsilonP_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().primitiveField()));
#else
            gMax
            (
                mag
                (
                    DEpsilonP_.internalField()
                  - DEpsilonP_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
#endif
    }
}


void Foam::linearElasticMisesPlastic::updateTotalFields()
{
    // Count cells actively yielding
    int numCellsYielding = 0;

#ifdef OPENFOAM_NOT_EXTEND
    forAll(activeYield_.primitiveField(), celli)
    {
        if (DEpsilonPEq_.primitiveField()[celli] > SMALL)
        {
            activeYield_.primitiveFieldRef()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.primitiveFieldRef()[celli] = 0.0;
        }
    }
#else
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
#endif

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchi)
    {
        if (!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
                {
#ifdef OPENFOAM_NOT_EXTEND
                    activeYield_.boundaryFieldRef()[patchi][facei] = 1.0;
#else
                    activeYield_.boundaryField()[patchi][facei] = 1.0;
#endif
                }
                else
                {
#ifdef OPENFOAM_NOT_EXTEND
                    activeYield_.boundaryFieldRef()[patchi][facei] = 0.0;
#else
                    activeYield_.boundaryField()[patchi][facei] = 0.0;
#endif
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    const int nTotalCells = returnReduce(mesh().nCells(), sumOp<int>());

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << nl
        << "    " << numCellsYielding << " cells ("
        << 100.0*scalar(numCellsYielding)/scalar(nTotalCells)
        << "% of the cells in this material) are actively yielding"
        << nl << endl;

    // Write out magnitude of plastic strain
    // if (mesh().time().outputTime())
    // {
    //     volScalarField epsilonPMag
    //     (
    //         IOobject
    //         (
    //             "epsilonPMag",
    //             mesh().time().timeName(),
    //             mesh(),
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         sqrt((2.0/3.0)*magSqr(dev(epsilonP_)))
    //     );

    //     epsilonPMag.write();
    // }
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
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon()))));

    // Take reference to internal fields
#ifdef OPENFOAM_NOT_EXTEND
    const symmTensorField& DEpsilonPI = DEpsilonP_.primitiveField();
    const symmTensorField& plasticNI = plasticN_.primitiveField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().primitiveField();
    const scalarField& epsilonEqI = epsilonEq.primitiveField();
#else
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();
#endif

    // Calculate error field
    const symmTensorField DEpsilonPErrorI
    (
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL)
    );

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
