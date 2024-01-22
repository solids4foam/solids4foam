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

#include "neoHookeanMisesLemaitreDamage.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanMisesLemaitreDamage, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanMisesLemaitreDamage, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar neoHookeanMisesLemaitreDamage::LoopTol_ = 1e-4;

    // Maximum number of iterations for Newton loop
    label neoHookeanMisesLemaitreDamage::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar neoHookeanMisesLemaitreDamage::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar neoHookeanMisesLemaitreDamage::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor Foam::neoHookeanMisesLemaitreDamage::logm(const symmTensor& T)
{
    // Finds the matrix log of the tensor T

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    Eigen::Vector3d LogEVals;
    LogEVals(0) = Foam::log(EVals(0));
    LogEVals(1) = Foam::log(EVals(1));
    LogEVals(2) = Foam::log(EVals(2));

    Eigen::Matrix3d D = LogEVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Eigen::Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}
Foam::symmTensor Foam::neoHookeanMisesLemaitreDamage::expm(const symmTensor& T)
{
    // Calculate exponetial of a tensor

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    Eigen::Vector3d ExpEVals;
    ExpEVals(0) = Foam::exp(EVals(0));
    ExpEVals(1) = Foam::exp(EVals(1));
    ExpEVals(2) = Foam::exp(EVals(2));

    Eigen::Matrix3d D = ExpEVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Eigen::Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}

void Foam::neoHookeanMisesLemaitreDamage::decompose
(
  const symmTensor& tau,
  symmTensor& tauPositive,
  symmTensor& tauNegative
)
{

    // Convert tau to a Matrix3d
    Eigen::Matrix3d TM(3,3);
    TM  <<
        tau.xx(), tau.xy(), tau.xz(),
        tau.xy(), tau.yy(), tau.yz(),
        tau.xz(), tau.yz(), tau.zz();

    //Compute Eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    //Initialise Vectors which will hold positive and negative Eigenvalues respectively
    Eigen::Vector3d EPosVals=EVals;
    Eigen::Vector3d ENegVals=EVals;

    //Set positive and negative Eignvalue vectors
    for (int i=0;i<3;i++)
    {
        if (EVals(i)>0)
        {
            ENegVals(i)=0.0;
        }
        else
        {
            EPosVals(i)=0.0;
        }
    }

    // Create 3x3 Diagonal Matrices
    Eigen::Matrix3d P = EPosVals.asDiagonal();
    Eigen::Matrix3d N = ENegVals.asDiagonal();


    //Compute positive and negative matrices
    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();
    Eigen::Matrix3d resultPos = EVecs*P*EVecs.inverse();
    Eigen::Matrix3d resultNeg = EVecs*N*EVecs.inverse();


    //Fill tauPositive and tauNegative
    tauPositive.xx()=resultPos(0,0);
    tauPositive.xy()=resultPos(0,1);
    tauPositive.xz()=resultPos(0,2);
    tauPositive.yy()=resultPos(1,1);
    tauPositive.yz()=resultPos(1,2);
    tauPositive.zz()=resultPos(2,2);

    tauNegative.xx()=resultNeg(0,0);
    tauNegative.xy()=resultNeg(0,1);
    tauNegative.xz()=resultNeg(0,2);
    tauNegative.yy()=resultNeg(1,1);
    tauNegative.yz()=resultNeg(1,2);
    tauNegative.zz()=resultNeg(2,2);

}

Foam::scalar Foam::neoHookeanMisesLemaitreDamage::macauley
(
  scalar x
)
{
    //macauley brackets

    if (x<0)
    {
        x=0;
    }
    return x;
}

void Foam::neoHookeanMisesLemaitreDamage::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanMisesLemaitreDamage::makeRelF()")
            << "pointer already set" << abort(FatalError);
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


Foam::volTensorField& Foam::neoHookeanMisesLemaitreDamage::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}

void Foam::neoHookeanMisesLemaitreDamage::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanMisesLemaitreDamage::makeJ()")
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


Foam::volScalarField& Foam::neoHookeanMisesLemaitreDamage::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}


void Foam::neoHookeanMisesLemaitreDamage::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanMisesLemaitreDamage::makeF()")
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


Foam::volTensorField& Foam::neoHookeanMisesLemaitreDamage::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


Foam::scalar Foam::neoHookeanMisesLemaitreDamage::calcZheta
(
    const symmTensor& S //deviatoric stress
)
{
    //von Mises equivalent stress
    scalar q = sqrt(3.0/2.0)*mag(S);

    //normalised third invarient
    scalar zheta = (27.0/2.0)*(det(S)/(pow(q,3.0)));

    //These if loops ensure that a floating point exception does not occur when
    //numerical error leads to a zheta value of e.g. 1.000001
    if (zheta > 1.0) {zheta = 1;}
    if (zheta < -1.0) {zheta = -1;}

//    scalar lodeAngle=1-(2.0/3.14)*(acos(zheta));

//    if (lodeAngle>1.0)
//    {
//        lodeAngle=1.0;
//    }
//    if (lodeAngle<-1.0)
//    {
//        lodeAngle=-1.0;
//    }

    return zheta;

}

void Foam::neoHookeanMisesLemaitreDamage::damage
(
    const scalar& DLambda, //plastic multiplier
    const scalar& DEpsilonPEq, //increment of equivalent palstic strain
    const scalar& epsilonPEqOld, //old value of the equivalent plastic strain
    const scalar& triaxiality,  //triaxiality value
    scalar& damage,  //damage
    scalar& damageOld,  //old value of the damage
    const scalar& damageNonLocal,  //non-local damage value
    const symmTensor& tau  //stress tensor
)
{
    scalar epsilonPEq = epsilonPEqOld + DEpsilonPEq;
    scalar sigmaEq = sqrt(3.0/2.0)*mag(dev(tau))/(1 - damageNonLocal);

    // Initialise total energy release rates
    scalar Y = 0;
    scalar Ystar = 0;

    if (triaxiality> (-1.0/3.0) && epsilonPEq > epsilonD_.value())
    {
         // denominator
         const scalar denom =
             (2.0/3.0)*(1.0 + nu_.value())
           + 3.0*(1.0 - 2.0* nu_.value())*pow(triaxiality, 2.0);

         Y= -(pow(sigmaEq, 2.0)/(2.0*E_.value()))*denom;

         scalar sDenom = 0; //denominator for the damage law

         if (lemaitreType_=="classic")
         {
             sDenom = s0_.value();
         }
         else if (lemaitreType_ == "Malcher")
         {
             scalar zheta = calcZheta(dev(tau));
             sDenom =
                 s33_.value()
                /(
                    3.0*mag(triaxiality)
                  + (s33_.value()/s00_.value())*(1 - pow(zheta,2.0))
                 );
         }

         Ystar = pow(-Y/sDenom, b_.value())*DEpsilonPEq;

         // Update damage
         damage = damageOld + Ystar/(1 - damageNonLocal);

         if (damage > damageC_)
         {
             damage = 0.99;
         }
    }
}

void Foam::neoHookeanMisesLemaitreDamage::smallStrainReturnMap
(
    const symmTensor& Ee, //elastic strain
    scalar& sigmaY, //yield stress
    scalar& DLambda, //plasticMultplier
    scalar& damage,  //damage
    const scalar& damageOld,  //old damage value
    const scalar& damageNonLocal,  //non-local damage
    const scalar& epsilonPEqOld, //old equivalent plastic strain
    scalar& activeYield  // check if cell is actively yielding
) const
{
    int i = 0;
    sigmaY = stressPlasticStrainSeries_(epsilonPEqOld);
    scalar qTrial = sqrt(3.0/2.0)*mag(2*mu_.value()*(dev(Ee)));
    scalar fTrial = qTrial - sigmaY;
    DLambda = 0;
    scalar dfTrial = fTrial;

    if (fTrial > 0.0)
    {
        // Newton loop for calcualting plastic increment
        do
        {
            fTrial=
                qTrial - 3*mu_.value()*DLambda/(1 - damageNonLocal)
                - stressPlasticStrainSeries_(epsilonPEqOld+DLambda);

            scalar dSigmaY =
                (
                    stressPlasticStrainSeries_(epsilonPEqOld + DLambda + 1e-8)
                  - stressPlasticStrainSeries_(epsilonPEqOld + DLambda)
                )/1e-8;

            dfTrial = -3*mu_.value()/(1 - damageNonLocal) - dSigmaY;

            DLambda -= fTrial/dfTrial;

            if (i == MaxNewtonIter_)
            {
                WarningIn("neoHookeanElasticMisesPlasticDamage::newtonLoop()")
                    << "Plasticity Newton loop not converging" << endl;
            }
        }
        while ((mag(fTrial) > LoopTol_) && ++i < MaxNewtonIter_);

        activeYield = 1.0;
    }
    else
    {
        DLambda = 0;
        damage = damageOld;
        activeYield = 0.0;
    }
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanMisesLemaitreDamage::Ibar
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

#ifdef OPENFOAM_NOT_EXTEND
    volScalarField& Ibar = tIbar.ref();
#else
    volScalarField& Ibar = tIbar();
#endif

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar;
    const symmTensorField devBEbarI = devBEbar;

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
#ifdef OPENFOAM_NOT_EXTEND
            scalarField& IbarP = Ibar.boundaryFieldRef()[patchI];
#else
            scalarField& IbarP = Ibar.boundaryField()[patchI];
#endif
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


Foam::tmp<Foam::surfaceScalarField> Foam::neoHookeanMisesLemaitreDamage::Ibar
(
    const surfaceSymmTensorField& devBEbar
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

    tmp<surfaceScalarField> tIbar
    (
        new surfaceScalarField
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

#ifdef OPENFOAM_NOT_EXTEND
    surfaceScalarField& Ibar = tIbar.ref();
#else
    surfaceScalarField& Ibar = tIbar();
#endif

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar;
    const symmTensorField devBEbarI = devBEbar;

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
#ifdef OPENFOAM_NOT_EXTEND
            scalarField& IbarP = Ibar.boundaryFieldRef()[patchI];
#else
            scalarField& IbarP = Ibar.boundaryField()[patchI];
#endif
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

#ifdef FOAM_EXTEND
    Ibar.correctBoundaryConditions();
#endif

    return tIbar;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanMisesLemaitreDamage::neoHookeanMisesLemaitreDamage
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    lemaitreType_(dict.lookupOrDefault<string>("lemaitreType", "classic")),
    stressPlasticStrainSeries_(dict),
    charLength_("zero", dimLength, 0.0),
    epsilonD_("zero", dimless, 0.0),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    b_("zero", dimless, 0.0),
    s0_("zero", dimPressure, 0.0),
    s00_("zero", dimPressure, 0.0),
    s33_("zero", dimPressure, 0.0),
    damageC_(dict.lookupOrDefault<scalar>("damageC", 0.99)),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    relFPtr_(NULL),
    JPtr_(NULL),
    FPtr_(NULL),
    damage_
    (
        IOobject
        (
            "damage",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    damageNonLocal_
    (
        IOobject
        (
            "damageNonLocal",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0),
        "zeroGradient"
    ),
    gradDamageNonLocal_
    (
        IOobject
        (
            "grad(damageNonLocal)",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimless/dimLength, vector::zero)
    ),
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
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    smoothPressure_(dict.lookupOrDefault<Switch>("smoothPressure", true)),
    tau_
    (
        IOobject
        (
            "tau",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
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
            "initialYieldStress", dimPressure, 0.0
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
    epsilonEl_
    (
        IOobject
        (
            "epsilonEl",
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
    sigmaPositive_
    (
        IOobject
        (
            "sigmaPositive",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    sigmaNegative_
    (
        IOobject
        (
            "sigmaNegative",
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
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
    ),
    triaxiality_
    (
        IOobject
        (
            "triaxiality",
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
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
    bEbar_.oldTime();
    damage_.oldTime();
        epsilonEl_.oldTime();


    Info<< "    smoothPressure: " << smoothPressure_ << nl
        << "    updateBEbarConsistent: " << updateBEbarConsistent_ << endl;

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
        if (planeStress())
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
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
            "neoHookeanMisesLemaitreDamage::neoHookeanMisesLemaitreDamage::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    if (lemaitreType_=="classic")
    {
        s0_ = dimensionedScalar(dict.lookup("s0"));
        b_ = dimensionedScalar(dict.lookup("b"));
        epsilonD_ = dimensionedScalar(dict.lookup("epsilonD"));
        charLength_ = dimensionedScalar(dict.lookup("charLength"));
    }
    else if (lemaitreType_=="Malcher")
    {
        s00_ = dimensionedScalar(dict.lookup("s00"));
        s33_ = dimensionedScalar(dict.lookup("s33"));
        b_ = dimensionedScalar(dict.lookup("b"));
        epsilonD_ = dimensionedScalar(dict.lookup("epsilonD"));
        charLength_ = dimensionedScalar(dict.lookup("charLength"));
    }

    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanMisesLemaitreDamage::~neoHookeanMisesLemaitreDamage()
{
    deleteDemandDrivenData(relFPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(FPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanMisesLemaitreDamage::rho() const
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
Foam::neoHookeanMisesLemaitreDamage::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    const volScalarField muBar(Ibar*mu_);

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*muBar*DLambda_/magSTrial));

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


void Foam::neoHookeanMisesLemaitreDamage::correct(volSymmTensorField& sigma)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the relative deformation gradient
        relF() = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF() & F().oldTime();
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            F() = F().oldTime() + gradDD.T();

            // Update the relative deformation gradient
            relF() = F() & inv(F().oldTime());
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
    }
    else
    {
        FatalErrorIn
        (
            type() + "::correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }
    DEpsilonP_.storePrevIter();

    // Take references to the internal fields for efficiency
    const tensorField& FI = F();
    symmTensorField epsilonElOldI = epsilonEl_.oldTime();
    symmTensorField& epsilonElI = epsilonEl_;
    // symmTensorField& sigmaPositiveI =sigmaPositive_;
    // symmTensorField& sigmaNegativeI =sigmaNegative_;
    symmTensorField& DEpsilonPI = DEpsilonP_;
    symmTensorField& tauI = tau_;
    symmTensorField& sigmaI = sigma;
    tensorField  relFI = relF();
    // scalarField& DSigmaYI = DSigmaY_;
    scalarField& epsilonPEqOldI = epsilonPEq_.oldTime();
    scalarField& DEpsilonPEqI = DEpsilonPEq_;
    scalarField& DLambdaI = DLambda_;
    scalarField& sigmaYI = sigmaY_;
    scalarField& activeYieldI = activeYield_;

    scalarField& damageNonLocalI = damageNonLocal_;
    scalarField& damageI = damage_;
    scalarField damageOldI = damage_.oldTime();

    // Initialise variables for plasticity model
    symmTensor BeOld = I; // old value of the left Cachy-Green tensor
    symmTensor Be = I; // left Cachy-Green tensor
    symmTensor Ee = symmTensor::zero;  // elastic strain tensor
    tensor relf = I; // relative deformation gradient
    symmTensor plasticN = symmTensor::zero; // normal to the plastic yielding

    // Calculate DLambda and plasticN
    forAll(relFI, cellI)
    {
        // Pre-processing step
        relf = relFI[cellI];
        BeOld = expm(2.0*epsilonElOldI[cellI]);
        Be = symm(relf & BeOld & relf.T());
        Ee = 0.5*logm(Be);

        if (mag(Ee) > 0)
        {
            plasticN = dev(Ee)/(mag(dev(Ee)));
        }
        else
        {
            plasticN = symmTensor::zero;
        }

        smallStrainReturnMap
        (
            Ee,
            sigmaYI[cellI],
            DLambdaI[cellI],
            damageI[cellI],
            damageOldI[cellI],
            damageNonLocalI[cellI],
            epsilonPEqOldI[cellI],
            activeYieldI[cellI]
       );

       // Post-processing step
       DEpsilonPI[cellI] =
           sqrt(3.0/2.0)*DLambdaI[cellI]/(1-damageNonLocalI[cellI])*(plasticN);
       epsilonElI[cellI] = Ee - DEpsilonPI[cellI];
       tauI[cellI] =
           (1 - damageNonLocalI[cellI])
          *(
              2*mu_.value()*dev(epsilonElI[cellI])
            + K_.value()*tr(epsilonElI[cellI])*I
          );
       sigmaI[cellI] = (1/det(FI[cellI]))*tauI[cellI];

       DEpsilonPEqI[cellI] = DLambdaI[cellI];

       // If cell undergoing plastic deformation calculate damage
       // PC: should be check "> 1.0 - SMALL" here?
       if (activeYieldI[cellI] == 1.0)
       {
             scalar sigmaEq = sqrt(3.0/2.0)*mag(dev(tauI[cellI]))/(1-damageNonLocalI[cellI]);
             scalar triaxiality = (1.0/3.0)*(tr(tauI[cellI])/(1-damageNonLocalI[cellI]))/sigmaEq;

             damage
             (
                 DLambdaI[cellI],
                 DEpsilonPEqI[cellI],
                 epsilonPEqOldI[cellI],
                 triaxiality,
                 damageI[cellI],
                 damageOldI[cellI],
                 damageNonLocalI[cellI],
                 tauI[cellI]
             );
        }

    }

    forAll(F().boundaryField(), patchI)
    {

        const tensorField& FP = F().boundaryField()[patchI];
#ifdef OPENFOAM_NOT_EXTEND
        symmTensorField epsilonElOldP = epsilonEl_.oldTime().boundaryFieldRef()[patchI];
        symmTensorField& epsilonElP = epsilonEl_.boundaryFieldRef()[patchI];
        symmTensorField& DEpsilonPP = DEpsilonP_.boundaryFieldRef()[patchI];
        symmTensorField& tauP = tau_.boundaryFieldRef()[patchI];
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
        tensorField  relFP = relF().boundaryFieldRef()[patchI];
        scalarField& epsilonPEqOldP = epsilonPEq_.oldTime().boundaryFieldRef()[patchI];
        scalarField& DEpsilonPEqP = DEpsilonPEq_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryFieldRef()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryFieldRef()[patchI];

        scalarField& damageNonLocalP = damageNonLocal_.boundaryFieldRef()[patchI];
        scalarField& damageP = damage_.boundaryFieldRef()[patchI];
        scalarField damageOldP = damage_.oldTime().boundaryFieldRef()[patchI];
#else
        symmTensorField epsilonElOldP = epsilonEl_.oldTime().boundaryField()[patchI];
        symmTensorField& epsilonElP = epsilonEl_.boundaryField()[patchI];
        symmTensorField& DEpsilonPP = DEpsilonP_.boundaryField()[patchI];
        symmTensorField& tauP = tau_.boundaryField()[patchI];
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
        tensorField  relFP = relF().boundaryField()[patchI];
        scalarField& epsilonPEqOldP = epsilonPEq_.oldTime().boundaryField()[patchI];
        scalarField& DEpsilonPEqP = DEpsilonPEq_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        scalarField& activeYieldP = activeYield_.boundaryField()[patchI];

        scalarField& damageNonLocalP = damageNonLocal_.boundaryField()[patchI];
        scalarField& damageP = damage_.boundaryField()[patchI];
        scalarField damageOldP = damage_.oldTime().boundaryField()[patchI];
#endif

        // Calculate DLambda and plasticN
        forAll(FP, faceI)
        {
            // Pre-processing step
            relf = relFP[faceI];
            BeOld = expm(2.0*epsilonElOldP[faceI]);
            Be = symm(relf & BeOld & relf.T());
            Ee = 0.5*logm(Be);

            if (mag(Ee) > 0)
            {
                plasticN = dev(Ee)/(mag(dev(Ee)));
            }
            else
            {
                plasticN = symmTensor::zero;
            }

            smallStrainReturnMap
            (
                Ee,
                sigmaYP[faceI],
                DLambdaP[faceI],
                damageP[faceI],
                damageOldP[faceI],
                damageNonLocalP[faceI],
                epsilonPEqOldP[faceI],
                activeYieldP[faceI]
            );

            // Post-processing step
            DEpsilonPP[faceI] =
                sqrt(3.0/2.0)*DLambdaP[faceI]*(plasticN)
               /(1-damageNonLocalP[faceI]);

            epsilonElP[faceI] = Ee - DEpsilonPP[faceI];

            tauP[faceI] =
                (1-damageNonLocalP[faceI])
               *(
                   2*mu_.value()*dev(epsilonElP[faceI])
                 + K_.value()*tr(epsilonElP[faceI])*I
               );

            sigmaP[faceI] = (1/det(FP[faceI]))* tauP[faceI];

            DEpsilonPEqP[faceI] = DLambdaP[faceI];

            // If cell undergoing plastic deformation calculate damage
            if (activeYieldP[faceI] == 1.0)
            {
                scalar sigmaEq =
                    sqrt(3.0/2.0)*mag(dev(tauP[faceI]))/(1 - damageNonLocalP[faceI]);
                scalar triaxiality =
                    (1.0/3.0)*(tr(tauP[faceI])/(1 - damageNonLocalP[faceI]))/sigmaEq;

                damage
                (
                    DLambdaP[faceI],
                    DEpsilonPEqP[faceI],
                    epsilonPEqOldP[faceI],
                    triaxiality,
                    damageP[faceI],
                    damageOldP[faceI],
                    damageNonLocalP[faceI],
                    tauP[faceI]
               );
           }
        }
    }

    sigmaHyd_ = (1.0/3.0)*tr(sigma);

    // Calculate non-local damage
    fvScalarMatrix damageEqn
    (
        fvm::Sp(1.0, damageNonLocal_)
      - fvm::laplacian(pow(charLength_, 2.0), damageNonLocal_)
     == damage_
    );

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<scalar>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    damageEqn.solve();

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<scalar>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    gradDamageNonLocal_ = fvc::grad(damageNonLocal_);
}


void Foam::neoHookeanMisesLemaitreDamage::correct(surfaceSymmTensorField& sigma)
{
   notImplemented("wip");
}


Foam::scalar Foam::neoHookeanMisesLemaitreDamage::residual()
{
    // Calculate residual based on change in plastic strain increment
    return
#ifdef OPENFOAM_NOT_EXTEND
    (
        max
        (
            mag
            (
                DEpsilonP_.internalField()
              - DEpsilonP_.prevIter().internalField()
            )()
        )/max(SMALL + mag(DEpsilonP_.prevIter().internalField()))
    ).value();
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


void Foam::neoHookeanMisesLemaitreDamage::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonP_ += DEpsilonP_;

    triaxiality_=sigmaHyd_/(sqrt(3.0/2.0)*mag(dev(tau_)));

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.internalField(), celli)
    {
        if (DEpsilonPEq_.internalField()[celli] > SMALL)
        {
#ifdef OPENFOAM_NOT_EXTEND
            activeYield_.primitiveFieldRef()[celli] = 1.0;
#else
            activeYield_.internalField()[celli] = 1.0;
#endif
            numCellsYielding++;
        }
        else
        {
#ifdef OPENFOAM_NOT_EXTEND
            activeYield_.primitiveFieldRef()[celli] = 0.0;
#else
            activeYield_.internalField()[celli] = 0.0;
#endif
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

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


Foam::scalar Foam::neoHookeanMisesLemaitreDamage::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformatio gradient
    if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    {
        F() = fvc::average(relFf()) & F().oldTime();
    }
    else
    {
        F() = relF() & F().oldTime();
    }

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon(0.5*log(symm(F().T() & F())));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon))));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();

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
                "Foam::scalar Foam::neoHookeanMisesLemaitreDamage::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is lover 50 times larger "
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
