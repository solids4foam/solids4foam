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

\*---------------------------------------------------------------------------*/
#include <iostream>
#include<stdio.h>
#include "neoHookeanMisesGTNDamage.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include <Eigen/Dense>
#include "zeroGradientFvPatchFields.H"

using namespace Eigen;
using namespace std;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanMisesGTNDamage, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanMisesGTNDamage, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar neoHookeanMisesGTNDamage::LoopTol_ = 1e-4;

    // Maximum number of iterations for Newton loop
    label neoHookeanMisesGTNDamage::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar neoHookeanMisesGTNDamage::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar neoHookeanMisesGTNDamage::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::symmTensor Foam::neoHookeanMisesGTNDamage::logm(const symmTensor& T)
{
    // Finds the matrix log of the tensor T

    // Convert T to a MatrixXd
    Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Vector3d EVals;
    EVals = es.eigenvalues();

    Vector3d LogEVals;
    LogEVals(0) = Foam::log(EVals(0));
    LogEVals(1) = Foam::log(EVals(1));
    LogEVals(2) = Foam::log(EVals(2));

    Matrix3d D = LogEVals.asDiagonal();

    Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}
Foam::symmTensor Foam::neoHookeanMisesGTNDamage::expm(const symmTensor& T)
{
    // Calculate exponetial of a tensor

    // Convert T to a MatrixXd
    Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Vector3d EVals;
    EVals = es.eigenvalues();

    Vector3d ExpEVals;
    ExpEVals(0) = Foam::exp(EVals(0));
    ExpEVals(1) = Foam::exp(EVals(1));
    ExpEVals(2) = Foam::exp(EVals(2));

    Matrix3d D = ExpEVals.asDiagonal();

    Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}
void Foam::neoHookeanMisesGTNDamage::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanMisesGTNDamage::makeRelF()")
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


Foam::volTensorField& Foam::neoHookeanMisesGTNDamage::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}

void Foam::neoHookeanMisesGTNDamage::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanMisesGTNDamage::makeJ()")
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


Foam::volScalarField& Foam::neoHookeanMisesGTNDamage::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}

void Foam::neoHookeanMisesGTNDamage::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::neoHookeanMisesGTNDamage::makeF()")
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


Foam::volTensorField& Foam::neoHookeanMisesGTNDamage::F1()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}

void Foam::neoHookeanMisesGTNDamage::calculatef
(
    const scalar& DHyd, //specherical plastic strain increment
    const scalar& DqEpsilonP, //deviatroic plastic strain increment
    scalar& f, //porosity
    scalar& fStar, //effective porosity
    const scalar& fOld, //old prosity value
    const scalar& fNonLocal, //non-local porosity from previous iteration
    const scalar& epsilonPEqOld, //old equivalent plastic strain value
    const scalar& DEpsilonPEq, //equivalent plastic strain increment
    const symmTensor& DEpsilonP, //icnrement of plastic strain
    const symmTensor& devT,  //deviatoric stress tensor
    const scalar& p //pressure
)
{
    scalar A = 0;

    // Calaculate void growth due to nucleation if pressure is positive
    if (p > 0)
    {
        A =
            (fN_.value()/(sN_.value()*sqrt(2*3.14)))
           *pow
            (
                2.7182818285,
               -0.5*pow
                (
                    (
                        epsilonPEqOld + DEpsilonPEq - epsilonN_.value()
                    )/sN_.value(),
                    2.0
                )
            );
    }

    scalar fShear = 0;

    //option to calculate shear effects
    if (includeShear_)
    {
        scalar detS = det(devT);
        scalar q1 = sqrt(3.0/2.0)*mag(devT);
        fShear =
            fStar*(1 - pow(27*detS/(2.0*pow(q1,3.0)), 2.0))
           *(devT && DEpsilonP)/q1;
    }

    if (fShear < 0)
    {
        fShear = 0;
    }

    f = fOld + (1 - f)*DHyd + A*DEpsilonPEq + kw_.value()*fShear;

    fStar = fNonLocal;

    if (fNonLocal >= fC_.value() && includeCoalescence_)
    {
        fStar =
            fC_.value()
          + (f - fC_.value())
           *(1/q1_.value() - fC_.value())/(fF_.value() - fC_.value());
    }
}


void Foam::neoHookeanMisesGTNDamage::smallStrainReturnMap
(
    const scalar& pTrial, //trial pressure
    const symmTensor& sTrial, //trial deviatoric sress
    const scalar& qTrial, //trial equivalent plastic stress
    scalar& DHyd, //specherical plastic strain increment
    scalar& DqEpsilonP, //deviatoric plastic strain increment
    const scalar& fNonLocal, //non-local porosity
    const scalar& fStar,  //effective porosity
    const scalar& epsilonPEqOld, //old equivalent plastic strain value
    scalar& DEpsilonPEq, //equivalent plastic strain increment
    symmTensor& plasticN //normal to the yield surface
) const
{
    // Initialise pressure and equivalent stress
    scalar p = 0.0;
    scalar q = 0.0;

    // Set yield stress
    scalar sigmaY = stressPlasticStrainSeries_(epsilonPEqOld);

    if (qTrial > 0)
    {
        plasticN = (3.0/2.0)*(sTrial/qTrial);
    }

    // GTN yield condition
    scalar fTrial =
        pow(qTrial/sigmaY, 2.0)
      + 2*q1_.value()*fStar*cosh((3*pTrial*q2_.value())/(2*sigmaY))
      - (1 + q3_.value()*pow(fStar,2.0));

    if (debug)
    {
        Info<< "pTrial: " << pTrial << nl
            << "qTrial: "<< qTrial << nl
            << "epsilonPEqOld1: "<< epsilonPEqOld << nl
            << "fNonLocal: " << fNonLocal << nl
            << "fTrial: " << fTrial << endl;
    }

    if (fTrial > 0.0)
    {
        // Set deviatoric, speherical and equivalent plastic strain increments
        DqEpsilonP = 0;
        DHyd = 0;
        DEpsilonPEq = 0;

        // Newton loop to solve for speherical, deviatoric and equivalent plastic strain increments
        for (int j = 0; j <= 1000; j++)
        {
            // Update pressure and equivalent stress
            p = pTrial - K_.value()*DHyd;
            q = qTrial - 3*mu_.value()*DqEpsilonP;

            // Update yield stress
            sigmaY = stressPlasticStrainSeries_(epsilonPEqOld+DEpsilonPEq);

            // Initialse fcosh and fsinh - these are components of the GTN yield equation
            scalar fCosh = cosh(3*q2_.value()*p/(2*sigmaY));
            scalar fSinh = sinh(3*q2_.value()*p/(2*sigmaY));

            // Derivaitive of the yield stress
            scalar dSigmaY =
                (
                    stressPlasticStrainSeries_(epsilonPEqOld+DEpsilonPEq + 1e-8)
                  - stressPlasticStrainSeries_(epsilonPEqOld+DEpsilonPEq)
                )/1e-8;

            // These derivatives are components of the derivatives which fill
            // the Jacobian matrix
            scalar df1dq = 2.0*(q/pow(sigmaY,2.0));
            scalar df1dp = 3*fStar*q1_.value()*q2_.value()*fSinh/sigmaY;
            scalar df1dSigmaY =- 2.0*pow(q,2.0)/pow(sigmaY,3.0) - df1dp*p/sigmaY;
            scalar df1dpdp = 9.0*q1_.value()*pow(q2_.value(),2.0)*fStar*fCosh/(2.0*pow(sigmaY,2.0));
            scalar df1dpdSigmaY = -9.0*q1_.value()*pow(q2_.value(),2.0)*
                                   fStar*p*fCosh/(2.0*pow(sigmaY,3.0))
                                  -3*fStar*q1_.value()*q2_.value()*fSinh/pow(sigmaY,2.0);
            scalar df1dqdSigmaY = -4.0*(q/pow(sigmaY,3.0));

            // Calcualte derivatives which constitute the Jacobian matrix
            scalar a1 = -3.0*mu_.value()*df1dq;
            scalar a2 = -K_.value()*df1dp;
            scalar a3 = df1dSigmaY*dSigmaY;

            scalar b1 = df1dp+3.0*DqEpsilonP*mu_.value()*2.0/pow(sigmaY,2.0);
            scalar b2 = -DqEpsilonP*df1dpdp*K_.value() - df1dq;
            scalar b3 = DqEpsilonP*df1dpdSigmaY*dSigmaY - DHyd*df1dqdSigmaY*dSigmaY;

            scalar c1 = -q+3.0*mu_.value()*DqEpsilonP;
            scalar c2 = -p+K_.value()*DHyd;
            scalar c3 = (1-fNonLocal)*(sigmaY+dSigmaY*DEpsilonPEq);

            if (debug)
            {
                Info<< "a1= " << a1 << nl
                    << "a2= " << a2 << nl
                    << "a3= " << a3 << nl
                    << "b1= " << b1 << nl
                    << "b2= " << b2 << nl
                    << "b3= " << b3 << nl
                    << "c1= " << c1 << nl
                    << "c2= " << c2 << nl
                    << "c3= " << c3 << endl;
           }

            // Calculate values of fucnctions f1, f2 and f3
            scalar f1val = pow(q/sigmaY, 2.0) + 2*q1_.value()*fStar*fCosh - (1+q3_.value()*pow(fStar,2.0));
            scalar f2val = DqEpsilonP*df1dp-DHyd*df1dq;
            scalar f3val = (1-fNonLocal)*sigmaY*DEpsilonPEq - q*DqEpsilonP-p*DHyd;

            //create 3x3 Jacobian matrix
            Matrix3d J;
            J   << a1, a2, a3,
                   b1, b2, b3,
                   c1, c2, c3;

            if (debug)
            {
                cout<< "The inverse of J is:" << J.inverse() << nl
                    << "J test:" << J.inverse()*J << nl << nl;
            }

            // Create 3 dimensional vector
            Vector3d fx(f1val,f2val,f3val);

            // Find inverse of matrix
            Matrix3d Jinv = J.inverse();

            // Gauss-newton step
            Vector3d hgn = -Jinv*fx;

            // Update values of deviatoric, spherical and equivalent plastic
            // strain increments
            DqEpsilonP = DqEpsilonP + (hgn(0));
            DHyd = DHyd + (hgn(1));
            DEpsilonPEq = DEpsilonPEq + (hgn(2));

            if (debug)
            {
                Info<< "f1val= " << f1val << nl
                    << "f2val= " << f2val << nl
                    << "f3val= " << f3val << nl
                    << "h0: " << hgn(0) << nl
                    << "h1: " << hgn(1) << nl
                    << "h2: " << hgn(2)<< nl
                    << "DqEpsilonP: " << DqEpsilonP << nl
                    << "DHyd: " << DHyd << nl
                    << "DEpsilonPEq: " << DEpsilonPEq<< endl;
            }

            if (j == 1000)
            {
                WarningIn("neoHookeanMisesGTNDamage::newtonLoop()")
                    << "Plasticity Newton loop not converging" << endl;
            }

            if
            (
                mag(f1val) < LoopTol_
             && mag(f2val) < LoopTol_
             && mag(f3val) < LoopTol_
            )
            {
                break;
            }
        }
    }
    else
    {
        DqEpsilonP = 0;
        DHyd = 0;
        DEpsilonPEq = 0;
    }
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanMisesGTNDamage::Ibar
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


Foam::tmp<Foam::surfaceScalarField> Foam::neoHookeanMisesGTNDamage::Ibar
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
Foam::neoHookeanMisesGTNDamage::neoHookeanMisesGTNDamage
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    relaxationFactor_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    includeShear_
    (
        dict.lookupOrDefault<Switch>("includeShear", false)
    ),
    includeCoalescence_
    (
        dict.lookupOrDefault<Switch>("includeCoalescence", false)
    ),
    charLength_("zero", dimLength, 0.0),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    q1_("zero", dimless, 0.0),
    q2_("zero", dimless, 0.0),
    q3_("zero", dimless, 0.0),
    f0_("zero", dimless, 0.0),
    fN_("zero", dimless, 0.0),
    epsilonN_("zero", dimless, 0.0),
    sN_("zero", dimless, 0.0),
    fC_("zero", dimless, 0.0),
    fU_("zero", dimless, 0.0),
    fF_("zero", dimless, 0.0),
    kw_("zero", dimless, 0.0),
    relFPtr_(NULL),
    JPtr_(NULL),
    FPtr_(NULL),
    stressPlasticStrainSeries_(dict),
    f_
    (
        IOobject
        (
            "f",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    fStar_
    (
        IOobject
        (
            "fStar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    fNonLocal_
    (
        IOobject
        (
            "fNonLocal",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0),
        "zeroGradient"
    ),
    gradfNonLocal_
    (
        IOobject
        (
            "grad(fNonLocal)",
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
    //nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
  // Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
    bEbar_.oldTime();
    epsilonEl_.oldTime();
    f_.oldTime();
    fStar_.oldTime();

    //read the caracteristic length
    charLength_ = dimensionedScalar(dict.lookup("charLength"));

    //read GTN parameters
    q1_ = dimensionedScalar(dict.lookup("q1"));
    q2_ = dimensionedScalar(dict.lookup("q2"));
    q3_ = dimensionedScalar(dict.lookup("q3"));
    fN_ = dimensionedScalar(dict.lookup("fN"));
    epsilonN_ = dimensionedScalar(dict.lookup("epsilonN"));
    sN_ = dimensionedScalar(dict.lookup("sN"));
    fC_ = dimensionedScalar(dict.lookup("fC"));
    f0_ = dimensionedScalar(dict.lookup("f0"));
    fU_ = dimensionedScalar(dict.lookup("fU"));
    fF_ = dimensionedScalar(dict.lookup("fF"));
    kw_ = dimensionedScalar(dict.lookup("kw"));

    //intial porosity
    f_=f0_;


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
            "neoHookeanMisesGTNDamage::neoHookeanMisesGTNDamage::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }



    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanMisesGTNDamage::~neoHookeanMisesGTNDamage()
{
    deleteDemandDrivenData(relFPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(FPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanMisesGTNDamage::rho() const
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
Foam::neoHookeanMisesGTNDamage::impK() const
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


void Foam::neoHookeanMisesGTNDamage::correct(volSymmTensorField& sigma)
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
        F1() = relF() & F1().oldTime();
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
            F1() = F1().oldTime() + gradDD.T();

            // Update the relative deformation gradient
            relF() = F1() & inv(F1().oldTime());
        }
        else
        {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            // Update the total deformation gradient
            F1() = I + gradD.T();

            // Update the relative deformation gradient
            relF() = F1() & inv(F1().oldTime());
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

    // These fields are useful for information and debugging purposes but are not
    // used in this version of the class
    volScalarField dHydField(epsilonPEq_);
    volScalarField dqField(epsilonPEq_);
    scalarField& dHydFieldI = dHydField;
    scalarField& dqFieldI = dqField;

    // Take references to the internal fields for efficiency
    const tensorField& FI = F1();
    symmTensorField epsilonElOldI(epsilonEl_.oldTime());
    symmTensorField& epsilonElI = epsilonEl_;
    symmTensorField& DEpsilonPI = DEpsilonP_;
    symmTensorField& tauI = tau_;
    symmTensorField& sigmaI = sigma;
    tensorField  relFI = relF();

    scalarField epsilonPEqOldI(epsilonPEq_.oldTime());
    scalarField& DEpsilonPEqI = DEpsilonPEq_;
    scalarField fOldI(f_.oldTime());
    scalarField& fI = f_;
    scalarField& fStarI = fStar_;
    scalarField& fStarOldI = fStar_.oldTime();
    scalarField& fNonLocalI = fNonLocal_;

    symmTensor BeOld = I; // old value of the left Cachy-Green tensor
    symmTensor Be = I; // left Cachy-Green tensor
    symmTensor Ee = symmTensor::zero; // elastic strain tensor
    tensor relf = I; // relative deformation gradient
    symmTensor plasticN = symmTensor::zero; // normal to the plastic yielding

    symmTensor sTrial = I; // initialise deviatoric stress
    scalar DHyd, DqEpsilonP = 0.0; // initialise spherical and deviatoric plastic strain
    scalar pTrial, p = 0.0; // initialise trial pressure and pressure
    scalar q, qTrial = 0.0; // initialise equivalent stress and trial equivalent stress

    // Calculate DHyd, DqEpsilonP, DEpsilonPEq and plasticN
    forAll(relFI, cellI)
    {

        //pre-processing step
        relf = relFI[cellI];
        BeOld = expm(2.0*epsilonElOldI[cellI]);
        Be = symm(relf & BeOld & relf.T());
        Ee = 0.5*logm(Be);

        pTrial = K_.value()*tr(Ee);
        sTrial = 2*mu_.value()*dev(Ee);
        qTrial = sqrt(3.0/2.0)*mag(sTrial);

        smallStrainReturnMap
        (
            pTrial,
            sTrial,
            qTrial,
            DHyd,
            DqEpsilonP,
            fNonLocalI[cellI],
            fStarI[cellI],
            epsilonPEqOldI[cellI],
            DEpsilonPEqI[cellI],
            plasticN
        );

        //post-processing step
        DEpsilonPI[cellI] = (1.0/3.0)*DHyd*I + DqEpsilonP*plasticN;
        epsilonElI[cellI] = Ee - DEpsilonPI[cellI];
        p = pTrial - K_.value()*DHyd;
        q = qTrial - 3*mu_.value()*DqEpsilonP;

        tauI[cellI] = p*I + (sTrial - 2*mu_.value()*dev(DEpsilonPI[cellI]));
        sigmaI[cellI] = (1/det(FI[cellI]))*tauI[cellI];
        dHydFieldI[cellI] = DHyd;
        dqFieldI[cellI] = DqEpsilonP;


        //update porosity if cell is undergoing plastic yielding
        if (DEpsilonPEqI[cellI] > 0)
        {
            symmTensor tDev = (q/qTrial)*sTrial;
            calculatef
            (
                DHyd,
                DqEpsilonP,
                fI[cellI],
                fStarI[cellI],
                fOldI[cellI],
                fNonLocalI[cellI],
                epsilonPEqOldI[cellI],
                DEpsilonPEqI[cellI],
                DEpsilonPI[cellI],
                tDev,
                p
            );
         }
         else
         {
             fStarI[cellI] = fStarOldI[cellI];
             fI[cellI] = fOldI[cellI];
         }
    }

    forAll(F1().boundaryField(), patchI)
    {
#ifdef OPENFOAM_NOT_EXTEND
        scalarField& dHydFieldP = dHydField.boundaryFieldRef()[patchI];
        scalarField& dqFieldP = dqField.boundaryFieldRef()[patchI];

        const tensorField& FP = F1().boundaryFieldRef()[patchI];
        symmTensorField epsilonElOldP = epsilonEl_.oldTime().boundaryFieldRef()[patchI];
        symmTensorField& epsilonElP = epsilonEl_.boundaryFieldRef()[patchI];
        symmTensorField& DEpsilonPP = DEpsilonP_.boundaryFieldRef()[patchI];
        symmTensorField tauP = tau_.boundaryFieldRef()[patchI];
        symmTensorField& sigmaP = sigma.boundaryFieldRef()[patchI];
        tensorField  relFP = relF().boundaryFieldRef()[patchI];

        scalarField epsilonPEqOldP = epsilonPEq_.oldTime().boundaryFieldRef()[patchI];
        scalarField& DEpsilonPEqP = DEpsilonPEq_.boundaryFieldRef()[patchI];
        scalarField fOldP = f_.oldTime().boundaryFieldRef()[patchI];
        scalarField& fP = f_.boundaryFieldRef()[patchI];
        scalarField& fStarP = fStar_.boundaryFieldRef()[patchI];
        scalarField& fStarOldP = fStar_.oldTime().boundaryFieldRef()[patchI];
        scalarField& fNonLocalP = fNonLocal_.boundaryFieldRef()[patchI];
#else
        scalarField& dHydFieldP = dHydField.boundaryField()[patchI];
        scalarField& dqFieldP = dqField.boundaryField()[patchI];

        const tensorField& FP = F1().boundaryField()[patchI];
        symmTensorField epsilonElOldP = epsilonEl_.oldTime().boundaryField()[patchI];
        symmTensorField& epsilonElP = epsilonEl_.boundaryField()[patchI];
        symmTensorField& DEpsilonPP = DEpsilonP_.boundaryField()[patchI];
        symmTensorField tauP = tau_.boundaryField()[patchI];
        symmTensorField& sigmaP = sigma.boundaryField()[patchI];
        tensorField  relFP = relF().boundaryField()[patchI];

        scalarField epsilonPEqOldP = epsilonPEq_.oldTime().boundaryField()[patchI];
        scalarField& DEpsilonPEqP = DEpsilonPEq_.boundaryField()[patchI];
        scalarField fOldP = f_.oldTime().boundaryField()[patchI];
        scalarField& fP = f_.boundaryField()[patchI];
        scalarField& fStarP = fStar_.boundaryField()[patchI];
        scalarField& fStarOldP = fStar_.oldTime().boundaryField()[patchI];
        scalarField& fNonLocalP = fNonLocal_.boundaryField()[patchI];
#endif

        // Calculate DHyd, DqEpsilonP, DEpsilonPEq and plasticN
        forAll(FP, faceI)
        {
            // pre-processing step
            relf = relFP[faceI];
            BeOld = expm(2.0*epsilonElOldP[faceI]);
            Be = symm(relf & BeOld & relf.T());
            Ee = 0.5*logm(Be);

            pTrial = K_.value()*tr(Ee);
            sTrial = 2*mu_.value()*dev(Ee);
            qTrial = sqrt(3.0/2.0)*mag(sTrial);

            smallStrainReturnMap
            (
                pTrial,
                sTrial,
                qTrial,
                DHyd,
                DqEpsilonP,
                fNonLocalP[faceI],
                fStarP[faceI],
                epsilonPEqOldP[faceI],
                DEpsilonPEqP[faceI],
                plasticN
            );

            //post-processing step
            DEpsilonPP[faceI] = (1.0/3.0)*DHyd*I + DqEpsilonP*plasticN;
            epsilonElP[faceI] = Ee - (1.0/3.0)*DHyd*I - DqEpsilonP*plasticN;
            p = pTrial - K_.value()*DHyd;
            q = qTrial - 3*mu_.value()*DqEpsilonP;

            tauP[faceI] = p*I + (sTrial-2*mu_.value()*dev(DEpsilonPP[faceI]));
            sigmaP[faceI] = p*I + (sTrial-2*mu_.value()*dev(DEpsilonPP[faceI]));

            dHydFieldP[faceI] = DHyd;
            dqFieldP[faceI] = DqEpsilonP;

            //update porosity if cell is undergoing plastic yielding
            if (DEpsilonPEqP[faceI]>0)
            {
                symmTensor tDev = (q/qTrial)*sTrial;
                calculatef
                (
                    DHyd,
                    DqEpsilonP,
                    fP[faceI],
                    fStarP[faceI],
                    fOldP[faceI],
                    fNonLocalP[faceI],
                    epsilonPEqOldP[faceI],
                    DEpsilonPEqP[faceI],
                    DEpsilonPP[faceI],
                    tDev,
                    p
                );
            }
            else
            {
                fStarP[faceI] = fStarOldP[faceI];
                fP[faceI] = fOldP[faceI];
            }

        }
    }


    // Calcualte non-local porosity
    fvScalarMatrix fEqn
    (
        fvm::Sp(1.0, fNonLocal_)
      - fvm::laplacian(pow(charLength_, 2.0), fNonLocal_)
     == f_
    );

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<scalar>::debug = 0;
#else
    blockLduMatrix::debug = 0;
#endif

    fEqn.solve();

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<scalar>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    gradfNonLocal_ = fvc::grad(fNonLocal_);
}


void Foam::neoHookeanMisesGTNDamage::correct(surfaceSymmTensorField& sigma)
{
   notImplemented("Surface field version of correct not yet implemented");
}


Foam::scalar Foam::neoHookeanMisesGTNDamage::residual()
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


void Foam::neoHookeanMisesGTNDamage::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonP_ += DEpsilonP_;

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


Foam::scalar Foam::neoHookeanMisesGTNDamage::newDeltaT()
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

    F1() = relF() & F1().oldTime();

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon(0.5*log(symm(F1().T() & F1())));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon))));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP_;
    const symmTensorField& plasticNI = plasticN_;
    const symmTensorField& plasticNIold = plasticN_.oldTime();
    const scalarField& epsilonEqI = epsilonEq;

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
                "Foam::scalar Foam::neoHookeanMisesGTNDamage::newDeltaT()"
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
