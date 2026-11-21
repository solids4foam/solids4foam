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

#include "linearElasticGTNPlasticMechanicalConstitutiveLaw.H"
#include "misesPlasticDiagnostics.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "SubList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        linearElasticGTNPlasticMechanicalConstitutiveLaw, 1
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        linearElasticGTNPlasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


const Foam::word
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::nonLocalFieldName =
    "voidFractionNonLocal";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::yieldStress
(
    const scalar epsilonPEq
) const
{
    return stressPlasticStrainSeries_(max(epsilonPEq, SMALL));
}


Foam::scalar
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::effectiveVoidFraction
(
    const scalar f
) const
{
    if (coalescence_ && f > fC_)
    {
        return fC_ + (fu_ - fC_)*(f - fC_)/(fF_ - fC_);
    }

    return f;
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::returnMap
(
    const symmTensor& epsilonETrial,
    const scalar epsilonPEq0,
    const scalar f,
    const scalar maxMagDEpsilon,
    pointResult& result
) const
{
    const scalar mu = mu_.value();
    const scalar kappa = kappa_.value();
    const scalar fStar = effectiveVoidFraction(f);

    // Trial stress
    const symmTensor sTrial(2.0*mu*dev(epsilonETrial));
    const scalar magSTrial = mag(sTrial);
    const scalar qTrial = sqrt(3.0/2.0)*magSTrial;
    const scalar pTrial = kappa*tr(epsilonETrial);

    const scalar sigmaY0 = yieldStress(epsilonPEq0);
    const scalar PhiTrial =
        sqr(qTrial/sigmaY0)
      + 2.0*q1_*fStar*cosh(1.5*q2_*pTrial/sigmaY0)
      - 1.0 - q3_*sqr(fStar);

    // The unknowns: the equivalent deviatoric, volumetric and matrix
    // equivalent plastic strain increments
    scalar DEpsilonQ = 0;
    scalar DEpsilonV = 0;
    scalar DEpsilonPEq = 0;

    if (PhiTrial > SMALL)
    {
        // The residuals of the yield condition, the normality of the flow,
        //     DEpsilonV dPhi/dq - DEpsilonQ dPhi/dp = 0
        // multiplied by sigmaY, and the work equivalence, divided by sigmaY;
        // the last two are strains. With the Jacobian, if asked for
        auto residual = [&]
        (
            const scalar dq,
            const scalar dv,
            const scalar dr,
            vector& R,
            tensor* JPtr
        )
        {
            const scalar q = qTrial - 3.0*mu*dq;
            const scalar p = pTrial - kappa*dv;
            const scalar sigmaY = yieldStress(epsilonPEq0 + dr);

            // Bounded against overflow of the hyperbolic functions, which
            // only a wild Newton iterate reaches
            const scalar a = min(max(1.5*q2_*p/sigmaY, -300.0), 300.0);
            const scalar ch = cosh(a);
            const scalar sh = sinh(a);
            const scalar c = 3.0*q1_*q2_*fStar;

            R.x() = sqr(q/sigmaY) + 2.0*q1_*fStar*ch - 1.0 - q3_*sqr(fStar);
            R.y() = 2.0*dv*q/sigmaY - c*dq*sh;
            R.z() = (1.0 - f)*dr - (q*dq + p*dv)/sigmaY;

            if (JPtr)
            {
                tensor& J = *JPtr;

                const scalar H =
                    (
                        yieldStress(epsilonPEq0 + dr + finiteDiff_)
                      - yieldStress(epsilonPEq0 + dr)
                    )/finiteDiff_;

                const scalar dadv = -1.5*q2_*kappa/sigmaY;
                const scalar dadr = -a*H/sigmaY;

                J.xx() = -6.0*mu*q/sqr(sigmaY);
                J.xy() = 2.0*q1_*fStar*sh*dadv;
                J.xz() =
                    -2.0*sqr(q)*H/pow3(sigmaY) + 2.0*q1_*fStar*sh*dadr;

                J.yx() = -6.0*mu*dv/sigmaY - c*sh;
                J.yy() = 2.0*q/sigmaY - c*dq*ch*dadv;
                J.yz() = -2.0*dv*q*H/sqr(sigmaY) - c*dq*ch*dadr;

                J.zx() = -(q - 3.0*mu*dq)/sigmaY;
                J.zy() = -(p - kappa*dv)/sigmaY;
                J.zz() = (1.0 - f) + (q*dq + p*dv)*H/sqr(sigmaY);
            }
        };

        auto merit = [&](const vector& R)
        {
            return
                sqr(R.x())
              + sqr(R.y()/maxMagDEpsilon)
              + sqr(R.z()/maxMagDEpsilon);
        };

        // The equivalent stress cannot be returned past zero
        const scalar maxDEpsilonQ = qTrial/(3.0*mu);

        vector R(vector::zero);
        tensor J(tensor::zero);
        residual(DEpsilonQ, DEpsilonV, DEpsilonPEq, R, &J);

        label iter = 0;

        do
        {
            if (mag(det(J)) < VSMALL)
            {
                WarningInFunction
                    << "Singular Jacobian in the GTN return mapping" << endl;
                break;
            }

            const vector dx(-(inv(J) & R));

            // Backtrack until the residual falls
            const scalar merit0 = merit(R);
            scalar alpha = 1.0;
            scalar dq = 0;
            scalar dv = 0;
            scalar dr = 0;
            vector RNew(vector::zero);

            for (label i = 0; i < 30; ++i)
            {
                dq = min(max(DEpsilonQ + alpha*dx.x(), 0.0), maxDEpsilonQ);
                dv = DEpsilonV + alpha*dx.y();
                dr = max(DEpsilonPEq + alpha*dx.z(), 0.0);

                residual(dq, dv, dr, RNew, nullptr);

                if (merit(RNew) < merit0)
                {
                    break;
                }

                alpha /= 2;
            }

            const scalar magStep =
                max
                (
                    max(mag(dq - DEpsilonQ), mag(dv - DEpsilonV)),
                    mag(dr - DEpsilonPEq)
                );

            DEpsilonQ = dq;
            DEpsilonV = dv;
            DEpsilonPEq = dr;

            residual(DEpsilonQ, DEpsilonV, DEpsilonPEq, R, &J);

            if
            (
                magStep/maxMagDEpsilon < loopTol_
             && mag(R.x()) < loopTol_
             && mag(R.y())/maxMagDEpsilon < loopTol_
             && mag(R.z())/maxMagDEpsilon < loopTol_
            )
            {
                break;
            }
        } while (++iter < maxNewtonIter_);

        if (iter >= maxNewtonIter_)
        {
            WarningInFunction
                << "GTN return mapping Newton loop not converging" << endl;
        }
    }

    // The flow is normal to the yield surface, which is isotropic: the
    // deviatoric part lies along the trial deviatoric stress
    symmTensor n(symmTensor::zero);
    if (magSTrial > SMALL)
    {
        n = sTrial/magSTrial;
    }

    result.DEpsilonP =
        (DEpsilonV/3.0)*symmTensor(I) + (sqrt(3.0/2.0)*DEpsilonQ)*n;
    result.epsilonE = epsilonETrial - result.DEpsilonP;
    result.sigma =
        2.0*mu*dev(result.epsilonE)
      + kappa*tr(result.epsilonE)*symmTensor(I);
    result.DEpsilonV = DEpsilonV;
    result.DEpsilonQ = DEpsilonQ;
    result.DEpsilonPEq = DEpsilonPEq;
    result.p = pTrial - kappa*DEpsilonV;
    result.q = qTrial - 3.0*mu*DEpsilonQ;
    result.f = f;
    result.fYield = f;
    result.sigmaY = yieldStress(epsilonPEq0 + DEpsilonPEq);

    result.theta = 1.0;
    if (qTrial > SMALL)
    {
        result.theta = min(max(result.q/qTrial, 0.0), 1.0);
    }
}


Foam::scalar
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::voidFractionUpdate
(
    const pointResult& result,
    const scalar f0,
    const scalar epsilonPEq0
) const
{
    if (result.DEpsilonPEq <= 0 && result.DEpsilonV == 0)
    {
        return f0;
    }

    // Strain-controlled nucleation, under a tensile mean stress
    scalar A = 0;

    if (fN_ > 0 && result.p > 0)
    {
        const scalar epsilonPEq = epsilonPEq0 + result.DEpsilonPEq;

#ifdef OPENFOAM_NOT_EXTEND
        const scalar pi = constant::mathematical::pi;
#else
        const scalar pi = mathematicalConstant::pi;
#endif

        A =
            fN_/(sN_*sqrt(2.0*pi))
           *exp(-0.5*sqr((epsilonPEq - epsilonN_)/sN_));
    }

    // The shear term, kw f omega DEpsilonQ, which is linear in f
    scalar shear = 0;

    if (kw_ > 0 && result.q > SMALL)
    {
        const scalar xi = 27.0*det(dev(result.sigma))/(2.0*pow3(result.q));
        const scalar omega = min(max(1.0 - sqr(xi), 0.0), 1.0);

        shear = kw_*omega*result.DEpsilonQ;
    }

    // f = f0 + (1 - f) DEpsilonV + A DEpsilonPEq + kw f omega DEpsilonQ
    const scalar denom = max(1.0 + result.DEpsilonV - shear, SMALL);

    return max((f0 + result.DEpsilonV + A*result.DEpsilonPEq)/denom, 0.0);
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::localUpdate
(
    const symmTensor& epsilonETrial,
    const scalar epsilonPEq0,
    const scalar f0,
    const scalar maxMagDEpsilon,
    pointResult& result
) const
{
    // A failed point stays failed
    if (f0 >= fFailed_)
    {
        returnMap(epsilonETrial, epsilonPEq0, fFailed_, maxMagDEpsilon, result);
        return;
    }

    // At the void fraction the step started with. An elastic step there is
    // a solution: it changes nothing
    returnMap(epsilonETrial, epsilonPEq0, f0, maxMagDEpsilon, result);

    const scalar R0 = f0 - voidFractionUpdate(result, f0, epsilonPEq0);

    if (mag(R0) < voidFractionTol_)
    {
        return;
    }

    // The residual f - f(f0, increments at f) is negative at f0 when the
    // voids grow, and positive when they close
    scalar a = 0;
    scalar b = 0;
    scalar Ra = 0;
    scalar Rb = 0;

    pointResult resultBound;

    if (R0 < 0)
    {
        // The root nearest f0, which the step reaches continuously: bracket
        // it by doubling the growth from f0, starting from the explicit
        // growth, -R0. A root further on may exist where nearer ones do,
        // so the bracket [f0, fFailed] is not enough. If the residual is
        // still negative at the failed void fraction, the step takes the
        // point past it, and it fails
        a = f0;
        Ra = R0;
        scalar Df = -R0;

        while (true)
        {
            b = min(f0 + 2*Df, fFailed_);

            returnMap
            (
                epsilonETrial, epsilonPEq0, b, maxMagDEpsilon, resultBound
            );

            Rb = b - voidFractionUpdate(resultBound, f0, epsilonPEq0);

            if (Rb >= 0)
            {
                break;
            }

            if (b >= fFailed_)
            {
                result = resultBound;
                return;
            }

            a = b;
            Ra = Rb;
            Df = b - f0;
        }
    }
    else
    {
        // With no voids the flow is isochoric, so the void fraction can only
        // nucleate there, and the residual is negative
        returnMap(epsilonETrial, epsilonPEq0, 0, maxMagDEpsilon, resultBound);

        a = 0;
        Ra = -voidFractionUpdate(resultBound, f0, epsilonPEq0);
        b = f0;
        Rb = R0;

        if (Ra >= 0)
        {
            result = resultBound;
            return;
        }
    }

    // Root in [a, b], by regula falsi with the Illinois modification
    label side = 0;

    for (label iter = 0; iter < maxVoidFractionIter_; ++iter)
    {
        const scalar f = (a*Rb - b*Ra)/(Rb - Ra);

        returnMap(epsilonETrial, epsilonPEq0, f, maxMagDEpsilon, result);

        const scalar R = f - voidFractionUpdate(result, f0, epsilonPEq0);

        if (mag(R) < voidFractionTol_ || (b - a) < voidFractionTol_)
        {
            return;
        }

        if (R < 0)
        {
            a = f;
            Ra = R;

            if (side == -1)
            {
                Rb /= 2;
            }
            side = -1;
        }
        else
        {
            b = f;
            Rb = R;

            if (side == 1)
            {
                Ra /= 2;
            }
            side = 1;
        }
    }

    WarningInFunction
        << "Coupled void fraction update not converging" << endl;
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::nonLocalUpdate
(
    const symmTensor& epsilonETrial,
    const scalar epsilonPEq0,
    const scalar f0,
    const scalar fNonLocal,
    const scalar maxMagDEpsilon,
    pointResult& result
) const
{
    // The non-local void fraction is a smoothed local one, so it lies in
    // [0, fFailed] up to the tolerance of its linear solver
    const scalar fYield = min(max(fNonLocal, 0.0), fFailed_);

    returnMap(epsilonETrial, epsilonPEq0, fYield, maxMagDEpsilon, result);

    // Explicit in the local void fraction, which the stress does not
    // depend on
    scalar f = fFailed_;

    if (f0 < fFailed_)
    {
        f = min(voidFractionUpdate(result, f0, epsilonPEq0), fFailed_);
    }

    result.f = f;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::
linearElasticGTNPlasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_("mu", dimPressure, 0.0),
    kappa_("kappa", dimPressure, 0.0),
    stressPlasticStrainSeries_(dict),
    q1_(dict.lookupOrDefault<scalar>("q1", 1.5)),
    q2_(dict.lookupOrDefault<scalar>("q2", 1.0)),
    q3_(dict.lookupOrDefault<scalar>("q3", sqr(q1_))),
    f0_(readScalar(dict.lookup("f0"))),
    fN_(dict.lookupOrDefault<scalar>("fN", 0.0)),
    sN_(dict.lookupOrDefault<scalar>("sN", 0.0)),
    epsilonN_(dict.lookupOrDefault<scalar>("epsilonN", 0.0)),
    kw_(dict.lookupOrDefault<scalar>("kw", 0.0)),
    coalescence_(dict.lookupOrDefault<Switch>("coalescence", false)),
    fC_(0),
    fF_(0),
    fu_(0),
    fFailed_(0),
    nonLocal_(false),
    characteristicLength_(0),
    loopTol_(dict.lookupOrDefault<scalar>("NewtonLoopTol", 1e-10)),
    maxNewtonIter_(dict.lookupOrDefault<label>("NewtonMaxIter", 100)),
    finiteDiff_(dict.lookupOrDefault<scalar>("NewtonFiniteDiffEps", 0.25e-6)),
    voidFractionTol_(dict.lookupOrDefault<scalar>("voidFractionTol", 1e-12)),
    maxVoidFractionIter_
    (
        dict.lookupOrDefault<label>("voidFractionMaxIter", 100)
    )
{
    if (E_.dimensions() != dimPressure || nu_.dimensions() != dimless)
    {
        FatalIOErrorInFunction(dict)
            << "E must have dimensions " << dimPressure << " and nu must be "
            << "dimensionless" << exit(FatalIOError);
    }

    if (nu_.value() >= 0.5 - SMALL || nu_.value() <= -1.0)
    {
        FatalIOErrorInFunction(dict)
            << "Poisson's ratio nu = " << nu_.value() << " must lie in "
            << "(-1, 0.5)" << exit(FatalIOError);
    }

    // Note: the planeStress entry is injected into this dictionary by the
    // mechanicalConstitutiveLawManager
    if (Switch(dict.lookupOrDefault<Switch>("planeStress", false)))
    {
        FatalIOErrorInFunction(dict)
            << "Not implemented for planeStress. If needed, you can solve the "
            << "case in 3-D and set the back to a symmetryPlane and the front "
            << "to traction-free" << exit(FatalIOError);
    }

    mu_ = E_/(2.0*(1.0 + nu_));
    kappa_ = E_/(3.0*(1.0 - 2.0*nu_));

    if (q1_ <= 0 || q2_ <= 0 || q3_ <= 0 || q3_ > sqr(q1_))
    {
        FatalIOErrorInFunction(dict)
            << "q1 and q2 must be positive, and q3 must lie in (0, q1^2]"
            << exit(FatalIOError);
    }

    // The effective void fraction at which the yield surface shrinks to a
    // point: the smaller root of 2 q1 f* - 1 - q3 f*^2
    fu_ = (q1_ - sqrt(sqr(q1_) - q3_))/q3_;

    // A failed point keeps about 1% of its strength, as for q3 = q1^2 the
    // yield stress in shear is (1 - q1 f*) sigmaY
    const scalar fStarFailed = 0.99*fu_;
    fFailed_ = fStarFailed;

    if (coalescence_)
    {
        fC_ = readScalar(dict.lookup("fC"));
        fF_ = readScalar(dict.lookup("fF"));

        if (fC_ <= 0 || fF_ <= fC_ || fC_ >= fStarFailed)
        {
            FatalIOErrorInFunction(dict)
                << "The coalescence void fractions must satisfy 0 < fC < fF "
                << "and fC < 0.99 fu = " << fStarFailed << exit(FatalIOError);
        }

        fFailed_ = fC_ + (fStarFailed - fC_)*(fF_ - fC_)/(fu_ - fC_);
    }

    if (fFailed_ >= 1)
    {
        FatalIOErrorInFunction(dict)
            << "The void fraction of a failed point, " << fFailed_ << ", "
            << "must be below 1" << exit(FatalIOError);
    }

    if (f0_ < 0 || f0_ >= fFailed_)
    {
        FatalIOErrorInFunction(dict)
            << "The initial void fraction f0 = " << f0_ << " must lie in "
            << "[0, " << fFailed_ << ")" << exit(FatalIOError);
    }

    if (fN_ < 0 || kw_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "fN and kw must not be negative" << exit(FatalIOError);
    }

    if (fN_ > 0 && (sN_ <= 0 || !dict.found("epsilonN")))
    {
        FatalIOErrorInFunction(dict)
            << "Nucleation, with fN > 0, needs a positive sN and epsilonN"
            << exit(FatalIOError);
    }

    // The regularisation
    const word nonLocal(dict.lookupOrDefault<word>("nonLocal", "none"));

    if (nonLocal == "classic")
    {
        nonLocal_ = true;

        const dimensionedScalar l(dict.lookup("characteristicLength"));

        if (l.dimensions() != dimLength || l.value() < 0)
        {
            FatalIOErrorInFunction(dict)
                << "characteristicLength must not be negative, with "
                << "dimensions " << dimLength << exit(FatalIOError);
        }

        characteristicLength_ = l.value();
    }
    else if (nonLocal != "none")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown nonLocal " << nonLocal << ". Valid forms are none "
            << "and classic" << exit(FatalIOError);
    }

    Info<< "    GTN: q1 = " << q1_ << ", q2 = " << q2_ << ", q3 = " << q3_
        << ", f0 = " << f0_ << ", fN = " << fN_;

    if (fN_ > 0)
    {
        Info<< ", sN = " << sN_ << ", epsilonN = " << epsilonN_;
    }

    Info<< ", kw = " << kw_ << ", coalescence = " << Switch(coalescence_);

    if (coalescence_)
    {
        Info<< ", fC = " << fC_ << ", fF = " << fF_;
    }

    Info<< ", failed void fraction = " << fFailed_
        << ", nonLocal = " << nonLocal;

    if (nonLocal_)
    {
        Info<< ", characteristicLength = " << characteristicLength_;
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    typedef mechanicalConstitutiveLawStateSpec::stateRole role;

    spec.addSymmTensor("epsilonP", role::persistent, symmTensor::zero);

    // The equivalent plastic strain of the matrix, which drives the
    // hardening
    spec.addScalar("epsilonPEq", role::persistent, 0.0);

    spec.addScalar("sigmaY", role::persistent, stressPlasticStrainSeries_(0.0));

    spec.addScalar("voidFraction", role::persistent, f0_);
}


Foam::scalar
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::
smallStrainConvergenceScale
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawState&
) const
{
    const UIndirectList<tensor>& gradD = kin.gradD();

    scalar maxMagEpsilon = 0;

    forAll(gradD, i)
    {
        maxMagEpsilon = max(maxMagEpsilon, mag(symm(gradD[i])));
    }

    return maxMagEpsilon;
}


Foam::wordList
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::
requiredScalarInputs() const
{
    if (!nonLocal_)
    {
        return wordList();
    }

    return wordList(1, nonLocalFieldName);
}


Foam::List<Foam::mechanicalConstitutiveLaw::implicitGradientField>
Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::
implicitGradientFields() const
{
    if (!nonLocal_)
    {
        return List<implicitGradientField>();
    }

    List<implicitGradientField> fields(1);
    fields[0].name = nonLocalFieldName;
    fields[0].source = "voidFraction";
    fields[0].length = characteristicLength_;

    return fields;
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    const UIndirectList<tensor>& gradD = kin.gradD();

    const Field<symmTensor>& epsilonP0 = state.getSymmTensorField0("epsilonP");
    Field<symmTensor>& epsilonP = state.symmTensorField("epsilonP");

    // The trial elastic strain is the strain less the plastic strain the
    // step started with
    Field<symmTensor> epsilonETrial(gradD.size());
    forAll(epsilonETrial, i)
    {
        epsilonETrial[i] = symm(gradD[i]) - epsilonP0[i];
    }

    Field<symmTensor> epsilonE(gradD.size());

    evaluateElasticPredictor
    (
        epsilonETrial, inputs, state, epsilonE, response
    );

    forAll(epsilonP, i)
    {
        epsilonP[i] = symm(gradD[i]) - epsilonE[i];
    }

    // No analytical fourth-order tangent has been derived for this law
    fourthOrderByFiniteDifferenceOnly(kin, inputs, state, response);
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::
evaluateElasticPredictor
(
    const UList<symmTensor>& epsilonETrial,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    UList<symmTensor>& epsilonE,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const scalar mu = mu_.value();
    const scalar kappa = kappa_.value();

    const Field<scalar>& epsilonPEq0 = state.getScalarField0("epsilonPEq");
    Field<scalar>& epsilonPEq = state.scalarField("epsilonPEq");
    Field<scalar>& sigmaY = state.scalarField("sigmaY");

    const Field<scalar>& voidFraction0 = state.getScalarField0("voidFraction");
    Field<scalar>& voidFraction = state.scalarField("voidFraction");

    const bool needScalarTan = response.wantsScalarTangent();
    const bool deviatoricTan =
        response.tangentReq() == tangentRequest::scalarDeviatoric;

    // Normalisation for the Newton residual: the scale the manager
    // supplies, or the largest trial elastic strain where it supplies none
    scalar maxMagDEpsilon = inputs.convergenceScale();

    if (maxMagDEpsilon <= 0)
    {
        forAll(epsilonETrial, i)
        {
            maxMagDEpsilon = max(maxMagDEpsilon, mag(epsilonETrial[i]));
        }
    }

    maxMagDEpsilon = max(maxMagDEpsilon, SMALL);

    const scalarField* fNonLocalPtr = nullptr;

    if (nonLocal_)
    {
        fNonLocalPtr = &inputs.getScalar(nonLocalFieldName);
    }

    pointResult result;

    forAll(sigma, i)
    {
        if (fNonLocalPtr)
        {
            nonLocalUpdate
            (
                epsilonETrial[i], epsilonPEq0[i], voidFraction0[i],
                (*fNonLocalPtr)[i], maxMagDEpsilon, result
            );
        }
        else
        {
            localUpdate
            (
                epsilonETrial[i], epsilonPEq0[i], voidFraction0[i],
                maxMagDEpsilon, result
            );
        }

        epsilonE[i] = result.epsilonE;
        sigma[i] = result.sigma;
        epsilonPEq[i] = epsilonPEq0[i] + result.DEpsilonPEq;
        sigmaY[i] = result.sigmaY;
        voidFraction[i] = result.f;

        if (needScalarTan)
        {
            // The elastic stiffness, scaled by the return in its deviatoric
            // part, which is stiffer than the true tangent and so aids
            // convergence only
            if (deviatoricTan)
            {
                response.scalarTangent()[i] = result.theta*(4.0/3.0)*mu;
            }
            else
            {
                response.scalarTangent()[i] =
                    result.theta*(4.0/3.0)*mu + kappa;
            }
        }
    }
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::endTimeStep
(
    const mechanicalConstitutiveLawState& state,
    const scalar time,
    const label timeIndex,
    DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics
) const
{
    appendMisesPlasticDiagnostics(state, diagnostics);

    const Field<scalar>& voidFraction = state.scalarField("voidFraction");

    scalar maxVoidFraction = 0;
    label nFailed = 0;

    forAll(voidFraction, i)
    {
        maxVoidFraction = max(maxVoidFraction, voidFraction[i]);

        if (voidFraction[i] >= fFailed_)
        {
            ++nFailed;
        }
    }

    typedef mechanicalConstitutiveLawDiagnostic diagnostic;

    diagnostics.append
    (
        diagnostic
        (
            "max void fraction",
            maxVoidFraction,
            diagnostic::combineOperation::maximum
        )
    );

    diagnostics.append
    (
        diagnostic
        (
            "failed integration points",
            scalar(nFailed),
            diagnostic::combineOperation::sum
        )
    );
}


void Foam::linearElasticGTNPlasticMechanicalConstitutiveLaw::reportDiagnostics
(
    const UList<mechanicalConstitutiveLawDiagnostic>& diagnostics
) const
{
    if (!debug || diagnostics.size() != 5)
    {
        return;
    }

    reportMisesPlasticDiagnostics
    (
        debug, SubList<mechanicalConstitutiveLawDiagnostic>(diagnostics, 3)
    );

    Info<< "Max void fraction = " << diagnostics[3].value()
        << ", failed integration points = " << label(diagnostics[4].value())
        << endl;
}


// ************************************************************************* //
