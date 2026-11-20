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

#include "linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw.H"
#include "misesPlasticDiagnostics.H"
#include "addToRunTimeSelectionTable.H"
#include "eig3.H"
#include "SubList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw, 1
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


const Foam::scalar
Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
Dfailed = 0.99;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar
Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
yieldStress
(
    const scalar r
) const
{
    return stressPlasticStrainSeries_(max(r, SMALL));
}


void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
returnMap
(
    const symmTensor& epsilonETrial,
    const scalar epsilonPEq0,
    const scalar D,
    const scalar maxMagDEpsilon,
    pointResult& result
) const
{
    const scalar mu = mu_.value();
    const scalar kappa = kappa_.value();
    const scalar omega = 1.0 - D;

    // Trial deviatoric effective stress
    const symmTensor sTrial(2.0*mu*dev(epsilonETrial));
    const scalar magSTrial = mag(sTrial);
    const scalar qTrial = sqrt(3.0/2.0)*magSTrial;

    const scalar sigmaY0 = yieldStress(epsilonPEq0);
    const scalar fTrial = qTrial - sigmaY0;

    scalar DLambda = 0;
    scalar sigmaY = sigmaY0;

    if (fTrial > SMALL)
    {
        // The effective equivalent stress returns by 3 mu dLambda/(1 - D):
        //     qTrial - 3 mu DLambda/(1 - D) - sigmaY(r0 + DLambda) = 0
        if (nonLinearPlasticity_)
        {
            scalar f = fTrial;
            label iter = 0;

            do
            {
                const scalar dSigmaY =
                    (
                        yieldStress(epsilonPEq0 + DLambda + finiteDiff_)
                      - yieldStress(epsilonPEq0 + DLambda)
                    )/finiteDiff_;

                const scalar dfdLambda = -3.0*mu/omega - dSigmaY;

                const scalar dLambda = -f/dfdLambda;
                DLambda += dLambda;

                f =
                    qTrial - 3.0*mu*DLambda/omega
                  - yieldStress(epsilonPEq0 + DLambda);

                if (mag(dLambda)/maxMagDEpsilon < loopTol_)
                {
                    break;
                }
            } while (++iter < maxNewtonIter_);

            if (iter >= maxNewtonIter_)
            {
                WarningInFunction
                    << "Plasticity Newton loop not converging" << endl;
            }
        }
        else
        {
            DLambda = fTrial/(3.0*mu/omega + Hp_);
        }

        sigmaY = yieldStress(epsilonPEq0 + DLambda);
    }

    // Flow direction and plastic strain increment
    symmTensor n(symmTensor::zero);
    if (magSTrial > SMALL)
    {
        n = sTrial/magSTrial;
    }

    result.DEpsilonP = (sqrt(3.0/2.0)*DLambda/omega)*n;
    result.epsilonE = epsilonETrial - result.DEpsilonP;
    result.sigma =
        omega
       *(
            2.0*mu*dev(result.epsilonE)
          + kappa*tr(result.epsilonE)*symmTensor(I)
        );
    result.DLambda = DLambda;
    result.D = D;
    result.sigmaY = sigmaY;

    result.theta = 1.0;
    if (magSTrial > SMALL)
    {
        result.theta =
            max(1.0 - 2.0*mu*mag(result.DEpsilonP)/magSTrial, 0.0);
    }
}


Foam::scalar
Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
energyReleaseRate
(
    const symmTensor& sigma,
    const scalar D
) const
{
    const scalar E = E_.value();
    const scalar nu = nu_.value();
    const scalar trSigma = tr(sigma);

    if (h_ >= 1.0 - SMALL)
    {
        return
           -((1.0 + nu)*(sigma && sigma) - nu*sqr(trSigma))
           /(2.0*E*sqr(1.0 - D));
    }

    // The tensile and compressive parts, from the principal stresses. The
    // eigenvectors are the rows of V
    tensor V(tensor::zero);
    vector s(vector::zero);
    eig3::eigen_decomposition(sigma, V, s);

    symmTensor sigmaPos(symmTensor::zero);
    symmTensor sigmaNeg(symmTensor::zero);

    for (direction i = 0; i < 3; ++i)
    {
        const vector ei(i == 0 ? V.x() : (i == 1 ? V.y() : V.z()));

        if (s[i] > 0)
        {
            sigmaPos += s[i]*sqr(ei);
        }
        else
        {
            sigmaNeg += s[i]*sqr(ei);
        }
    }

    const scalar YPos =
        ((1.0 + nu)*(sigmaPos && sigmaPos) - nu*sqr(max(trSigma, 0.0)))
       /(2.0*E*sqr(1.0 - D));

    const scalar YNeg =
        h_*((1.0 + nu)*(sigmaNeg && sigmaNeg) - nu*sqr(max(-trSigma, 0.0)))
       /(2.0*E*sqr(1.0 - h_*D));

    return -YPos - YNeg;
}


Foam::scalar
Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
damageIncrement
(
    const pointResult& result,
    const scalar epsilonPEq0
) const
{
    // No plastic flow, or not yet past the damage threshold
    if
    (
        result.DLambda <= 0
     || epsilonPEq0 + result.DLambda <= epsilonD_
    )
    {
        return 0;
    }

    const symmTensor& sigma = result.sigma;
    const symmTensor s(dev(sigma));
    const scalar sigmaEq = sqrt(3.0/2.0)*mag(s);

    if (sigmaEq < SMALL)
    {
        return 0;
    }

    // Stress triaxiality, the same for the effective stress as for the
    // stress, since (1 - D) scales both parts alike
    const scalar eta = tr(sigma)/(3.0*sigmaEq);

    if (triaxialityCutOff_ && eta <= -1.0/3.0)
    {
        return 0;
    }

    scalar S = S0_;

    if (lemaitreType_ == lemaitreType::Malcher)
    {
        // Normalised third invariant, clipped against round-off
        const scalar xi =
            max(min((27.0/2.0)*det(s)/pow3(sigmaEq), 1.0), -1.0);

        S =
            S0Tension_
           /(3.0*mag(eta) + (S0Tension_/S0Shear_)*(1.0 - sqr(xi)));
    }
    else if (lemaitreType_ == lemaitreType::triaxiality)
    {
        // Positive, since the cut-off has removed eta <= -1/3
        S = 2.0*S0_/(1.0 + 3.0*eta);
    }

    const scalar Y = energyReleaseRate(sigma, result.D);

    return
        result.DLambda/(1.0 - result.D)*pow(max(-Y, 0.0)/S, b_);
}


void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
localUpdate
(
    const symmTensor& epsilonETrial,
    const scalar epsilonPEq0,
    const scalar D0,
    const scalar maxMagDEpsilon,
    pointResult& result
) const
{
    // A failed point stays failed
    if (D0 >= Dc_)
    {
        returnMap(epsilonETrial, epsilonPEq0, Dfailed, maxMagDEpsilon, result);
        return;
    }

    // At the damage the step started with. No damage increment there means
    // none in the bracket: an elastic step stays elastic at any damage, and
    // a greater damage only reduces the plastic multiplier, so a step below
    // the threshold stays below it
    returnMap(epsilonETrial, epsilonPEq0, D0, maxMagDEpsilon, result);

    const scalar DD0 = damageIncrement(result, epsilonPEq0);

    if (DD0 <= 0)
    {
        return;
    }

    // The residual D - D0 - dD(D) is negative at D0. If it is still negative
    // at Dc, the step takes the damage past Dc and the point fails
    pointResult resultC;
    returnMap(epsilonETrial, epsilonPEq0, Dc_, maxMagDEpsilon, resultC);

    scalar Ra = -DD0;
    scalar Rb = Dc_ - D0 - damageIncrement(resultC, epsilonPEq0);

    if (Rb < 0)
    {
        returnMap(epsilonETrial, epsilonPEq0, Dfailed, maxMagDEpsilon, result);
        return;
    }

    // Root in [D0, Dc], by regula falsi with the Illinois modification
    scalar a = D0;
    scalar b = Dc_;
    label side = 0;

    for (label iter = 0; iter < maxDamageIter_; ++iter)
    {
        const scalar D = (a*Rb - b*Ra)/(Rb - Ra);

        returnMap(epsilonETrial, epsilonPEq0, D, maxMagDEpsilon, result);

        const scalar R = D - D0 - damageIncrement(result, epsilonPEq0);

        if (mag(R) < damageTol_ || (b - a) < damageTol_)
        {
            return;
        }

        if (R < 0)
        {
            a = D;
            Ra = R;

            if (side == -1)
            {
                Rb /= 2;
            }
            side = -1;
        }
        else
        {
            b = D;
            Rb = R;

            if (side == 1)
            {
                Ra /= 2;
            }
            side = 1;
        }
    }

    WarningInFunction
        << "Coupled damage update not converging" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw
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
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    Hp_(0),
    lemaitreType_(lemaitreType::classic),
    S0_(0),
    S0Tension_(0),
    S0Shear_(0),
    b_(dict.lookupOrDefault<scalar>("b", 1.0)),
    epsilonD_(dict.lookupOrDefault<scalar>("epsilonD", 0.0)),
    Dc_(dict.lookupOrDefault<scalar>("Dc", Dfailed)),
    h_(dict.lookupOrDefault<scalar>("crackClosure", 1.0)),
    triaxialityCutOff_
    (
        dict.lookupOrDefault<Switch>("triaxialityCutOff", true)
    ),
    loopTol_(dict.lookupOrDefault<scalar>("NewtonLoopTol", 1e-8)),
    maxNewtonIter_(dict.lookupOrDefault<label>("NewtonMaxIter", 200)),
    finiteDiff_(dict.lookupOrDefault<scalar>("NewtonFiniteDiffEps", 0.25e-6)),
    damageTol_(dict.lookupOrDefault<scalar>("damageTol", 1e-12)),
    maxDamageIter_(dict.lookupOrDefault<label>("damageMaxIter", 100))
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

    if (stressPlasticStrainSeries_.size() == 2)
    {
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

    // The damage denominator
    const word type(dict.lookupOrDefault<word>("lemaitreType", "classic"));

    if (type == "classic" || type == "triaxiality")
    {
        lemaitreType_ =
            type == "classic"
          ? lemaitreType::classic
          : lemaitreType::triaxiality;

        const dimensionedScalar S0(dict.lookup("S0"));

        if (S0.dimensions() != dimPressure || S0.value() <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "S0 must be positive, with dimensions " << dimPressure
                << exit(FatalIOError);
        }

        S0_ = S0.value();
    }
    else if (type == "Malcher")
    {
        lemaitreType_ = lemaitreType::Malcher;

        const dimensionedScalar S0Tension(dict.lookup("S0Tension"));
        const dimensionedScalar S0Shear(dict.lookup("S0Shear"));

        if
        (
            S0Tension.dimensions() != dimPressure
         || S0Shear.dimensions() != dimPressure
         || S0Tension.value() <= 0
         || S0Shear.value() <= 0
        )
        {
            FatalIOErrorInFunction(dict)
                << "S0Tension and S0Shear must be positive, with dimensions "
                << dimPressure << exit(FatalIOError);
        }

        S0Tension_ = S0Tension.value();
        S0Shear_ = S0Shear.value();
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown lemaitreType " << type << ". Valid types are "
            << "classic, Malcher and triaxiality" << exit(FatalIOError);
    }

    if (lemaitreType_ == lemaitreType::triaxiality && !triaxialityCutOff_)
    {
        FatalIOErrorInFunction(dict)
            << "lemaitreType triaxiality needs the triaxiality cut-off: its "
            << "denominator 2*S0/(1 + 3*eta) is singular at eta = -1/3"
            << exit(FatalIOError);
    }

    if (b_ <= 0 || epsilonD_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "b must be positive and epsilonD must not be negative"
            << exit(FatalIOError);
    }

    if (Dc_ <= 0 || Dc_ > Dfailed)
    {
        FatalIOErrorInFunction(dict)
            << "The critical damage Dc = " << Dc_ << " must lie in (0, "
            << Dfailed << "]" << exit(FatalIOError);
    }

    if (h_ < 0 || h_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "The crack closure parameter " << h_ << " must lie in [0, 1]"
            << exit(FatalIOError);
    }

    Info<< "    Lemaitre damage: " << type << ", b = " << b_
        << ", epsilonD = " << epsilonD_ << ", Dc = " << Dc_
        << ", crackClosure = " << h_
        << ", triaxialityCutOff = " << Switch(triaxialityCutOff_) << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    typedef mechanicalConstitutiveLawStateSpec::stateRole role;

    spec.addSymmTensor("epsilonP", role::persistent, symmTensor::zero);

    // The accumulated plastic multiplier, which drives the hardening
    spec.addScalar("epsilonPEq", role::persistent, 0.0);

    spec.addScalar("sigmaY", role::persistent, stressPlasticStrainSeries_(0.0));

    spec.addScalar("damage", role::persistent, 0.0);
}


Foam::scalar
Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
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


void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
evaluate
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


void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
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

    const Field<scalar>& damage0 = state.getScalarField0("damage");
    Field<scalar>& damage = state.scalarField("damage");

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

    pointResult result;

    forAll(sigma, i)
    {
        localUpdate
        (
            epsilonETrial[i], epsilonPEq0[i], damage0[i], maxMagDEpsilon,
            result
        );

        epsilonE[i] = result.epsilonE;
        sigma[i] = result.sigma;
        epsilonPEq[i] = epsilonPEq0[i] + result.DLambda;
        sigmaY[i] = result.sigmaY;
        damage[i] = result.D;

        if (needScalarTan)
        {
            // The elastic stiffness, degraded by the damage and scaled by
            // the return
            const scalar omega = 1.0 - result.D;

            if (deviatoricTan)
            {
                response.scalarTangent()[i] =
                    omega*result.theta*(4.0/3.0)*mu;
            }
            else
            {
                response.scalarTangent()[i] =
                    omega*(result.theta*(4.0/3.0)*mu + kappa);
            }
        }
    }
}


void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
endTimeStep
(
    const mechanicalConstitutiveLawState& state,
    const scalar time,
    const label timeIndex,
    DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics
) const
{
    appendMisesPlasticDiagnostics(state, diagnostics);

    const Field<scalar>& damage = state.scalarField("damage");

    scalar maxDamage = 0;
    label nFailed = 0;

    forAll(damage, i)
    {
        maxDamage = max(maxDamage, damage[i]);

        if (damage[i] >= Dc_)
        {
            ++nFailed;
        }
    }

    typedef mechanicalConstitutiveLawDiagnostic diagnostic;

    diagnostics.append
    (
        diagnostic
        (
            "max damage", maxDamage, diagnostic::combineOperation::maximum
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


void Foam::linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw::
reportDiagnostics
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

    Info<< "Max damage = " << diagnostics[3].value()
        << ", failed integration points = " << label(diagnostics[4].value())
        << endl;
}


// ************************************************************************* //
