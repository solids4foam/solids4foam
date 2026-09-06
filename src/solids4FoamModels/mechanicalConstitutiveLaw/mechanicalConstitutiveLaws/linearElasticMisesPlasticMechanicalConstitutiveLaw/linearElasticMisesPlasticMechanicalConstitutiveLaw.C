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

#include "linearElasticMisesPlasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        linearElasticMisesPlasticMechanicalConstitutiveLaw, 1
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        linearElasticMisesPlasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::yieldStress
(
    const scalar epsilonPEq
) const
{
    return stressPlasticStrainSeries_(max(epsilonPEq, SMALL));
}

inline Foam::scalar
Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::yieldFunction
(
    const scalar epsilonPEq0,
    const scalar qTrial,
    const scalar DLambda,
    const scalar muBar
) const
{
    // q = sqrt(3/2) |s|
    // For radial return: q = qTrial - 3*muBar*DLambda
    // with epsilonPEq = epsilonPEq0 + sqrt(2/3)*DLambda
    // DLambda is the equivalent plastic strain increment itself, so the
    // hardening curve is evaluated at epsilonPEq0 + DLambda and the return
    // mapping reduces q by 3*mu*DLambda. See the derivation at the closed-form
    // branch in evaluate()
    const scalar sigmaY = yieldStress(epsilonPEq0 + DLambda);

    return qTrial - 3.0*muBar*DLambda - sigmaY;
}


inline void Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::newtonLoop
(
    scalar& DLambda,
    scalar& curSigmaY,
    const scalar epsilonPEq0,
    const scalar qTrial,
    const scalar muBar,
    const scalar maxMagDEpsilon
) const
{
    label iter = 0;

    scalar f = yieldFunction(epsilonPEq0, qTrial, DLambda, muBar);
    scalar residual = 1.0;

    do
    {
        const scalar fStep =
            yieldFunction(epsilonPEq0, qTrial, DLambda + finiteDiff_, muBar);

        const scalar dfdLambda = (fStep - f)/finiteDiff_;

        if (mag(dfdLambda) < SMALL)
        {
            FatalErrorInFunction
                << "df/dLambda ~ 0 in plasticity Newton loop"
                << exit(FatalError);
        }

        const scalar dLambda = f/dfdLambda;
        DLambda -= dLambda;

        residual = dLambda/(maxMagDEpsilon + SMALL);

        f = yieldFunction(epsilonPEq0, qTrial, DLambda, muBar);

    } while (mag(residual) > loopTol_ && ++iter < maxNewtonIter_);

    if (iter >= maxNewtonIter_)
    {
        WarningInFunction
            << "Plasticity Newton loop not converging" << nl;
    }

    // At the same equivalent plastic strain the yield function was solved
    // with, and that the strain update below uses. This law's DLambda is the
    // equivalent plastic strain increment itself, where the legacy law's is
    // the plastic multiplier and increments it by sqrt(2/3)*DLambda; this line
    // kept the legacy form and so applied that factor a second time, storing a
    // yield stress from a plastic strain the solve never visited
    curSigmaY = yieldStress(epsilonPEq0 + DLambda);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::
linearElasticMisesPlasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    mu_("mu", E_.dimensions(), 0.0),
    kappa_("kappa", E_.dimensions(), 0.0),
    stressPlasticStrainSeries_(dict),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    Hp_(0),
    loopTol_(dict.lookupOrDefault<scalar>("NewtonLoopTol", 1e-8)),
    maxNewtonIter_(dict.lookupOrDefault<label>("NewtonMaxIter", 200)),
    finiteDiff_(dict.lookupOrDefault<scalar>("NewtonFiniteDiffEps", 0.25e-6))
{
    // The material may be given either as E and nu or as mu and K, matching
    // the legacy law, so that an existing case dictionary needs no change.
    // The legacy law tries E and nu first, so this does too
    if (dict.found("E") && dict.found("nu"))
    {
        E_ = dimensionedScalar(dict.lookup("E"));
        nu_ = dimensionedScalar(dict.lookup("nu"));
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        const dimensionedScalar muIn(dict.lookup("mu"));
        const dimensionedScalar KIn(dict.lookup("K"));

        if (muIn.dimensions() != dimPressure || KIn.dimensions() != dimPressure)
        {
            FatalIOErrorInFunction(dict)
                << "The shear modulus mu and bulk modulus K must both have "
                << "dimensions " << dimPressure
                << exit(FatalIOError);
        }

        // Invert mu = E/(2*(1 + nu)) and K = E/(3*(1 - 2*nu))
        E_ = 9.0*KIn*muIn/(3.0*KIn + muIn);
        nu_ = (3.0*KIn - 2.0*muIn)/(2.0*(3.0*KIn + muIn));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Specify the elastic properties either as 'E' and 'nu' or as "
            << "'mu' and 'K'."
            << exit(FatalIOError);
    }

    if (E_.dimensions() != dimPressure)
    {
        FatalIOErrorInFunction(dict)
            << "Young's modulus E has incorrect dimensions. "
            << "Expected " << dimPressure << " but got "
            << E_.dimensions()
            << exit(FatalIOError);
    }

    if (nu_.dimensions() != dimless)
    {
        FatalIOErrorInFunction(dict)
            << "Poisson's ratio nu must be dimensionless. "
            << "Got " << nu_.dimensions()
            << exit(FatalIOError);
    }

    if (nu_.value() >= 0.5 - SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Poisson's ratio nu = " << nu_.value()
            << " is too close to 0.5 for linear elasticity. "
            << "This leads to an ill-conditioned bulk modulus."
            << exit(FatalIOError);
    }
    else if (nu_.value() <= -1.0)
    {
        FatalIOErrorInFunction(dict)
            << "Invalid Poisson's ratio nu = " << nu_.value()
            << ". Expected -1 <= nu for linear elasticity."
            << exit(FatalIOError);
    }

    // Plane stress is not supported, matching the legacy linearElasticMisesPlastic
    // Note: the planeStress entry is injected into this dictionary by the
    // mechanicalConstitutiveLawManager from the top-level entry in
    // mechanicalProperties; it is not given by the user in this sub-dictionary
    if (Switch(dict.lookupOrDefault<Switch>("planeStress", false)))
    {
        FatalIOErrorInFunction(dict)
            << "Not implemented for planeStress. If needed, you can solve the "
            << "case in 3-D and set the back to a symmetryPlane and the front "
            << "to traction-free" << exit(FatalIOError);
    }

    // Set mu and kappa
    mu_ = E_/(2.0*(1.0 + nu_));
    kappa_ = E_/(3.0*(1 - 2*nu_));

    // Set Hp for linear hardening
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
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // This law's history, declared rather than created in initialiseState so
    // that the manager knows it exists. Nothing can be written for a restart
    // that has not been declared, and these are exactly the fields a restart
    // has to carry: a plastic law resumed without them starts again from an
    // unyielded state.
    //
    // The defaults are the values initialiseState used to assign, so a cold
    // start is unchanged
    spec.addSymmTensor
    (
        "epsilonP",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        symmTensor::zero
    );

    spec.addScalar
    (
        "epsilonPEq",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        0.0
    );

    spec.addScalar
    (
        "sigmaY",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        stressPlasticStrainSeries_(0.0)
    );
}


void Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    // Output
    UIndirectList<symmTensor>& sigma = response.stress();

    // Kinematics
    const UIndirectList<tensor>& gradD = kin.gradD();

    // Material constants
    const scalar mu = mu_.value();
    const scalar kappa = kappa_.value();

    // State (current + old-time)
    Field<symmTensor>& epsilonP = state.symmTensorField("epsilonP");
    const Field<symmTensor>& epsilonP0 =
        state.getSymmTensorField0("epsilonP");

    Field<scalar>& epsilonPEq = state.scalarField("epsilonPEq");
    const Field<scalar>& epsilonPEq0 =
        state.getScalarField0("epsilonPEq");

    // Optional yield stress state (handy for debugging/post-proc)
    Field<scalar>& sigmaY = state.scalarField("sigmaY");

    // Optional scalar tangent
    // Whether the caller wants the deviatoric stress and the volumetric
    // response separately, as a mixed displacement-pressure formulation does.
    // The return mapping only ever touches the deviatoric trial stress, and
    // J2 plasticity has no plastic volume change, so the volumetric term is
    // already a separate quantity here rather than something to extract
    const bool wantsSplit = response.wantsVolumetricSplit();
    UIndirectList<scalar>* volumetricPtr =
        wantsSplit ? &response.volumetric() : nullptr;

    const bool needScalarTan = response.wantsScalarTangent();
    UIndirectList<scalar>* scalarTanPtr = nullptr;
    if (needScalarTan)
    {
        scalarTanPtr = &response.scalarTangent();
    }

    // Normalisation for Newton residual
    // Use max equivalent strain (small strain) as a scale
    scalar maxMagDEpsilon = SMALL;
    forAll(sigma, i)
    {
        const symmTensor eps = symm(gradD[i]);
        maxMagDEpsilon = max(maxMagDEpsilon, mag(eps));
    }

    // Loop over integration points
    forAll(sigma, i)
    {
        const symmTensor eps = symm(gradD[i]);
        const scalar trEps = tr(eps);

        const symmTensor eDev = dev(eps);

        // Elastic predictor (trial deviatoric stress)
        const symmTensor sTrial = 2.0*mu*(eDev - dev(epsilonP0[i]));

        const scalar magSTrial = mag(sTrial);

        // qTrial = sqrt(3/2) |s_trial|
        const scalar qTrial = sqrt(3.0/2.0)*magSTrial;

        // Yield stress from old eq plastic strain
        // (we could store this for convenience)
        const scalar sigmaYTrial = yieldStress(epsilonPEq0[i]);

        // fTrial = qTrial - sigmaY
        const scalar fTrial = qTrial - sigmaYTrial;

        if (fTrial <= SMALL)
        {
            // Elastic
            epsilonP[i] = epsilonP0[i];
            epsilonPEq[i] = epsilonPEq0[i];
            sigmaY[i] = sigmaYTrial;

            // Full stress = hydro + deviatoric
            sigma[i] = sTrial + kappa*trEps*I;

            if (wantsSplit)
            {
                (*volumetricPtr)[i] = kappa*trEps;
                sigma[i] = sTrial;
            }

            if (needScalarTan)
            {
                // Small-strain scalar tangent options
                // - scalar: bulk + deviatoric contribution (usual "2mu+lambda"
                //   style)
                // - scalarDeviatoric: deviatoric-only (useful for mixed u-p
                //   solve)
                // Here, represent them in (kappa, mu) form:
                // lambda = kappa - 2/3 mu
                // 2mu + lambda = 4/3 mu + kappa
                if (response.tangentReq() == tangentRequest::scalarDeviatoric)
                {
                    // Scalar Laplacian surrogate for div(dev(sigma)), which is
                    // mu*lap(D) + (1/3)*mu*grad(div(D)).
                    // This is the theta = 1 limit of the plastic branch below
                    (*scalarTanPtr)[i] = (4.0/3.0)*mu;
                }
                else
                {
                    (*scalarTanPtr)[i] = (4.0/3.0)*mu + kappa;
                }
            }

            continue;
        }

        // Plastic step: compute return direction
        symmTensor n = symmTensor::zero;
        if (magSTrial > SMALL)
        {
            n = sTrial/magSTrial;
        }
        else
        {
            // Very unlikely if fTrial>0, but be defensive
            n = symmTensor::zero;
        }

        // Solve for DLambda and updated sigmaY
        scalar dLambda = 0.0;
        scalar curSigmaY = sigmaYTrial;

        if (nonLinearPlasticity_)
        {
            // Newton loop (table-based nonlinear hardening)
            newtonLoop
            (
                dLambda,
                curSigmaY,
                epsilonPEq0[i],
                qTrial,
                mu,
                maxMagDEpsilon
            );
        }
        else
        {
            // Linear / perfect plasticity (closed-form radial return)
            //
            // Yield function is written in terms of the equivalent (von Mises)
            // stress:
            //
            //     q = sqrt(3/2) * |s|
            //
            // For associative J2 plasticity with plastic multiplier increment
            // Δλ, the radial return mapping gives:
            //
            //     q = qTrial - 3 μ Δλ
            //
            // The equivalent plastic strain evolves as:
            //
            //     εp_eq = εp_eq0 + Δλ
            //
            // For linear isotropic hardening:
            //
            //     σY = σY0 + Hp Δλ
            //
            // Substituting into the consistency condition f = q − σY = 0:
            //
            //     qTrial
            //   − 3 μ Δλ
            //   − (σY0 + Hp Δλ) = 0
            //
            // which gives the closed-form solution:
            //
            //     Δλ = (qTrial − σY0) / (3 μ + Hp)
            //
            // Perfect plasticity is recovered by setting Hp = 0.
            const scalar denom = 3.0*mu + Hp_;
            dLambda = fTrial/max(denom, SMALL);
            curSigmaY = sigmaYTrial + Hp_*dLambda;
        }

        // Plastic flow, with dLambda the equivalent plastic strain increment:
        //
        //     DEpsilonP  = (3/2) dLambda s/q = sqrt(3/2) dLambda n
        //     |DEpsilonP| = sqrt(3/2) dLambda
        //     DEpsilonPEq = sqrt(2/3) |DEpsilonP| = dLambda
        //
        // which is self-consistent, unlike the previous form that combined
        // DEpsilonP = (3/2) dLambda n with DEpsilonPEq = sqrt(2/3) dLambda and
        // so disagreed with itself by a factor of 3/2
        const symmTensor DEpsilonP = sqrt(3.0/2.0)*dLambda*n;

        epsilonP[i] = epsilonP0[i] + DEpsilonP;
        epsilonPEq[i] = epsilonPEq0[i] + dLambda;
        sigmaY[i] = curSigmaY;

        // Returned deviatoric stress:
        //     s = sTrial - 2 mu DEpsilonP = sTrial - 2 mu sqrt(3/2) dLambda n
        // so that |s| = |sTrial| - 2 mu sqrt(3/2) dLambda and hence
        //     q = qTrial - 3 mu dLambda
        const symmTensor s = sTrial - (2.0*mu*sqrt(3.0/2.0)*dLambda)*n;

        // Total small strain stress
        sigma[i] = s + kappa*trEps*I;

        if (wantsSplit)
        {
            (*volumetricPtr)[i] = kappa*trEps;
            sigma[i] = s;
        }

        if (needScalarTan)
        {
            // Approximate scalar tangent
            scalar theta = 1.0;
            if (magSTrial > SMALL)
            {
                theta = 1.0 - (2.0*mu*sqrt(3.0/2.0)*dLambda)/magSTrial;
                theta = max(theta, 0.0);
            }

            if (response.tangentReq() == tangentRequest::scalarDeviatoric)
            {
                (*scalarTanPtr)[i] = theta*(4.0/3.0)*mu;
            }
            else
            {
                (*scalarTanPtr)[i] = theta*(4.0/3.0)*mu + kappa;
            }
        }
    }

    // Fourth-order tangent.
    // There is no analytical consistent tangent for this law yet, but the
    // finite-difference tangent of the base class is well defined for any law
    // and is evaluated against a shadow state, so it neither disturbs the
    // return mapping just performed nor the history it started from
    if (response.tangentReq() == tangentRequest::fourthOrderFiniteDifference)
    {
        finiteDifferenceFourthOrder(kin, inputs, state, response);
    }
    else if (response.tangentReq() == tangentRequest::fourthOrder)
    {
        FatalErrorInFunction
            << "An analytical fourth-order tangent is not implemented for "
            << type() << "." << nl
            << "Use 'fourthOrderFiniteDifference' to obtain one by finite "
            << "differences." << exit(FatalError);
    }
}


void Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::endTimeStep
(
    const mechanicalConstitutiveLawState& state,
    const scalar time,
    const label timeIndex
) const
{
    const Field<scalar>& epsilonPEq = state.scalarField("epsilonPEq");
    const Field<scalar>& epsilonPEq0 =
        state.getScalarField0("epsilonPEq");

    label nYielding = 0;
    scalar curDEpsilonPEq = 0.0;
    scalar maxDEpsilonPEq = 0.0;
    forAll(epsilonPEq, i)
    {
        curDEpsilonPEq = epsilonPEq[i] - epsilonPEq0[i];

        if (curDEpsilonPEq > SMALL)
        {
            ++nYielding;

            maxDEpsilonPEq = max(maxDEpsilonPEq, curDEpsilonPEq);
        }
    }

    reduce(nYielding, sumOp<label>());
    reduce(maxDEpsilonPEq, maxOp<scalar>());

    const int nTotalCells = returnReduce(epsilonPEq.size(), sumOp<int>());

    if (debug)
    {
        Info<< nl << "Max DEpsilonPEq is " << maxDEpsilonPEq << nl
            << "Number of yielding integration points = " << nYielding
            << "/" << nTotalCells
            << nl << endl;
    }
}


// ************************************************************************* //
