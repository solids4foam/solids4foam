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
    const scalar sigmaY =
        yieldStress(epsilonPEq0 + sqrt(2.0/3.0)*DLambda);

    return qTrial - 3.0*muBar*sqrt(3.0/2.0)*DLambda - sigmaY;
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

    curSigmaY =
        yieldStress(epsilonPEq0 + sqrt(2.0/3.0)*DLambda);
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
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_("mu", E_.dimensions(), 0.0),
    kappa_("kappa", E_.dimensions(), 0.0),
    stressPlasticStrainSeries_(dict),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    Hp_(0),
    loopTol_(dict.lookupOrDefault<scalar>("NewtonLoopTol", 1e-8)),
    maxNewtonIter_(dict.lookupOrDefault<label>("NewtonMaxIter", 200)),
    finiteDiff_(dict.lookupOrDefault<scalar>("NewtonFiniteDiffEps", 0.25e-6))
{
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

void Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::initialiseState
(
    mechanicalConstitutiveLawState& state
) const
{
    // Total plastic strain
    state.symmTensorField("epsilonP") = symmTensor::zero;

    // Total plastic strain (old time)
    state.symmTensorField0("epsilonP") = symmTensor::zero;

    // Equivalent plastic strain
    state.scalarField("epsilonPEq") = 0.0;

    // Equivalent plastic strain (old time)
    state.scalarField0("epsilonPEq") = 0.0;

    // Yield stress
    state.scalarField("sigmaY") = stressPlasticStrainSeries_(0.0);
}


void Foam::linearElasticMisesPlasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
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
            //     q = qTrial - 3 μ sqrt(3/2) Δλ
            //
            // The equivalent plastic strain evolves as:
            //
            //     εp_eq = εp_eq0 + sqrt(2/3) Δλ
            //
            // For linear isotropic hardening:
            //
            //     σY = σY0 + Hp * (sqrt(2/3) Δλ)
            //
            // Substituting into the consistency condition f = q − σY = 0:
            //
            //     qTrial
            //   − 3 μ sqrt(3/2) Δλ
            //   − (σY0 + Hp sqrt(2/3) Δλ) = 0
            //
            // which gives the closed-form solution:
            //
            //     Δλ = (qTrial − σY0)
            //           / (3 μ sqrt(3/2) + Hp sqrt(2/3))
            //
            // Perfect plasticity is recovered by setting Hp = 0.
            const scalar denom = 3.0*mu*sqrt(3.0/2.0) + Hp_*sqrt(2.0/3.0);
            dLambda = fTrial/max(denom, SMALL);
            curSigmaY = sigmaYTrial + Hp_*sqrt(2.0/3.0)*dLambda;
        }

        // Plastic strain increment = (3/2)*dLambda*n  (standard J2 small-strain)
        const symmTensor DEpsilonP = 1.5*dLambda*n;

        epsilonP[i] = epsilonP0[i] + DEpsilonP;
        epsilonPEq[i] = epsilonPEq0[i] + sqrt(2.0/3.0)*dLambda;
        sigmaY[i] = curSigmaY;

        // Returned deviatoric stress:
        // s = sTrial - 2 mu DEpsilonP
        // with DEpsilonP = 3/2 dLambda n => s = sTrial - 3 mu dLambda n
        const symmTensor s = sTrial - (3.0*mu*dLambda)*n;

        // Total small strain stress
        sigma[i] = s + kappa*trEps*I;

        if (needScalarTan)
        {
            // Approximate scalar tangent
            scalar theta = 1.0;
            if (magSTrial > SMALL)
            {
                theta = 1.0 - (3.0*mu*dLambda)/magSTrial;
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
        finiteDifferenceFourthOrder(kin, state, response);
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
