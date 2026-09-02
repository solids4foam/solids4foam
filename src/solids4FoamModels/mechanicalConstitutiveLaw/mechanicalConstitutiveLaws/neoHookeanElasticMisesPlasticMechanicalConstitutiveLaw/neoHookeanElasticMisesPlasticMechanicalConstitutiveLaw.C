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

#include "neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw, 1
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * * //

inline Foam::scalar
Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::yieldStress
(
    const scalar epsilonPEq,
    const scalar J
) const
{
    // The curve is Cauchy stress against true strain, and the return mapping
    // is written in Kirchhoff stress, so scale by J
    return J*stressPlasticStrainSeries_(max(epsilonPEq, SMALL));
}


inline Foam::scalar
Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::yieldFunction
(
    const scalar epsilonPEqOld,
    const scalar magSTrial,
    const scalar dLambda,
    const scalar muBar,
    const scalar J
) const
{
    // f = |s| - sqrt(2/3) sigmaY, with |s| = |sTrial| - 2 muBar dLambda,
    // and the yield stress evaluated at the updated plastic strain
    return
        magSTrial - 2.0*muBar*dLambda
      - sqrt(2.0/3.0)
       *yieldStress(epsilonPEqOld + sqrt(2.0/3.0)*dLambda, J);
}


inline void
Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::newtonLoop
(
    scalar& dLambda,
    scalar& curSigmaY,
    const scalar epsilonPEqOld,
    const scalar magSTrial,
    const scalar muBar,
    const scalar J,
    const scalar maxMagDEpsilon
) const
{
    label i = 0;
    scalar f = yieldFunction(epsilonPEqOld, magSTrial, dLambda, muBar, J);
    scalar residual = 1.0;

    do
    {
        // First-order numerical derivative: two function evaluations per
        // iteration, following Hauser (2009)
        const scalar fStep =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, dLambda + finiteDiff_, muBar, J
            );

        const scalar fDeriv = (fStep - f)/finiteDiff_;

        residual = f/fDeriv;
        dLambda -= residual;

        // Normalise with respect to the largest strain increment, so the
        // tolerance means the same thing at every integration point
        residual /= maxMagDEpsilon;

        f = yieldFunction(epsilonPEqOld, magSTrial, dLambda, muBar, J);

        if (i == maxNewtonIter_)
        {
            WarningInFunction
                << "Plasticity Newton loop is not converging" << endl;
        }
    }
    while ((mag(residual) > loopTol_) && ++i < maxNewtonIter_);

    // Back to a Cauchy yield stress
    curSigmaY =
        yieldStress(epsilonPEqOld + sqrt(2.0/3.0)*dLambda, J)/J;
}


inline Foam::scalar
Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::IbarConsistent
(
    const symmTensor& devBEbar
) const
{
    // Simo & Hughes take the trial value for the spherical part, which does
    // not give det(bEbar) == 1. Rubin and Attia (1996), in the form
    // recommended by Hollenstein et al. (2013), obtain it from a cubic.
    // Translated from Rubin's SmoothMultiPhase subroutine
    const scalar detdevBepr = det(devBEbar);
    const scalar dotprod = devBEbar && devBEbar;
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

    return alpha1/3.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::
neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw
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
    finiteDiff_(dict.lookupOrDefault<scalar>("NewtonFiniteDiffEps", 0.25e-6)),
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>("updateBEbarConsistent", true)
    )
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

    // Plane stress is not supported, matching the legacy neoHookeanElasticMisesPlastic
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

void Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::initialiseState
(
    mechanicalConstitutiveLawState& state
) const
{
    // Isochoric elastic left Cauchy-Green tensor. The undeformed value is the
    // identity, not zero, so both times must be set explicitly
    state.symmTensorField("bEbar") = symmTensor(I);
    state.symmTensorField0("bEbar") = symmTensor(I);

    // Equivalent plastic strain
    state.scalarField("epsilonPEq") = 0.0;
    state.scalarField0("epsilonPEq") = 0.0;

    // Yield stress, kept for post-processing
    state.scalarField("sigmaY") = stressPlasticStrainSeries_(0.0);
}


void Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    // Output
    UIndirectList<symmTensor>& sigma = response.stress();

    // Kinematics
    const UIndirectList<tensor>& F = kin.F();
    const UIndirectList<tensor>& Finv0 = kin.Finv0();
    const UIndirectList<scalar>& J = kin.J();
    const UIndirectList<scalar>& J0 = kin.J0();

    // Material constants
    const scalar mu = mu_.value();
    const scalar kappa = kappa_.value();

    // State: bEbar is the history this formulation turns on
    Field<symmTensor>& bEbar = state.symmTensorField("bEbar");
    const Field<symmTensor>& bEbar0 = state.getSymmTensorField0("bEbar");

    Field<scalar>& epsilonPEq = state.scalarField("epsilonPEq");
    const Field<scalar>& epsilonPEq0 = state.getScalarField0("epsilonPEq");

    Field<scalar>& sigmaY = state.scalarField("sigmaY");

    // Optional scalar tangent
    const bool needScalarTan = response.wantsScalarTangent();
    UIndirectList<scalar>* scalarTanPtr = nullptr;
    if (needScalarTan)
    {
        scalarTanPtr = &response.scalarTangent();
    }

    const label nIP = sigma.size();

    // The Newton residual is normalised by the largest trial elastic strain,
    // so that the tolerance means the same at every integration point. It is
    // a tolerance scale only, so evaluating it over this law's own points
    // rather than the whole mesh is sufficient
    scalar maxMagBE = SMALL;
    if (nonLinearPlasticity_)
    {
        for (label i = 0; i < nIP; ++i)
        {
            const scalar relJi = J[i]/J0[i];
            const tensor relFbari
            (
                pow(relJi, -1.0/3.0)*(F[i] & Finv0[i])
            );

            maxMagBE =
                max
                (
                    maxMagBE,
                    mag(symm(relFbari & bEbar0[i] & relFbari.T()))
                );
        }

        reduce(maxMagBE, maxOp<scalar>());
        maxMagBE = max(maxMagBE, SMALL);
    }

    for (label i = 0; i < nIP; ++i)
    {
        const scalar Ji = J[i];
        const scalar relJi = Ji/J0[i];

        // Relative deformation gradient with the volumetric part removed
        const tensor relFbari(pow(relJi, -1.0/3.0)*(F[i] & Finv0[i]));

        // Trial isochoric elastic left Cauchy-Green tensor
        const symmTensor bEbarTrial
        (
            symm(relFbari & bEbar0[i] & relFbari.T())
        );

        // Trial deviatoric Kirchhoff stress
        const symmTensor sTrial(mu*dev(bEbarTrial));

        const scalar Ibar = tr(bEbarTrial)/3.0;
        const scalar muBar = Ibar*mu;

        const scalar magSTrial = mag(sTrial);

        // Yield stress at the start of the step, as a Cauchy stress
        const scalar sigmaYOld =
            stressPlasticStrainSeries_(max(epsilonPEq0[i], SMALL));

        // Trial yield function, in Kirchhoff stress
        const scalar fTrial =
            magSTrial - sqrt(2.0/3.0)*Ji*sigmaYOld;

        // Return direction
        symmTensor plasticN(symmTensor::zero);
        if (magSTrial > SMALL)
        {
            plasticN = sTrial/magSTrial;
        }

        scalar dLambda = 0.0;
        scalar curSigmaY = sigmaYOld;

        if (fTrial >= SMALL)
        {
            if (nonLinearPlasticity_)
            {
                newtonLoop
                (
                    dLambda,
                    curSigmaY,
                    epsilonPEq0[i],
                    magSTrial,
                    muBar,
                    Ji,
                    maxMagBE
                );
            }
            else
            {
                // Linear hardening has a closed-form return
                dLambda = fTrial/(2.0*muBar);

                if (mag(Hp_) > SMALL)
                {
                    dLambda /= 1.0 + Hp_/(3.0*muBar);
                    curSigmaY = sigmaYOld + sqrt(2.0/3.0)*dLambda*Hp_;
                }
            }
        }

        // Plastic strain increment
        const scalar dEpsilonPEq = sqrt(2.0/3.0)*dLambda;
        const symmTensor dEpsilonP(Ibar*dLambda*plasticN);

        // Deviatoric Kirchhoff stress
        const symmTensor sDev(sTrial - 2.0*mu*dEpsilonP);

        // Update the isochoric elastic left Cauchy-Green tensor
        const symmTensor devBEbar(sDev/mu);

        if (updateBEbarConsistent_)
        {
            bEbar[i] = devBEbar + IbarConsistent(devBEbar)*symmTensor(I);
        }
        else
        {
            bEbar[i] = devBEbar + Ibar*symmTensor(I);
        }

        // Hydrostatic stress, unsmoothed: see the note in the header
        const scalar p = 0.5*kappa*(sqr(Ji) - 1.0);

        sigma[i] = (1.0/Ji)*(p*symmTensor(I) + sDev);

        // Commit the history
        epsilonPEq[i] = epsilonPEq0[i] + dEpsilonPEq;
        sigmaY[i] = curSigmaY;

        if (needScalarTan)
        {
            // The legacy law scales the elastic stiffness down by how far the
            // return mapping moved the stress, which keeps the Laplacian
            // coefficient representative once a point is yielding
            const scalar scaleFactor =
                1.0 - 2.0*muBar*dLambda/max(magSTrial, SMALL);

            (*scalarTanPtr)[i] = scaleFactor*(4.0/3.0)*mu + kappa;
        }
    }

    // Fourth-order tangent by finite differences, as for the other
    // finite-strain laws. There is no analytical spatial tangent here
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


void Foam::neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw::endTimeStep
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
