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

#include "viscousHookeanElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "mat66.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscousHookeanElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        viscousHookeanElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * * //

Foam::word Foam::viscousHookeanElasticMechanicalConstitutiveLaw::arm
(
    const label i
)
{
    return "h" + Foam::name(i);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscousHookeanElasticMechanicalConstitutiveLaw::
viscousHookeanElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    EInf_(dict.lookup("EInfinity")),
    E_(dict.lookup("E")),
    tau_(dict.lookup("relaxationTimes")),
    nu_(dict.lookup("nu")),
    gammaInf_("gammaInf", dimless, 0.0),
    gamma_(E_.size(), 0.0),
    williamsLandelFerryShift_
    (
        dict.lookupOrDefault<Switch>("WilliamsLandelFerry", false)
    ),
    C1_
    (
        williamsLandelFerryShift_
      ? dimensionedScalar(dict.subDict("WilliamsLandelFerryCoeffs").lookup("C1"))
      : dimensionedScalar("C1", dimless, 0.0)
    ),
    C2_
    (
        williamsLandelFerryShift_
      ? dimensionedScalar(dict.subDict("WilliamsLandelFerryCoeffs").lookup("C2"))
      : dimensionedScalar("C2", dimTemperature, 0.0)
    ),
    Tref_
    (
        williamsLandelFerryShift_
      ? dimensionedScalar
        (
            dict.subDict("WilliamsLandelFerryCoeffs").lookup("Tref")
        )
      : dimensionedScalar("Tref", dimTemperature, 0.0)
    ),
    lambda_("lambda", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    kappa_("kappa", dimPressure, 0.0)
{
    if (EInf_.dimensions() != dimPressure)
    {
        FatalIOErrorInFunction(dict)
            << "The long-term modulus EInfinity has incorrect dimensions. "
            << "Expected " << dimPressure << " but got " << EInf_.dimensions()
            << exit(FatalIOError);
    }

    if (E_.size() != tau_.size())
    {
        FatalIOErrorInFunction(dict)
            << "The 'E' and 'relaxationTimes' lists must have the same "
            << "length: got " << E_.size() << " and " << tau_.size()
            << exit(FatalIOError);
    }

    if (E_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "At least one Maxwell arm must be given in 'E'."
            << exit(FatalIOError);
    }

    forAll(tau_, i)
    {
        if (tau_[i] < SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Relaxation time " << i << " is " << tau_[i]
                << ": all relaxation times must be positive."
                << exit(FatalIOError);
        }

        if (E_[i] < SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Maxwell arm modulus " << i << " is " << E_[i]
                << ": all moduli must be positive."
                << exit(FatalIOError);
        }
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
            << " is too close to 0.5: the bulk modulus would be "
            << "ill-conditioned."
            << exit(FatalIOError);
    }

    // Relative moduli, from the instantaneous modulus E0
    const dimensionedScalar E0
    (
        EInf_ + dimensionedScalar("sumE", dimPressure, sum(E_))
    );

    gammaInf_ = EInf_/E0;

    forAll(gamma_, i)
    {
        gamma_[i] = E_[i]/E0.value();
    }

    // The deviatoric response is instantaneous and relaxes; the volumetric
    // response is the long-term one. This follows the legacy law exactly,
    // including that lambda is built from EInfinity rather than E0
    const Switch planeStress(dict.lookupOrDefault<Switch>("planeStress", false));

    mu_ = E0/(2.0*(1.0 + nu_));

    if (planeStress)
    {
        lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - nu_));
    }
    else
    {
        lambda_ = nu_*EInf_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    }

    kappa_ = lambda_ + (2.0/3.0)*gammaInf_*mu_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::viscousHookeanElasticMechanicalConstitutiveLaw::initialiseState
(
    mechanicalConstitutiveLawState& state
) const
{
    // Relaxing deviatoric stress, and one internal stress per Maxwell arm.
    // All start from an unstressed state, at both times
    state.symmTensorField("s") = symmTensor::zero;
    state.symmTensorField0("s") = symmTensor::zero;

    forAll(gamma_, i)
    {
        state.symmTensorField(arm(i)) = symmTensor::zero;
        state.symmTensorField0(arm(i)) = symmTensor::zero;
    }
}


void Foam::viscousHookeanElasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const UIndirectList<tensor>& gradD = kin.gradD();

    const scalar mu = mu_.value();
    const scalar kappa = kappa_.value();
    const scalar gammaInf = gammaInf_.value();

    // The time increment comes from the inputs object: it is not a
    // deformation measure, so it does not live in the kinematics
    const scalar dt = inputs.dt();

    // Relaxing deviatoric stress, current and old
    Field<symmTensor>& sDev = state.symmTensorField("s");
    const Field<symmTensor>& sDev0 = state.getSymmTensorField0("s");

    const label nArms = gamma_.size();

    // Temperature, only where the relaxation times are shifted with it
    const scalarField* TPtr = nullptr;
    if (williamsLandelFerryShift_)
    {
        // Required, so ask in the way that fails rather than reading zero
        TPtr = &inputs.getScalar("T");
    }

    const label nIP = sigma.size();

    for (label i = 0; i < nIP; ++i)
    {
        const symmTensor e(symm(gradD[i]));

        // Instantaneous deviatoric stress
        const symmTensor sNew(2.0*mu*dev(e));

        // Volumetric response is elastic and long-term
        symmTensor sigmaI(kappa*tr(e)*symmTensor(I) + gammaInf*sNew);

        const symmTensor dS(sNew - sDev0[i]);

        for (label m = 0; m < nArms; ++m)
        {
            Field<symmTensor>& h = state.symmTensorField(arm(m));
            const Field<symmTensor>& h0 = state.getSymmTensorField0(arm(m));

            // Shift the relaxation time with temperature where asked
            scalar tau = tau_[m];

            if (williamsLandelFerryShift_)
            {
                const scalar dT = (*TPtr)[i] - Tref_.value();

                const scalar aT =
                    pow(10.0, -C1_.value()*dT/(C2_.value() + dT));

                tau *= aT;
            }

            // Simo & Hughes (1998), eq. 10.3.12
            h[i] = exp(-dt/tau)*h0[i] + exp(-dt/(2.0*tau))*dS;

            sigmaI += gamma_[m]*h[i];
        }

        sigma[i] = sigmaI;
        sDev[i] = sNew;
    }

    // Scalar tangent, if asked for. It depends on the time increment through
    // the relaxation factor, so this law's Jacobian is a function of dt
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        scalar scaleFactor = gammaInf;
        forAll(gamma_, m)
        {
            scaleFactor += gamma_[m]*exp(-dt/(2.0*tau_[m]));
        }

        scalar Keff = 0.0;

        switch (response.tangentReq())
        {
            case tangentRequest::scalar:
                Keff = scaleFactor*2.0*mu + lambda_.value();
                break;

            case tangentRequest::scalarDeviatoric:
                Keff = scaleFactor*(4.0/3.0)*mu;
                break;

            default:
                break;
        }

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }

    // Fourth-order tangent by finite differences
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


// ************************************************************************* //
