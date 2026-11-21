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

#include "linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw.H"
#include "misesPlasticDiagnostics.H"
#include "addToRunTimeSelectionTable.H"
#include "SubList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw, 1
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


const Foam::word
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
phaseFieldName = "phaseField";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
yieldStress
(
    const scalar epsilonPEq
) const
{
    return stressPlasticStrainSeries_(max(epsilonPEq, SMALL));
}


Foam::scalar
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
plasticMultiplier
(
    const scalar qTrial,
    const scalar epsilonPEq0,
    const scalar maxMagDEpsilon
) const
{
    const scalar mu = mu_.value();
    const scalar fTrial = qTrial - yieldStress(epsilonPEq0);

    if (fTrial <= SMALL)
    {
        return 0;
    }

    if (!nonLinearPlasticity_)
    {
        return fTrial/(3.0*mu + Hp_);
    }

    //     qTrial - 3 mu DLambda - sigmaY(r0 + DLambda) = 0
    scalar DLambda = 0;
    scalar f = fTrial;
    label iter = 0;

    do
    {
        const scalar dSigmaY =
            (
                yieldStress(epsilonPEq0 + DLambda + finiteDiff_)
              - yieldStress(epsilonPEq0 + DLambda)
            )/finiteDiff_;

        const scalar dLambda = f/(3.0*mu + dSigmaY);
        DLambda += dLambda;

        f = qTrial - 3.0*mu*DLambda - yieldStress(epsilonPEq0 + DLambda);

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

    return DLambda;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw
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
    Gc_(0),
    l_(0),
    w0_(0),
    k_(dict.lookupOrDefault<scalar>("residualStiffness", 0.0)),
    loopTol_(dict.lookupOrDefault<scalar>("NewtonLoopTol", 1e-8)),
    maxNewtonIter_(dict.lookupOrDefault<label>("NewtonMaxIter", 200)),
    finiteDiff_(dict.lookupOrDefault<scalar>("NewtonFiniteDiffEps", 0.25e-6))
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

    const dimensionedScalar Gc(dict.lookup("Gc"));
    const dimensionedScalar l(dict.lookup("characteristicLength"));

    if
    (
        Gc.dimensions() != dimPressure*dimLength
     || Gc.value() <= 0
     || l.dimensions() != dimLength
     || l.value() <= 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "Gc and characteristicLength must be positive, with dimensions "
            << dimPressure*dimLength << " and " << dimLength
            << exit(FatalIOError);
    }

    Gc_ = Gc.value();
    l_ = l.value();

    if (dict.found("w0"))
    {
        const dimensionedScalar w0(dict.lookup("w0"));

        if (w0.dimensions() != dimPressure || w0.value() < 0)
        {
            FatalIOErrorInFunction(dict)
                << "w0 must not be negative, with dimensions " << dimPressure
                << exit(FatalIOError);
        }

        w0_ = w0.value();
    }

    if (k_ < 0 || k_ >= 1)
    {
        FatalIOErrorInFunction(dict)
            << "The residual stiffness " << k_ << " must lie in [0, 1)"
            << exit(FatalIOError);
    }

    Info<< "    Phase field fracture: Gc = " << Gc_
        << ", characteristicLength = " << l_ << ", w0 = " << w0_
        << ", residualStiffness = " << k_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
degradation
(
    const scalar d
) const
{
    return (1.0 - k_)*sqr(1.0 - d) + k_;
}


void
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    typedef mechanicalConstitutiveLawStateSpec::stateRole role;

    spec.addSymmTensor("epsilonP", role::persistent, symmTensor::zero);

    spec.addScalar("epsilonPEq", role::persistent, 0.0);

    spec.addScalar("sigmaY", role::persistent, stressPlasticStrainSeries_(0.0));

    // The history of the tensile elastic energy, and the plastic work
    spec.addScalar("elasticEnergyHistory", role::persistent, 0.0);

    spec.addScalar("plasticWork", role::persistent, 0.0);

    // The crack driving force A of the phase field equation, and its
    // coefficient 1 + A
    spec.addScalar("phaseFieldSource", role::persistent, 0.0);

    spec.addScalar("phaseFieldCoefficient", role::persistent, 1.0);
}


Foam::scalar
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
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
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
requiredScalarInputs() const
{
    return wordList(1, phaseFieldName);
}


Foam::List<Foam::mechanicalConstitutiveLaw::implicitGradientField>
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
implicitGradientFields() const
{
    // (1 + A) d - div((2 l)^2 grad(d)) = A
    List<implicitGradientField> fields(1);
    fields[0].name = phaseFieldName;
    fields[0].source = "phaseFieldSource";
    fields[0].length = 2.0*l_;
    fields[0].coefficient = "phaseFieldCoefficient";

    return fields;
}


void
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
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


void
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
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

    const Field<scalar>& H0 = state.getScalarField0("elasticEnergyHistory");
    Field<scalar>& H = state.scalarField("elasticEnergyHistory");

    const Field<scalar>& Wp0 = state.getScalarField0("plasticWork");
    Field<scalar>& Wp = state.scalarField("plasticWork");

    Field<scalar>& source = state.scalarField("phaseFieldSource");
    Field<scalar>& coefficient = state.scalarField("phaseFieldCoefficient");

    const scalarField& d = inputs.getScalar(phaseFieldName);

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

    forAll(sigma, i)
    {
        // The return of the undegraded material: the yield stress is
        // degraded as the stress is, so the degradation cancels
        const symmTensor sTrial(2.0*mu*dev(epsilonETrial[i]));
        const scalar magSTrial = mag(sTrial);
        const scalar qTrial = sqrt(3.0/2.0)*magSTrial;

        const scalar DLambda =
            plasticMultiplier(qTrial, epsilonPEq0[i], maxMagDEpsilon);

        symmTensor n(symmTensor::zero);
        if (magSTrial > SMALL)
        {
            n = sTrial/magSTrial;
        }

        epsilonE[i] = epsilonETrial[i] - (sqrt(3.0/2.0)*DLambda)*n;
        epsilonPEq[i] = epsilonPEq0[i] + DLambda;
        sigmaY[i] = yieldStress(epsilonPEq[i]);

        // The degraded stress, whose volumetric part is degraded only in
        // tension
        const scalar g = degradation(min(max(d[i], 0.0), 1.0));
        const scalar trEpsilonE = tr(epsilonE[i]);
        const scalar gV = trEpsilonE >= 0 ? g : 1.0;

        sigma[i] =
            g*2.0*mu*dev(epsilonE[i])
          + gV*kappa*trEpsilonE*symmTensor(I);

        // The crack driving force, of the undegraded energies
        const scalar psiPos =
            mu*(dev(epsilonE[i]) && dev(epsilonE[i]))
          + 0.5*kappa*sqr(max(trEpsilonE, 0.0));

        H[i] = max(H0[i], psiPos);
        Wp[i] = Wp0[i] + sigmaY[i]*DLambda;

        // The driving force is the energy release of g(d), whose slope
        // carries the factor 1 - k
        source[i] = 4.0*l_*(1.0 - k_)*(H[i] + max(Wp[i] - w0_, 0.0))/Gc_;
        coefficient[i] = 1.0 + source[i];

        if (needScalarTan)
        {
            scalar theta = 1.0;
            if (magSTrial > SMALL)
            {
                theta = max(1.0 - 3.0*mu*DLambda/qTrial, 0.0);
            }

            if (deviatoricTan)
            {
                response.scalarTangent()[i] = g*theta*(4.0/3.0)*mu;
            }
            else
            {
                response.scalarTangent()[i] =
                    g*theta*(4.0/3.0)*mu + gV*kappa;
            }
        }
    }
}


void
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
endTimeStep
(
    const mechanicalConstitutiveLawState& state,
    const scalar time,
    const label timeIndex,
    DynamicList<mechanicalConstitutiveLawDiagnostic>& diagnostics
) const
{
    appendMisesPlasticDiagnostics(state, diagnostics);

    const Field<scalar>& source = state.scalarField("phaseFieldSource");

    scalar maxSource = 0;

    forAll(source, i)
    {
        maxSource = max(maxSource, source[i]);
    }

    typedef mechanicalConstitutiveLawDiagnostic diagnostic;

    // The local phase field of a uniform field, A/(1 + A), which the phase
    // field approaches where the crack driving force is large
    diagnostics.append
    (
        diagnostic
        (
            "max local phase field",
            maxSource/(1.0 + maxSource),
            diagnostic::combineOperation::maximum
        )
    );
}


void
Foam::linearElasticMisesPlasticPhaseFieldFractureMechanicalConstitutiveLaw::
reportDiagnostics
(
    const UList<mechanicalConstitutiveLawDiagnostic>& diagnostics
) const
{
    if (!debug || diagnostics.size() != 4)
    {
        return;
    }

    reportMisesPlasticDiagnostics
    (
        debug, SubList<mechanicalConstitutiveLawDiagnostic>(diagnostics, 3)
    );

    Info<< "Max local phase field A/(1 + A) = " << diagnostics[3].value()
        << endl;
}


// ************************************************************************* //
