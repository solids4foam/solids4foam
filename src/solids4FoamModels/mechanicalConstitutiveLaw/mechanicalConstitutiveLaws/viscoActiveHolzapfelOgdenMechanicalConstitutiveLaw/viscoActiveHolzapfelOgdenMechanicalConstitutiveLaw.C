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

#include "viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar
Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::boundedExp
(
    const scalar x
)
{
    // exp() with its argument clamped, so that an extreme trial state in a
    // Newton line search produces a large, finite stress rather than inf
    return exp(min(x, scalar(500)));
}


inline void
Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::fibreSwitch
(
    const scalar I4,
    scalar& H,
    scalar& dH
) const
{
    switch (fibreSwitch_)
    {
        case fibreSwitchType::smooth:
        {
            // 1/(1 + exp(-k(x - 1))), written so that neither branch
            // overflows
            const scalar z = k_*(I4 - 1.0);
            if (z >= 0)
            {
                H = 1.0/(1.0 + exp(-z));
            }
            else
            {
                const scalar ez = exp(z);
                H = ez/(1.0 + ez);
            }
            dH = k_*H*(1.0 - H);
            break;
        }
        case fibreSwitchType::sharp:
        {
            H = (I4 > 1.0) ? 1.0 : 0.0;
            dH = 0.0;
            break;
        }
        default:
        {
            H = 1.0;
            dH = 0.0;
        }
    }
}


inline Foam::scalar
Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::dPsidI4
(
    const scalar I4,
    const scalar a,
    const scalar b
) const
{
    // Psi_4 = a/(2b) H(I4) (exp(b (I4 - 1)^2) - 1)
    scalar H = 0;
    scalar dH = 0;
    fibreSwitch(I4, H, dH);

    const scalar E4 = boundedExp(b*sqr(I4 - 1.0));

    return 0.5*a/b*dH*(E4 - 1.0) + a*H*(I4 - 1.0)*E4;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::
viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    a_(dict.lookup("a")),
    b_(readScalar(dict.lookup("b"))),
    af_(dict.lookup("af")),
    bf_(readScalar(dict.lookup("bf"))),
    as_(dict.lookup("as")),
    bs_(readScalar(dict.lookup("bs"))),
    afs_(dict.lookup("afs")),
    bfs_(readScalar(dict.lookup("bfs"))),
    bulkModulus_(dict.lookup("bulkModulus")),
    eta_
    (
        dict.found("eta")
      ? dimensionedScalar(dict.lookup("eta"))
      : dimensionedScalar("eta", dimPressure*dimTime, 0.0)
    ),
    fibreSwitch_(fibreSwitchType::smooth),
    k_(dict.lookupOrDefault<scalar>("k", 100.0)),
    tangentShearModulus_
    (
        dict.lookupOrDefault<scalar>
        (
            "tangentShearModulus",
            a_.value() + (af_.value() + as_.value() + afs_.value())/6.0
        )
    ),
    activeStressSeries_(),
    f0Default_
    (
        dict.lookupOrDefault<Switch>("uniformFibreField", false)
      ? vector(dict.lookup("f0"))
      : vector::zero
    ),
    s0Default_
    (
        dict.lookupOrDefault<Switch>("uniformFibreField", false)
      ? vector(dict.lookup("s0"))
      : vector::zero
    )
{
    const word switchName
    (
        dict.lookupOrDefault<word>("fibreSwitch", "smooth")
    );

    if (switchName == "smooth")
    {
        fibreSwitch_ = fibreSwitchType::smooth;
    }
    else if (switchName == "sharp")
    {
        fibreSwitch_ = fibreSwitchType::sharp;
    }
    else if (switchName == "none")
    {
        fibreSwitch_ = fibreSwitchType::none;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown fibreSwitch '" << switchName << "'. Valid options are "
            << "smooth, sharp and none." << exit(FatalIOError);
    }

    if (bulkModulus_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "bulkModulus must be positive: it is the penalty that keeps "
            << "this material near incompressible." << exit(FatalIOError);
    }

    if (b_ <= 0 || bf_ <= 0 || bs_ <= 0 || bfs_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "The exponents b, bf, bs and bfs must be positive."
            << exit(FatalIOError);
    }

    if (tangentShearModulus_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "tangentShearModulus must be positive." << exit(FatalIOError);
    }

    if (eta_.value() < 0)
    {
        FatalIOErrorInFunction(dict)
            << "The viscosity eta must not be negative." << exit(FatalIOError);
    }

    if (dict.found("activeStressSeries"))
    {
        activeStressSeries_ =
            interpolationTable<scalar>(dict.subDict("activeStressSeries"));
    }

    Info<< "    Visco-active Holzapfel-Ogden law: a = " << a_.value()
        << ", af = " << af_.value() << ", as = " << as_.value()
        << ", afs = " << afs_.value() << ", kappa = " << bulkModulus_.value()
        << ", eta = " << eta_.value() << ", fibreSwitch = " << switchName
        << ", tangentShearModulus = " << tangentShearModulus_;

    if (activeStressSeries_.size())
    {
        Info<< ", active stress from a table of "
            << activeStressSeries_.size() << " points";
    }
    else
    {
        Info<< ", passive";
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::activeStress
(
    const scalar t
) const
{
    if (activeStressSeries_.size())
    {
        return activeStressSeries_(t);
    }

    return 0.0;
}


void Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // Fibre and sheet directions in the reference configuration. Prescribed:
    // the case supplies them and this law only reads them
    spec.addVector
    (
        "f0",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        f0Default_
    );

    spec.addVector
    (
        "s0",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        s0Default_
    );
}


void Foam::viscoActiveHolzapfelOgdenMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();
    const UIndirectList<tensor>& F = kin.F();
    const UIndirectList<tensor>& F0 = kin.F0();
    const UIndirectList<scalar>& J = kin.J();

    // Prescribed fields are read at old time: see
    // mechanicalConstitutiveLawStateSpec
    const Field<vector>& f0 = state.getVectorField0("f0");
    const Field<vector>& s0 = state.getVectorField0("s0");

    const scalar a = a_.value();
    const scalar af = af_.value();
    const scalar as = as_.value();
    const scalar afs = afs_.value();
    const scalar kappa = bulkModulus_.value();
    const scalar eta = eta_.value();

    // The active stress is uniform in space: evaluate it once
    const scalar tau = activeStress(inputs.time());

    // The viscous term needs a time step; it is skipped before the first one
    const scalar dt = inputs.dt();
    const bool viscous = (eta > 0 && dt > VSMALL);

    // Whether the caller wants the volumetric response apart, as a mixed
    // displacement-pressure formulation does
    const bool wantsSplit = response.wantsVolumetricSplit();
    UIndirectList<scalar>* volumetricPtr =
        wantsSplit ? &response.volumetric() : nullptr;

    forAll(sigma, i)
    {
        const scalar Ji = J[i];

        if (Ji <= sqrt(SMALL))
        {
            FatalErrorInFunction
                << "Invalid deformation gradient determinant J = " << Ji
                << " at index " << i << '.' << exit(FatalError);
        }

        const scalar magF0 = mag(f0[i]);
        const scalar magS0 = mag(s0[i]);

        if (magF0 < SMALL || magS0 < SMALL)
        {
            FatalErrorInFunction
                << "The fibre or sheet direction has zero length at index "
                << i << '.' << nl
                << "This law needs both. Either supply f0 and s0 fields - "
                << "setEllipsoidFibreField writes them - or set "
                << "'uniformFibreField yes' and give 'f0' and 's0' in this "
                << "material's dictionary."
                << exit(FatalError);
        }

        const tensor& Fi = F[i];
        const scalar rJ = 1.0/Ji;

        // Isotropic part, on the isochoric invariant:
        //     sigma = (2/J) dPsi/dI1bar dev(bbar),  bbar = J^(-2/3) F.F^T
        const symmTensor bBar(pow(Ji, -2.0/3.0)*symm(Fi & Fi.T()));
        const scalar I1bar = tr(bBar);
        const scalar psi1 = 0.5*a*boundedExp(b_*(I1bar - 3.0));

        symmTensor sigmaI(2.0*rJ*psi1*dev(bBar));

        // Pushed-forward fibre and sheet directions
        const vector ff(Fi & (f0[i]/magF0));
        const vector ss(Fi & (s0[i]/magS0));

        // Fibre and sheet families, on the full invariants:
        //     sigma = (2/J) dPsi/dI4 (F f0) (x) (F f0)
        const scalar I4f = ff & ff;
        const scalar I4s = ss & ss;

        sigmaI += 2.0*rJ*dPsidI4(I4f, af, bf_)*sqr(ff);
        sigmaI += 2.0*rJ*dPsidI4(I4s, as, bs_)*sqr(ss);

        // Fibre-sheet coupling:
        //     sigma = (1/J) dPsi/dI8 ((F f0) (x) (F s0) + (F s0) (x) (F f0))
        const scalar I8fs = ff & ss;
        const scalar dPsidI8 = afs*I8fs*boundedExp(bfs_*sqr(I8fs));

        sigmaI += 2.0*rJ*dPsidI8*symm(ff*ss);

        // Volumetric penalty, dU/dJ of U = kappa/4 (J^2 - 1 - 2 ln J). Kept
        // apart until the end, so that a mixed formulation can be handed it
        // separately
        // Written exactly as the mixed solver's pressure equation writes it,
        // 0.5*K*(J^2 - 1)/J, rather than the equal 0.5*K*(J - 1/J): near
        // J = 1 the two differ by round-off relative to their tiny value,
        // and the solver's volumetric closure check compares them
        const scalar dUdJ = 0.5*kappa*(sqr(Ji) - 1.0)/Ji;

        // Active stress: S = tau f0 (x) f0, pushed forward
        if (mag(tau) > 0)
        {
            sigmaI += rJ*tau*sqr(ff);
        }

        // Viscous stress: S = eta Edot, with Edot over the time step
        if (viscous)
        {
            const tensor& F0i = F0[i];
            const symmTensor Edot
            (
                0.5*(symm(Fi.T() & Fi) - symm(F0i.T() & F0i))/dt
            );

            sigmaI += rJ*eta*symm(Fi & Edot & Fi.T());
        }

        if (wantsSplit)
        {
            // The mixed displacement-pressure formulation replaces dU/dJ by
            // its solved pressure, which satisfies the same penalty relation,
            // and adds it to the rest of the stress: the material is the same
            (*volumetricPtr)[i] = dUdJ;
            sigma[i] = sigmaI;
        }
        else
        {
            sigma[i] = sigmaI + dUdJ*symmTensor(I);
        }
    }

    // Scalar tangent, if asked for: a small-strain shear-modulus estimate.
    // It is not a tangent of this energy, but it is more than a
    // preconditioner: the solid model also scales its FV stabilisation by it,
    // and that stabilisation adds a stiffness of its own, of order
    // (scaleFactor * this) * h^2 * curvature. The relaxed matrix is very soft
    // (a = 59 Pa) and the fibres are slack in compression, so an estimate
    // far above the physical stiffness makes the wall artificially stiff.
    // The default is the linearised passive shear modulus,
    // a + (af + as + afs)/6, without the viscous eta/dt; it may be set
    // directly with tangentShearModulus
    const scalar mu = tangentShearModulus_;

    fillScalarTangent
    (
        response,
        (4.0/3.0)*mu + kappa,
        (4.0/3.0)*mu
    );

    // No analytical fourth-order tangent has been derived for this law
    fourthOrderByFiniteDifferenceOnly(kin, inputs, state, response);
}


// ************************************************************************* //
