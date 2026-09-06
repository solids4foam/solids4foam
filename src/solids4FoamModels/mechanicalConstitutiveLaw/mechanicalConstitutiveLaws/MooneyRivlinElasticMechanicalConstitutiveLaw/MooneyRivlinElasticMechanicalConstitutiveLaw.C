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

#include "MooneyRivlinElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MooneyRivlinElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        MooneyRivlinElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MooneyRivlinElasticMechanicalConstitutiveLaw::
MooneyRivlinElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    c10_(dict.lookup("c10")),
    c01_(dict.lookup("c01")),
    c11_(dict.lookup("c11")),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    kappa_("kappa", dimPressure, 0.0)
{
    if
    (
        c10_.dimensions() != dimPressure
     || c01_.dimensions() != dimPressure
     || c11_.dimensions() != dimPressure
    )
    {
        FatalIOErrorInFunction(dict)
            << "The Mooney-Rivlin coefficients c10, c01 and c11 must all have "
            << "dimensions " << dimPressure
            << exit(FatalIOError);
    }

    // Shear modulus, from the pure shear stress state, as in the legacy law
    mu_ = 2.0*(c10_ + c01_);

    // The bulk modulus is given directly or derived from nu.
    //
    // The legacy law takes K when it is present and ignores nu, so an existing
    // case dictionary carrying both must keep working and keep the same
    // meaning. Rejecting that combination here would be tidier but would break
    // dictionaries the legacy law accepts, so it is warned about instead
    const bool haveK = dict.found("K");
    const bool haveNu = dict.found("nu");

    if (!haveK && !haveNu)
    {
        FatalIOErrorInFunction(dict)
            << "Specify either the bulk modulus 'K' or Poisson's ratio 'nu'."
            << exit(FatalIOError);
    }

    if (haveK && haveNu)
    {
        WarningInFunction
            << "Both 'K' and 'nu' are given: K is used and nu is ignored, "
            << "as in the legacy MooneyRivlinElastic law." << endl;
    }

    if (haveK)
    {
        kappa_ = dimensionedScalar(dict.lookup("K"));

        if (kappa_.dimensions() != dimPressure)
        {
            FatalIOErrorInFunction(dict)
                << "The bulk modulus K must have dimensions " << dimPressure
                << exit(FatalIOError);
        }
    }
    else
    {
        nu_ = dimensionedScalar(dict.lookup("nu"));

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
                << " is too close to 0.5: the bulk modulus derived from it "
                << "would be ill-conditioned. Give K directly instead."
                << exit(FatalIOError);
        }

        // As in the legacy law: E follows from the small-strain limit
        E_ = 6.0*(c10_ + c01_);
        kappa_ = E_/(3.0*(1.0 - 2.0*nu_));
    }

    // Retained for the scalar tangent, which is expressed through lambda
    lambda_ = kappa_ - (2.0/3.0)*mu_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::MooneyRivlinElasticMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const UIndirectList<tensor>& F = kin.F();
    const UIndirectList<scalar>& J = kin.J();

    // Whether the caller wants the isochoric stress and the volumetric
    // response separately, as a mixed displacement-pressure formulation does
    const bool wantsSplit = response.wantsVolumetricSplit();
    UIndirectList<scalar>* volumetricPtr =
        wantsSplit ? &response.volumetric() : nullptr;

    const scalar c10 = c10_.value();
    const scalar c01 = c01_.value();
    const scalar c11 = c11_.value();
    const scalar kappaVal = kappa_.value();

    // Not merely a positivity check: the stress divides by J and forms
    // inv(isoB), so a J this close to zero is already meaningless
    // numerically. The legacy law imposes no such floor
    const scalar Jmin = sqrt(SMALL);

    forAll(sigma, i)
    {
        const tensor& Fi = F[i];
        const scalar Ji = J[i];

        if (Ji <= Jmin)
        {
            FatalErrorInFunction
                << "Deformation gradient determinant J = " << Ji
                << " at index " << i << " is not usable: J must be positive "
                << "and greater than " << Jmin << ", since the stress divides "
                << "by J. A J this small means the deformation has already "
                << "degenerated."
                << exit(FatalError);
        }

        // Isochoric left Cauchy-Green tensor and its invariants
        const symmTensor isoB(pow(Ji, -2.0/3.0)*symm(Fi & Fi.T()));
        const symmTensor sqrB(symm(isoB & isoB));

        const scalar I1 = tr(isoB);
        const scalar I2 = 0.5*(sqr(I1) - tr(sqrB));

        // Isochoric stress
        const symmTensor s
        (
            2.0*(c10 + c11*(I2 - 3.0))*isoB
          - 2.0*(c01 + c11*(I1 - 3.0))*inv(isoB)
        );

        // Hydrostatic stress, unsmoothed: see the note in the header. The
        // legacy law passes this through updateSigmaHyd, which may solve a
        // Laplacian over the field; that belongs to the solid model, not here
        const scalar p = 0.5*kappaVal*(sqr(Ji) - 1.0);

        sigma[i] = (1.0/Ji)*(dev(s) + p*I);

        // A caller that asked for the two apart gets them apart. dev(s)/J is
        // trace free, so for this law the split and a deviatoric projection of
        // the total agree; that is not true of every law, which is why the
        // caller asks rather than projecting
        if (wantsSplit)
        {
            (*volumetricPtr)[i] = p/Ji;
            sigma[i] = (1.0/Ji)*dev(s);
        }
    }

    // Scalar tangent: only if explicitly requested
    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        scalar Keff = 0.0;

        switch (response.tangentReq())
        {
            case tangentRequest::scalar:
                // Matches the implicit stiffness the legacy law passes to
                // updateSigmaHyd, which is (4/3)*mu + K
                Keff = (4.0/3.0)*mu_.value() + kappa_.value();
                break;

            case tangentRequest::scalarDeviatoric:
                // Scalar Laplacian surrogate for div(dev(sigma)), which is
                // mu*lap(D) + (1/3)*mu*grad(div(D))
                Keff = (4.0/3.0)*mu_.value();
                break;

            default:
                break;
        }

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }

    // Fourth-order tangent.
    // There is no analytical spatial tangent for this law yet, but the
    // finite-difference tangent of the base class is well defined for any
    // finite-strain law and is evaluated against a shadow state, so it leaves
    // neither the stress just computed nor the history it started from
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
