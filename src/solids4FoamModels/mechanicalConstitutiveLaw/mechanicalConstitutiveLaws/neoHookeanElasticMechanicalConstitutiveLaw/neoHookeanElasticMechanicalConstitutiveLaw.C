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

#include "neoHookeanElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        neoHookeanElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neoHookeanElasticMechanicalConstitutiveLaw::
neoHookeanElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", dimPressure, 0.0),
    mu_("mu", dimPressure, 0.0),
    kappa_("kappa", dimPressure, 0.0)
{
    // The material may be given either as E and nu or as mu and K, matching
    // the legacy neoHookeanElastic law, so that an existing case dictionary
    // needs no change. Exactly one of the two pairs must be present
    const bool haveENu = dict.found("E") && dict.found("nu");
    const bool haveMuK = dict.found("mu") && dict.found("K");

    if (haveENu == haveMuK)
    {
        FatalIOErrorInFunction(dict)
            << "Specify the elastic properties either as 'E' and 'nu' or as "
            << "'mu' and 'K', and not as both."
            << exit(FatalIOError);
    }

    if (haveMuK)
    {
        const dimensionedScalar mu(dict.lookup("mu"));
        const dimensionedScalar K(dict.lookup("K"));

        if (mu.dimensions() != dimPressure || K.dimensions() != dimPressure)
        {
            FatalIOErrorInFunction(dict)
                << "The shear modulus mu and bulk modulus K must both have "
                << "dimensions " << dimPressure
                << exit(FatalIOError);
        }

        // Invert mu = E/(2*(1 + nu)) and K = E/(3*(1 - 2*nu))
        E_ = 9.0*K*mu/(3.0*K + mu);
        nu_ = (3.0*K - 2.0*mu)/(2.0*(3.0*K + mu));
    }
    else
    {
        E_ = dimensionedScalar(dict.lookup("E"));
        nu_ = dimensionedScalar(dict.lookup("nu"));
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

    // Set lambda, mu and kappa
    // Note: the planeStress entry is injected into this dictionary by the
    // mechanicalConstitutiveLawManager from the top-level entry in
    // mechanicalProperties; it is not given by the user in this sub-dictionary
    const Switch planeStress(dict.lookupOrDefault<Switch>("planeStress", false));

    mu_ = E_/(2.0*(1.0 + nu_));

    if (planeStress)
    {
        lambda_ = E_*nu_/((1.0 + nu_)*(1.0 - nu_));
    }
    else
    {
        lambda_ = E_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    }

    // Note: for the three-dimensional case, this is equivalent to
    // E/(3*(1 - 2*nu))
    kappa_ = lambda_ + (2.0/3.0)*mu_;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::neoHookeanElasticMechanicalConstitutiveLaw::evaluate
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

    const scalar muVal = mu_.value();
    const scalar kappaVal = kappa_.value();

    const scalar Jmin = sqrt(SMALL);

    // Fast element-by-element approach
    forAll(sigma, i)
    {
        const scalar Ji = J[i];

        if (Ji <= Jmin)
        {
            FatalErrorInFunction
                << "Invalid deformation gradient determinant J = " << Ji
                << " at index " << i
                << ". J must be positive for log(J)."
                << exit(FatalError);
        }

        const symmTensor bEbar = pow(Ji, -2.0/3.0)*symm(F[i] & F[i].T());

        // Hydrostatic stress. The legacy neoHookeanElastic uses
        // 0.5*K*(J^2 - 1), and this must match it: the two volumetric
        // energies, 0.5*K*(J^2 - 1) and K*log(J), agree only to first order
        // in (J - 1), so they are indistinguishable near the identity and
        // diverge at finite strain. Using log(J) here made this law disagree
        // with the legacy one by about 1e-3 on neckingBar
        const scalar sigmaHyd = 0.5*kappaVal*(sqr(Ji) - 1.0);

        sigma[i] = (muVal/Ji)*dev(bEbar) + (sigmaHyd/Ji)*I;
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
