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
    rho_("rho", dict),
    E_("E", dict),
    nu_("nu", dict),
    lambda_("lambda", E_.dimensions(), 0.0),
    mu_("mu", E_.dimensions(), 0.0),
    kappa_("kappa", E_.dimensions(), 0.0)
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

    // Set lambda, mu and kappa
    lambda_ = E_*nu_/((1.0 + nu_)*(1.0 - 2.0*nu_));
    mu_ = E_/(2.0*(1.0 + nu_));
    kappa_ = E_/(3.0*(1 - 2*nu_));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::neoHookeanElasticMechanicalConstitutiveLaw::evaluate
(
    const finiteStrainMechanicalConstitutiveLawKinematics& kin,
    mechanicalConstitutiveLawState& /*state*/,
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

        sigma[i] = (muVal/Ji)*dev(bEbar) + (kappaVal/Ji)*log(Ji)*I;
    }

    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& impK = response.scalarTangent();
        const scalar Keff = (4.0/3.0)*muVal + kappaVal;

        forAll(impK, i)
        {
            impK[i] = Keff;
        }
    }
}


// ************************************************************************* //
