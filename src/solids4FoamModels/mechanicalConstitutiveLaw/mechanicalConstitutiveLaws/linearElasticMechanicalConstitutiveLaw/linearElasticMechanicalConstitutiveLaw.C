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

#include "linearElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "mat66.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        linearElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticMechanicalConstitutiveLaw::
linearElasticMechanicalConstitutiveLaw
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

void Foam::linearElasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    mechanicalConstitutiveLawState& /*state*/,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const UIndirectList<tensor>& gradD = kin.gradD();

    // Field approach: this creates intermiedate temporary fields
    // sigma = mu_*twoSymm(gradD) + lambda_*tr(gradD)*I;

    // Element-by-element approach: faster as it avoid intermediate fields
    const scalar muVal = mu_.value();
    const scalar lambdaVal = lambda_.value();
    forAll(sigma, i)
    {
        sigma[i] = muVal*twoSymm(gradD[i]) + lambdaVal*tr(gradD[i])*I;
    }

    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        scalar Keff = 0.0;

        switch (response.tangentReq())
        {
            case tangentRequest::scalar:
                Keff = 2.0*mu_.value() + lambda_.value();
                break;

            case tangentRequest::scalarDeviatoric:
                Keff = 2.0*mu_.value();
                break;

            default:
                break;
        }

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }
    else if (response.hasFourthOrderTangent())
    {
        UIndirectList<mat66>& C = response.fourthOrderTangent();

        forAll(C, i)
        {
            // Take a reference to C for the current integration point
            mat66& curC = C[i];

            // Important: zero everything first as mat66 is not initialised by
            // default
            curC.clear();

            // Normal-normal blocks
            curC(0,0) = lambdaVal + 2.0*muVal;
            curC(1,1) = lambdaVal + 2.0*muVal;
            curC(2,2) = lambdaVal + 2.0*muVal;

            curC(0,1) = lambdaVal;
            curC(0,2) = lambdaVal;
            curC(1,0) = lambdaVal;
            curC(1,2) = lambdaVal;
            curC(2,0) = lambdaVal;
            curC(2,1) = lambdaVal;

            // Shear terms (engineering shear)
            curC(3,3) = muVal;
            curC(4,4) = muVal;
            curC(5,5) = muVal;
        }
    }
}


// ************************************************************************* //
