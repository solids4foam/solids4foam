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
#include "Switch.H"

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
    rho_(dict.lookup("rho")),
    E_("E", dimPressure, 0.0),
    nu_("nu", dimless, 0.0),
    lambda_("lambda", E_.dimensions(), 0.0),
    mu_("mu", E_.dimensions(), 0.0),
    kappa_("kappa", E_.dimensions(), 0.0)
{
    // The material may be given either as E and nu or as mu and K, matching
    // the legacy linearElastic law, so an existing case dictionary needs no
    // change. The legacy law tries E and nu first, so this does too
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

    // Set mu, lambda and kappa
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

void Foam::linearElasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
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
    else if
    (
        response.tangentReq()
     == tangentRequest::fourthOrderFiniteDifference
    )
    {
        if (!response.hasFourthOrderTangent())
        {
            FatalErrorInFunction
                << "Finite difference fourth order tangent requested but the "
                << "response was not provided the with field"
                << exit(FatalError);
        }

        // Base-class finite difference implementation
        finiteDifferenceFourthOrder(kin, inputs, state, response);
    }
    else if (response.tangentReq() == tangentRequest::fourthOrder)
    {
        if (!response.hasFourthOrderTangent())
        {
            FatalErrorInFunction
                << "Fourth order tangent requested but the "
                << "response was not provided with the field"
                << exit(FatalError);
        }

        UIndirectList<mat66>& Cfield = response.fourthOrderTangent();

        // Define matrix indices for readability
        const label XX = symmTensor::XX;
        const label YY = symmTensor::YY;
        const label ZZ = symmTensor::ZZ;
        const label XY = symmTensor::XY;
        const label YZ = symmTensor::YZ;
        const label XZ = symmTensor::XZ;

        // Analytical calculate the tangent
        forAll(Cfield, i)
        {
            // Take a reference to C for the current integration point
            mat66& C = Cfield[i];

            // Important: zero everything first as mat66 is not initialised by
            // default
            C.clear();

            // Normal-normal blocks
            C(XX,XX) = lambdaVal + 2.0*muVal;
            C(YY,YY) = lambdaVal + 2.0*muVal;
            C(ZZ,ZZ) = lambdaVal + 2.0*muVal;

            C(XX,YY) = lambdaVal;
            C(XX,ZZ) = lambdaVal;
            C(YY,XX) = lambdaVal;
            C(YY,ZZ) = lambdaVal;
            C(ZZ,XX) = lambdaVal;
            C(ZZ,YY) = lambdaVal;

            // Shear terms (engineering shear)
            C(XY,XY) = muVal;
            C(YZ,YZ) = muVal;
            C(XZ,XZ) = muVal;
        }
    }
}


// ************************************************************************* //
