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

#include "isotropicFungElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isotropicFungElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        isotropicFungElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isotropicFungElasticMechanicalConstitutiveLaw::
isotropicFungElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    c_(dict.lookup("c")),
    a_(readScalar(dict.lookup("a"))),
    bulkModulus_(dict.lookup("bulkModulus")),
    mu_(0.5*c_*a_)
{
    if (bulkModulus_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "bulkModulus must be positive: it is the penalty that keeps "
            << "this material near incompressible." << exit(FatalIOError);
    }

    if (c_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "c must be positive: it is the stress scale "
            << "multiplying the exponential." << exit(FatalIOError);
    }

    if (a_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "a must be positive: it is the coefficient of "
            << "tr(Ebar & Ebar) in the exponent." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isotropicFungElasticMechanicalConstitutiveLaw::evaluate
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

    const scalar cVal = c_.value();
    const scalar bulkVal = bulkModulus_.value();

    forAll(sigma, i)
    {
        const scalar Ji = J[i];

        if (Ji <= sqrt(SMALL))
        {
            FatalErrorInFunction
                << "Invalid deformation gradient determinant J = " << Ji
                << " at index " << i << '.' << exit(FatalError);
        }

        // Green-Lagrange strain of the isochoric deformation, not of the
        // whole of it. Fbar = J^(-1/3)*F carries the shape change and none of
        // the volume change, so Q depends on shape alone and the volumetric
        // response below is the only place volume enters.
        //
        // Building Q from the full strain would make the energy coupled: the
        // deviatoric stress would then vary with J, and a mixed formulation
        // replacing the volumetric part would give a different material
        // rather than the same one solved differently. Written this way the
        // two formulations describe one material, and both reduce to the
        // published model in the incompressible limit it was defined for
        const tensor Fbar(pow(Ji, -1.0/3.0)*F[i]);
        const symmTensor E(0.5*(symm(Fbar.T() & Fbar) - I));
        const symmTensor sqrE(symm(E & E));

        const scalar Q = a_*tr(sqrE);

        const symmTensor dQdE(2*a_*E);

        // Second Piola-Kirchhoff stress, without the volumetric term
        const symmTensor S(dQdE*0.5*cVal*exp(Q));

        // Isochoric Cauchy stress. dev() applied after the push-forward is
        // DEV[] applied before it, since F & inv(C) & F.T() is I, and that
        // DEV[] is the chain rule through J^(-2/3), not a property S has
        const symmTensor s(dev(symm(Fbar & S & Fbar.T()))/Ji);

        // The volumetric response, dU/dJ, from the penalty that keeps this
        // material near incompressible: U(J) = 0.25*K*(J^2 - 1 - 2*log(J))
        const scalar dUdJ = 0.5*bulkVal*(sqr(Ji) - 1.0)/Ji;

        sigma[i] = s + dUdJ*I;

        // A caller that asked for the two apart gets them apart: the stress
        // above is already the isochoric part plus the volumetric one, so the
        // split costs a subtraction rather than a second evaluation
        if (wantsSplit)
        {
            (*volumetricPtr)[i] = dUdJ;
            sigma[i] = s;
        }
    }

    // Scalar tangent, if asked for: a small-strain estimate, which is a
    // preconditioner rather than a tangent of this energy.
    //
    // A mixed displacement-pressure formulation carries the volumetric
    // response in its own equation, so the deviatoric surrogate is the one for
    // div(dev(sigma)) alone. The bulk modulus here is a near-incompressibility
    // penalty, two orders above the shear modulus, so including it makes the
    // surrogate stiff enough that the linear solve does not converge at all
    fillScalarTangent
    (
        response,
        (4.0/3.0)*mu_.value() + bulkModulus_.value(),
        (4.0/3.0)*mu_.value()
    );

    // No analytical fourth-order tangent has been derived for this law
    fourthOrderByFiniteDifferenceOnly(kin, inputs, state, response);
}


// ************************************************************************* //
