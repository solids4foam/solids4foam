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

#include "GuccioneElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GuccioneElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        GuccioneElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GuccioneElasticMechanicalConstitutiveLaw::
GuccioneElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    k_(dict.lookup("k")),
    cf_(readScalar(dict.lookup("cf"))),
    ct_(readScalar(dict.lookup("ct"))),
    cfs_(readScalar(dict.lookup("cfs"))),
    bulkModulus_(dict.lookup("bulkModulus")),
    mu_(0.5*k_*(cf_ + cfs_ + ct_)/3.0),
    f0Default_
    (
        // A case whose fibres are the same everywhere may say so here. One
        // that does not is expected to supply an f0 field, and the default is
        // then a zero vector so that a missing field is a zero-length fibre
        // and refused, rather than silently making the material isotropic
        dict.lookupOrDefault<Switch>("uniformFibreField", false)
      ? vector(dict.lookup("f0"))
      : vector::zero
    )
{
    if (bulkModulus_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "bulkModulus must be positive: it is the penalty that keeps "
            << "this material near incompressible." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GuccioneElasticMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // The fibre direction. Prescribed rather than history: the case supplies
    // it, this law only reads it, and it does not change as the material
    // deforms - f0 is a direction in the reference configuration
    spec.addVector
    (
        "f0",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        f0Default_
    );
}


void Foam::GuccioneElasticMechanicalConstitutiveLaw::evaluate
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

    // Read at old time and through a const reference, so that the const
    // accessor is the one selected. A prescribed field is never written, so
    // its two times always hold the same value; taking the old one is what
    // makes this work inside a shadow state, which aliases its parent's
    // old-time fields but owns current-time fields that start empty
    const mechanicalConstitutiveLawState& cState = state;
    const Field<vector>& f0 = cState.vectorField0("f0");

    const scalar kVal = k_.value();
    const scalar bulkVal = bulkModulus_.value();

    // Grouped as the legacy law groups them, so that the two agree term for
    // term rather than only in exact arithmetic
    const scalar cI4 = cf_ - 2.0*cfs_ + ct_;
    const scalar cI5 = cfs_ - ct_;

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

        if (magF0 < SMALL)
        {
            FatalErrorInFunction
                << "The fibre direction has zero length at index " << i
                << '.' << nl
                << "This law needs a direction to be anisotropic about. "
                << "Either supply an f0 field - setFibreField writes one - or "
                << "set 'uniformFibreField yes' and give 'f0' in this "
                << "material's dictionary."
                << exit(FatalError);
        }

        // The structure tensor, from the normalised fibre direction. It is
        // normalised here rather than demanded normalised of the case
        const symmTensor f0f0(sqr(f0[i]/magF0));

        // Green-Lagrange strain
        const symmTensor E(0.5*(symm(F[i].T() & F[i]) - I));
        const symmTensor sqrE(symm(E & E));

        // Invariants: the first two of E, and the two formed with the fibre
        const scalar I1 = tr(E);
        const scalar I2 = 0.5*(sqr(I1) - tr(sqrE));
        const scalar I4 = E && f0f0;
        const scalar I5 = sqrE && f0f0;

        const scalar Q =
            ct_*sqr(I1) - 2.0*ct_*I2 + cI4*sqr(I4) + 2.0*cI5*I5;

        const symmTensor dQdE
        (
            2.0*ct_*E
          + 2.0*cI4*I4*f0f0
          + 2.0*cI5*symm((E & f0f0) + (f0f0 & E))
        );

        // Second Piola-Kirchhoff stress, without the volumetric term
        const symmTensor S(dQdE*0.5*kVal*exp(Q));

        // Push forward to Cauchy, and keep the deviatoric part. The volumetric
        // response is the penalty term below rather than anything in Q
        const symmTensor s(dev(symm(F[i] & S & F[i].T()))/Ji);

        const scalar sigmaHyd = 0.5*bulkVal*(sqr(Ji) - 1.0)/Ji;

        sigma[i] = s + sigmaHyd*I;
    }

    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        // The small-strain estimate the legacy law uses. It is a
        // preconditioner, not a tangent of this energy
        const scalar Keff = (4.0/3.0)*mu_.value() + bulkModulus_.value();

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }
}


// ************************************************************************* //
