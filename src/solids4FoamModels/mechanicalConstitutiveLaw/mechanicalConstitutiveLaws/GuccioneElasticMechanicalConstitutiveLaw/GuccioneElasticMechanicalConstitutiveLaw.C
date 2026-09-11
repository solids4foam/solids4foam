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

    // Read at old time: a prescribed field is never written, so its two times
    // always hold the same value, and the old-time one is what a shadow state
    // aliases. A tangent query evaluated into a shadow would find the
    // current-time field empty
    const Field<vector>& f0 = state.getVectorField0("f0");

    // Whether the caller wants the isochoric stress and the volumetric
    // response separately, as a mixed displacement-pressure formulation does
    const bool wantsSplit = response.wantsVolumetricSplit();
    UIndirectList<scalar>* volumetricPtr =
        wantsSplit ? &response.volumetric() : nullptr;

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

        // Green-Lagrange strain of the isochoric deformation, not of the
        // whole of it. Fbar = J^(-1/3)*F carries the shape change and none of
        // the volume change, so Q depends on shape alone and the volumetric
        // response below is the only place volume enters.
        //
        // The legacy law builds Q from the full strain, which makes its energy
        // coupled: its deviatoric stress then varies with J, and a mixed
        // formulation replacing the volumetric part gives a different material
        // rather than the same one solved differently. Written this way the
        // two formulations describe one material, and both reduce to the
        // published model in the incompressible limit it was defined for
        const tensor Fbar(pow(Ji, -1.0/3.0)*F[i]);
        const symmTensor E(0.5*(symm(Fbar.T() & Fbar) - I));
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

        // Push forward through the isochoric deformation and take the
        // deviatoric part:
        //
        //     sigma_iso = dev(Fbar & S & Fbar^T)/J
        //
        // The dev() is not tidying up, and it is not the claim that S is
        // deviatoric. Differentiating W(Cbar) with Cbar = J^(-2/3)*C gives
        //
        //     S_iso = J^(-2/3) DEV[S],   DEV[A] = A - (1/3)*(A && C)*inv(C)
        //
        // and the DEV[] is the chain rule through J^(-2/3), not a projection
        // that the energy already satisfies: S here has a non-zero trace and
        // a non-zero S && C.
        //
        // DEV[] and dev() are different operators - DEV[A] is C-orthogonal,
        // dev(a) is trace-free, and dev(F & S & F.T())/J is *not* the Cauchy
        // stress of DEV[S]. What makes the line below right is that they
        // correspond exactly under push-forward:
        //
        //     (1/J)*F & DEV[A] & F.T() = dev((1/J)*F & A & F.T())
        //
        // because F & inv(C) & F.T() == I identically. So applying dev() in
        // the current configuration, after the push-forward, is applying
        // DEV[] in the reference one. Applying dev() to S instead would not
        // be.
        //
        // Verified numerically: with the line as written, sigma agrees with a
        // central-difference dW/dF of the energy declared in the header to
        // 1e-8; with the dev() removed it is wrong by 16 to 37 per cent over
        // shear, uniaxial and large-strain states
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

    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        // The small-strain estimate the legacy law uses. It is a
        // preconditioner, not a tangent of this energy
        scalar Keff = 0.0;

        switch (response.tangentReq())
        {
            case tangentRequest::scalar:
                Keff = (4.0/3.0)*mu_.value() + bulkModulus_.value();
                break;

            case tangentRequest::scalarDeviatoric:
                // A mixed displacement-pressure formulation carries the
                // volumetric response in its own equation, so the Laplacian
                // surrogate here is the one for div(dev(sigma)) alone. The
                // bulk modulus here is a near-incompressibility penalty, two
                // orders above the shear modulus, so including it makes the
                // surrogate stiff enough that the linear solve does not
                // converge at all
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
}


// ************************************************************************* //
