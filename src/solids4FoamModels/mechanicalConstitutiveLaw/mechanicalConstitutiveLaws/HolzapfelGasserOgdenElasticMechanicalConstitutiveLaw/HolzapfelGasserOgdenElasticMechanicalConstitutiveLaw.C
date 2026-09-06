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

#include "HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw, 0
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw::
HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    mu_(dict.lookup("mu")),
    k1_(dict.lookup("k1")),
    k2_(readScalar(dict.lookup("k2"))),
    // Read in degrees, as the legacy law reads it, and kept in radians
    fibreAngle_(readScalar(dict.lookup("fibreAngle"))*M_PI/180.0),
    bulkModulus_(dict.lookup("bulkModulus")),
    EcDefault_
    (
        // A case whose fibres are the same everywhere may say so here. One
        // that does not is expected to supply Ec and Ea fields, and the
        // default is then zero so that a missing field is a zero-length
        // direction and refused, rather than silently isotropic
        dict.lookupOrDefault<Switch>("uniformLocalBasis", false)
      ? vector(dict.lookup("Ec"))
      : vector::zero
    ),
    EaDefault_
    (
        dict.lookupOrDefault<Switch>("uniformLocalBasis", false)
      ? vector(dict.lookup("Ea"))
      : vector::zero
    )
{
    if (bulkModulus_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "bulkModulus must be positive: it is the penalty that keeps "
            << "this material near incompressible. The legacy law has no "
            << "volumetric term at all and takes the whole spherical stress "
            << "from the solid model's pressure; here the pressure replaces "
            << "this response instead, so there has to be one."
            << exit(FatalIOError);
    }

    if (mu_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "mu must be positive." << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // The local basis. Prescribed rather than history: the case supplies it,
    // this law only reads it, and it does not change as the material deforms -
    // these are directions in the reference configuration
    spec.addVector
    (
        "Ec",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        EcDefault_
    );

    spec.addVector
    (
        "Ea",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        EaDefault_
    );
}


void Foam::HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw::evaluate
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
    // accessor is selected. A prescribed field is never written, so its two
    // times always hold the same value; taking the old one is what makes this
    // work inside a shadow state, which aliases its parent's old-time fields
    // and owns current-time fields that start empty
    const mechanicalConstitutiveLawState& cState = state;
    const Field<vector>& Ec = cState.vectorField0("Ec");
    const Field<vector>& Ea = cState.vectorField0("Ea");

    // Whether the caller wants the isochoric stress and the volumetric
    // response separately, as a mixed displacement-pressure formulation does
    const bool wantsSplit = response.wantsVolumetricSplit();
    UIndirectList<scalar>* volumetricPtr =
        wantsSplit ? &response.volumetric() : nullptr;

    const scalar muVal = mu_.value();
    const scalar k1Val = k1_.value();
    const scalar bulkVal = bulkModulus_.value();

    const scalar cosA = Foam::cos(fibreAngle_);
    const scalar sinA = Foam::sin(fibreAngle_);

    // The stress divides by J and the isochoric split forms J^(-1/3), so a J
    // this close to zero is already meaningless numerically
    const scalar Jmin = sqrt(SMALL);

    forAll(sigma, i)
    {
        const scalar Ji = J[i];

        if (Ji <= Jmin)
        {
            FatalErrorInFunction
                << "Deformation gradient determinant J = " << Ji
                << " at index " << i << " is not usable: J must be positive "
                << "and greater than " << Jmin << ", since the stress divides "
                << "by J."
                << exit(FatalError);
        }

        const scalar magEc = mag(Ec[i]);
        const scalar magEa = mag(Ea[i]);

        if (magEc < SMALL || magEa < SMALL)
        {
            FatalErrorInFunction
                << "The local basis has zero length at index " << i << '.'
                << nl << nl
                << "    Either supply Ec and Ea fields - calcLocCoordinates "
                << "writes them - or set 'uniformLocalBasis yes' and give "
                << "'Ec' and 'Ea' in this law's dictionary."
                << exit(FatalError);
        }

        const vector ec(Ec[i]/magEc);
        const vector ea(Ea[i]/magEa);

        // The two fibre families, wound symmetrically about the
        // circumferential direction
        const vector N4(cosA*ec + sinA*ea);
        const vector N6(-cosA*ec + sinA*ea);

        // The isochoric deformation. Every invariant below is taken on this,
        // which is what makes the stress a function of shape alone and lets
        // the volumetric response be replaced by a solved pressure
        const tensor Fbar(pow(Ji, -1.0/3.0)*F[i]);

        const symmTensor bbar(symm(Fbar & Fbar.T()));
        const symmTensor Cbar(symm(Fbar.T() & Fbar));

        const scalar I4 = Cbar && symm(N4*N4);
        const scalar I6 = Cbar && symm(N6*N6);

        const vector n4(Fbar & N4);
        const vector n6(Fbar & N6);

        // Fbar*Sbar*Fbar^T for each term of the energy: the ground substance
        // contributes mu*bbar, and each fibre family contributes its stress
        // along the deformed fibre direction
        symmTensor T(muVal*bbar);

        T += (2.0*k1Val*(I4 - 1.0)*exp(k2_*sqr(I4 - 1.0)))*symm(n4*n4);
        T += (2.0*k1Val*(I6 - 1.0)*exp(k2_*sqr(I6 - 1.0)))*symm(n6*n6);

        // The isochoric Cauchy stress
        const symmTensor s(dev(T)/Ji);

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

        // A preconditioner, not a tangent of this energy. The fibre stiffness
        // is left out of it deliberately: it varies by orders of magnitude
        // with the fibre stretch, and a Laplacian coefficient that swung with
        // it would help the solve less than a steady one
        scalar Keff = 0.0;

        switch (response.tangentReq())
        {
            case tangentRequest::scalar:
                Keff = (4.0/3.0)*muVal + bulkVal;
                break;

            case tangentRequest::scalarDeviatoric:
                // A mixed displacement-pressure formulation carries the
                // volumetric response in its own equation, so the surrogate
                // here is the one for div(dev(sigma)) alone. Including the
                // bulk modulus - a near-incompressibility penalty orders
                // above the shear modulus - makes it stiff enough that the
                // linear solve does not converge
                Keff = (4.0/3.0)*muVal;
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
