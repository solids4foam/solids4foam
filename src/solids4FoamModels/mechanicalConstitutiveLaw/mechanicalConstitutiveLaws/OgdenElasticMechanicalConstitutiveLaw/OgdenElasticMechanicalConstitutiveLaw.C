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

#include "OgdenElasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "eig3.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OgdenElasticMechanicalConstitutiveLaw, 0);
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        OgdenElasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OgdenElasticMechanicalConstitutiveLaw::
OgdenElasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    mu1_(dict.lookup("mu1")),
    mu2_(dict.lookup("mu2")),
    mu3_(dict.lookup("mu3")),
    alpha1_(dict.lookup("alpha1")),
    alpha2_(dict.lookup("alpha2")),
    alpha3_(dict.lookup("alpha3")),
    K_(dict.lookup("K")),
    mu_(mu1_ + mu2_ + mu3_)
{
    Info<< "    Ogden: mu = (" << mu1_.value() << " " << mu2_.value()
        << " " << mu3_.value() << "), alpha = (" << alpha1_.value()
        << " " << alpha2_.value() << " " << alpha3_.value()
        << "), K = " << K_.value() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::OgdenElasticMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    spec.addSymmTensor
    (
        "sigma0",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        symmTensor::zero
    );
}


void Foam::OgdenElasticMechanicalConstitutiveLaw::evaluate
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

    // Read at old time and through a const reference, so that a tangent query
    // evaluated into a shadow state sees the value rather than a silently
    // zero field. See linearElastic for the same reasoning
    const mechanicalConstitutiveLawState& cState = state;
    const Field<symmTensor>& sigma0 = cState.symmTensorField0("sigma0");

    const scalar mu1 = mu1_.value();
    const scalar mu2 = mu2_.value();
    const scalar mu3 = mu3_.value();
    const scalar alpha1 = alpha1_.value();
    const scalar alpha2 = alpha2_.value();
    const scalar alpha3 = alpha3_.value();
    const scalar KVal = K_.value();
    const scalar muVal = mu_.value();

    forAll(sigma, i)
    {
        const tensor& Fi = F[i];
        const scalar Ji = J[i];

        // Right Cauchy-Green tensor
        const symmTensor C(symm(Fi.T() & Fi));

        // Its eigenvalues are the squared principal stretches, and its
        // eigenvectors are stored in the rows
        tensor eigVec(tensor::zero);
        vector lambdaSqr(vector::zero);
        eig3().eigen_decomposition(C, eigVec, lambdaSqr);

        // Principal stresses from the principal stretches
        const scalar l1 = max(sqrt(lambdaSqr.x()), VSMALL);
        const scalar l2 = max(sqrt(lambdaSqr.y()), VSMALL);
        const scalar l3 = max(sqrt(lambdaSqr.z()), VSMALL);

        const scalar p1 =
            mu1*pow(l1, alpha1) + mu2*pow(l1, alpha2) + mu3*pow(l1, alpha3);
        const scalar p2 =
            mu1*pow(l2, alpha1) + mu2*pow(l2, alpha2) + mu3*pow(l2, alpha3);
        const scalar p3 =
            mu1*pow(l3, alpha1) + mu2*pow(l3, alpha2) + mu3*pow(l3, alpha3);

        const symmTensor prinStress(p1, 0, 0, p2, 0, p3);

        // Back to the global frame
        const symmTensor s(transform(eigVec.T(), prinStress));

        // The volumetric term, as the legacy law forms it
        const scalar sigmaHyd = 0.5*KVal*(sqr(Ji) - 1.0);

        sigma[i] =
            (1.0/Ji)
           *(
                dev(s - muVal*symmTensor(I))
              + sigmaHyd*symmTensor(I)
              + symm(Fi & sigma0[i] & Fi.T())
            );
    }

    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        const scalar Keff = (4.0/3.0)*muVal + KVal;

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }
}


// ************************************************************************* //
