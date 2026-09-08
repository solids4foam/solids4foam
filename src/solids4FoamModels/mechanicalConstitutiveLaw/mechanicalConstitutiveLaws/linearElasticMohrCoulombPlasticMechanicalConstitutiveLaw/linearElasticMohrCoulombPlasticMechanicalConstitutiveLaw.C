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

#include "linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "dimensionedSymmTensor.H"
#include "diagTensor.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "mathematicalConstants.H"
#else
    #include "mathematicalConstants.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw, 0
    );
    addToRunTimeSelectionTable
    (
        mechanicalConstitutiveLaw,
        linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw,
        mechanicalConstitutiveLaw
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw::calculateEigens
(
    vector& sigma_prin,
    tensor& ev,
    const symmTensor sigma
) const
{
    // Check for a zero stress tensor
    if (mag(sigma) < SMALL)
    {
        ev = I;

        return;
    }

    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if
    (
        (
            mag(sigma.xy()) + mag(sigma.xz()) + mag(sigma.xy())
          + mag(sigma.yz()) + mag(sigma.xz()) + mag(sigma.yz())
        ) < 1e-6*(mag(sigma.xx()) + mag(sigma.yy()) + mag(sigma.zz()))
        )
    {
        // Diagonal matrix
        i = sigma.xx();
        ii = sigma.yy();
        iii = sigma.zz();
    }
    else
    {
        const scalar a = -sigma.xx() - sigma.yy() - sigma.zz();

        const scalar b =
            sigma.xx()*sigma.yy()
          + sigma.xx()*sigma.zz()
          + sigma.yy()*sigma.zz()
          - sigma.xy()*sigma.xy()
          - sigma.xz()*sigma.xz()
          - sigma.yz()*sigma.yz();

        const scalar c =
          - sigma.xx()*sigma.yy()*sigma.zz()
          - sigma.xy()*sigma.yz()*sigma.xz()
          - sigma.xz()*sigma.xy()*sigma.yz()
          + sigma.xz()*sigma.yy()*sigma.xz()
          + sigma.xy()*sigma.xy()*sigma.zz()
          + sigma.xx()*sigma.yz()*sigma.yz();

        // If there is a zero root
        if (mag(c) < SMALL)
        {
            const scalar disc = Foam::max(sqr(a) - 4*b, 0.0);

            WarningIn("poroMohrCoulob::calculateEigens(...)")
                << "Stress tensor has a zero root!" << endl;

            const scalar q = -0.5*Foam::sqrt(max(scalar(0), disc));

            i = 0;
            ii = -0.5*a + q;
            iii = -0.5*a - q;
        }
        else
        {
            const scalar Q = (a*a - 3.0*b)/9.0;
            const scalar R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;

            const scalar R2 = sqr(R);
            const scalar Q3 = pow3(Q);

            // Three different real roots
            if (R2 < Q3)
            {
                const scalar sqrtQ = Foam::sqrt(Q);
                const scalar theta = Foam::acos(R/(Q*sqrtQ));

                const scalar m2SqrtQ = -2*sqrtQ;
                const scalar aBy3 = a/3;

                i = m2SqrtQ*Foam::cos(theta/3) - aBy3;
#ifdef OPENFOAM_NOT_EXTEND
                ii =
                    m2SqrtQ*Foam::cos((theta + constant::mathematical::twoPi)/3.0)
                  - aBy3;
                iii =
                    m2SqrtQ*Foam::cos((theta - constant::mathematical::twoPi)/3.0)
                  - aBy3;
#else
                ii =
                    m2SqrtQ*Foam::cos((theta + mathematicalConstant::twoPi)/3.0)
                  - aBy3;
                iii =
                    m2SqrtQ*Foam::cos((theta - mathematicalConstant::twoPi)/3.0)
                  - aBy3;
#endif
            }
            else
            {
                const scalar A = Foam::cbrt(R + Foam::sqrt(R2 - Q3));

                // Three equal real roots
                if (A < SMALL)
                {
                    const scalar root = -a/3;
                    i = root;
                    ii = root;
                    iii = root;
                }
                else
                {
                    // Complex roots
                    WarningIn
                    (
                        "linearElasticMohrCoulombPlastic::calculateEigens(...)"
                    )   << "Complex eigenvalues detected for symmTensor: "
                        << sigma << nl << "Setting roots to zero!" << endl;

                    i = 0;
                    ii = 0;
                    iii = 0;
                }
            }
        }
    }

    // Sort the eigenvalues into ascending order

    if (mag(i) > mag(ii))
    {
        Swap(i, ii);
    }

    if (mag(ii) > mag(iii))
    {
        Swap(ii, iii);
    }

    if (mag(i) > mag(ii))
    {
        Swap(i, ii);
    }

    sigma_prin[0] = i;
    sigma_prin[1] = ii;
    sigma_prin[2] = iii;

    for (int j = 0; j < 3; j++)
    {
        if (mag(sigma_prin[j]) < SMALL)
        {
            if (j == 0)
            {
                ev[j*3 + 0] = 1;
                ev[j*3 + 1] = 0;
                ev[j*3 + 2] = 0;
            }
            else if (j == 1)
            {
                ev[j*3 + 0] = 0;
                ev[j*3 + 1] = 1;
                ev[j*3 + 2] = 0;
            }
            else
            {
                ev[j*3 + 0] = 0;
                ev[j*3 + 1] = 0;
                ev[j*3 + 2] = 1;
            }
        }
        else
        {
            const symmTensor A = symmTensor(sigma - sigma_prin[j]*I);

            // Calculate the sub-determinants of the 3 components
            const scalar sd0 = A.yy()*A.zz() - A.yz()*A.yz();
            const scalar sd1 = A.xx()*A.zz() - A.xz()*A.xz();
            const scalar sd2 = A.xx()*A.yy() - A.xy()*A.xy();

            const scalar magSd0 = mag(sd0);
            const scalar magSd1 = mag(sd1);
            const scalar magSd2 = mag(sd2);

            // Evaluate the eigenvector using the largest sub-determinant
            if ((magSd0 > magSd1) && (magSd0 > magSd2) && (magSd0 > SMALL))
            {
                vector newEv =
                    vector
                    (
                        1,
                        (A.yz()*A.xz() - A.zz()*A.xy())/sd0,
                        (A.yz()*A.xy() - A.yy()*A.xz())/sd0
                    );
                newEv /= mag(newEv);

                ev[j*3] = newEv[0];
                ev[j*3+1] = newEv[1];
                ev[j*3+2] = newEv[2];
            }
            else if ((magSd1 > magSd2) && (magSd1 > SMALL))
            {
                vector newEv =
                    vector
                    (
                        (A.xz()*A.yz() - A.zz()*A.xy())/sd1,
                        1,
                        (A.xz()*A.xy() - A.xx()*A.yz())/sd1
                    );
                newEv /= mag(newEv);

                ev[j*3] = newEv[0];
                ev[j*3 + 1] = newEv[1];
                ev[j*3 + 2] = newEv[2];
            }
            else if (magSd2 > SMALL)
            {
                vector newEv =
                    vector
                    (
                        (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
                        (A.xy()*A.xz() - A.xx()*A.yz())/sd2,
                        1
                    );
                newEv /= mag(newEv);

                ev[j*3] = newEv[0];
                ev[j*3 + 1] = newEv[1];
                ev[j*3 + 2] = newEv[2];
            }
            else
            {
                WarningIn
                (
                    "linearElasticMohrCoulombPlastic::calculateEigens(...)"
                )   << "Strange things detected for stress tensor: "
                    << sigma << nl
                    << "Setting eigenvectors to (0, 0, 0)" << endl;

                ev[j*3] = 0;
                ev[j*3 + 1] = 0;
                ev[j*3 + 2] = 0;
            }
        }
    }
}


void Foam::linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw::calculateStress
(
    symmTensor& sigma,
    scalar& activeYield
) const
{
    // Principal stresses
    vector sigma_prin = vector::zero;

    // Principal directions
    tensor ev = tensor::zero;

    activeYield = 0;

    for (int i = 0; i < 6; i++)
    {
        if (sigma[i] > SMALL)
        {
            sigma[i] += (i + 1)*(1e-6)*sigma[i];
        }
    }

    calculateEigens(sigma_prin, ev, sigma);

    // Re-order the principal stresses
    scalar sigma1 = sigma_prin[0];
    scalar sigma2 = sigma_prin[1];
    scalar sigma3 = sigma_prin[2];

    label sigma1_po = 0;
    label sigma2_po = 1;
    label sigma3_po = 2;

    for (int i = 1; i < 3; i++)
    {
        if (sigma_prin[i] > sigma1)
        {
            sigma1 = sigma_prin[i];
            sigma1_po = i;
        }
    }

    for (int i = 1; i >= 0; i--)
    {
        if (sigma_prin[i] < sigma3)
        {
            sigma3 = sigma_prin[i];
            sigma3_po = i;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if ((i != sigma1_po) && (i != sigma3_po))
        {
            sigma2_po = i;
            sigma2 = sigma_prin[i];
        }
    }

    vector sigmaB = vector(sigma1, sigma2, sigma3);

    // Evaluate the yielding function f
    const scalar f = k_*sigmaB[0] - sigmaB[2] - 2*c_.value()*Foam::sqrt(k_);

    // Check if yielding
    if (f > SMALL)
    {
        // Determine the return type

        // Calculate the boundary planes
        const scalar t1 =
            (r_lg1_ & (invC_ & (sigma_prin - sigma_a_)))
            /(r_lg1_ & (invC_ & r_lf1_));

        const scalar t2 =
            (r_lg2_ & (invC_ & (sigma_prin - sigma_a_)))
            /(r_lg2_ & (invC_ & r_lf2_));

        const scalar p_12 = (rp_^r_lf1_) & (sigma_prin - sigma_a_);

        const scalar p_13 = (rp_^r_lf2_) & (sigma_prin - sigma_a_);

        // Compute the stress field based on the correct return type

        if ((p_12 >= 0) && (p_13 <= 0))
        {
            sigmaB -= f*rp_;
        }
        else if ((p_12 < 0) && (p_13 < 0))
        {
            sigmaB = t1*r_lf1_ + sigma_a_;
        }
        else if ((p_12 > 0) && (p_13 > 0))
        {
            sigmaB = t2*r_lf2_ + sigma_a_;
        }
        else if ((t1 > 0) && (t2 > 0))
        {
            sigmaB = sigma_a_;
        }

        sigma_prin[sigma1_po] = sigmaB[0];
        sigma_prin[sigma2_po] = sigmaB[1];
        sigma_prin[sigma3_po] = sigmaB[2];

        // Form the diagonal stress
        const diagTensor sigma_diag =
            diagTensor(sigma_prin[0], sigma_prin[1], sigma_prin[2]);

        // Transform the returned stress back into general space
        sigma = symm(ev.T() & sigma_diag & ev);

        activeYield = 1.0;
    }
    else
    {
        activeYield = 0.0;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw::
linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw
(
    const dictionary& dict
)
:
    mechanicalConstitutiveLaw(dict),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    lambda_("lambda", E_.dimensions(), 0.0),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_("K", E_.dimensions(), 0.0),
    varPhi_(dict.lookup("frictionAngle")),
    c_(dict.lookup("cohesion")),
    varPsi_(dict.lookup("dilationAngle")),
    k_(0.0),
    m_(0.0),
    a_(vector::zero),
    b_(vector::zero),
    C_(symmTensor::zero),
    invC_(symmTensor::zero),
    rp_(vector::zero),
    r_lf1_(vector::zero),
    r_lf2_(vector::zero),
    r_lg1_(vector::zero),
    r_lg2_(vector::zero),
    sigma_a_(vector::zero)
{
    // planeStress is injected into this dictionary by the manager from the
    // top-level entry in mechanicalProperties; it is not given by the user in
    // this sub-dictionary
    const Switch planeStress(dict.lookupOrDefault<Switch>("planeStress", false));

    lambda_ =
        planeStress
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_));

    K_ = lambda_ + (2.0/3.0)*mu_;

#ifdef OPENFOAM_NOT_EXTEND
    const scalar piBy180 = constant::mathematical::pi/180.0;
#else
    const scalar piBy180 = mathematicalConstant::pi/180.0;
#endif

    k_ =
    (
        (1 + sin(varPhi_*piBy180))/(1 - sin(varPhi_*piBy180))
    ).value();

    m_ =
    (
        (1 + sin(varPsi_*piBy180))/(1 - sin(varPsi_*piBy180))
    ).value();

    a_ = vector(k_, 0, -1);
    b_ = vector(m_, 0, -1);

    C_ = symmTensor
    (
        (K_ + (4.0/3.0)*mu_).value(),
        (K_ - (2.0/3.0)*mu_).value(),
        (K_ - (2.0/3.0)*mu_).value(),
        (K_ + (4.0/3.0)*mu_).value(),
        (K_ - (2.0/3.0)*mu_).value(),
        (K_ + (4.0/3.0)*mu_).value()
    );

    invC_ = inv(C_);
    rp_ = (C_ & b_)/(b_ & (C_ & a_));
    r_lf1_ = vector(1, 1, k_);
    r_lf2_ = vector(1, k_, k_);
    r_lg1_ = vector(1, 1, m_);
    r_lg2_ = vector(1, m_, m_);
    sigma_a_ = (2.0*c_*Foam::sqrt(k_)/(k_ - 1.0)).value()*vector::one;

    Info<< "    Mohr-Coulomb: frictionAngle " << varPhi_.value()
        << ", dilationAngle " << varPsi_.value()
        << ", cohesion " << c_.value() << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw::
declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    // The stress carried as a variation from the initial stress. The trial
    // stress is built on its previous value, so it is history
    spec.addSymmTensor
    (
        "deltaSigma",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        symmTensor::zero
    );

    // Whether the point is yielding. A diagnostic, but carried so that it
    // survives a step like the rest of the history
    spec.addScalar
    (
        "activeYield",
        mechanicalConstitutiveLawStateSpec::stateRole::persistent,
        0.0
    );

    // The initial stress, as for linearElastic
    spec.addSymmTensor
    (
        "sigma0",
        mechanicalConstitutiveLawStateSpec::stateRole::prescribed,
        symmTensor::zero
    );
}


void Foam::linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw::evaluate
(
    const smallStrainMechanicalConstitutiveLawKinematics& kin,
    const mechanicalConstitutiveLawInputs& inputs,
    mechanicalConstitutiveLawState& state,
    mechanicalConstitutiveLawResponse& response
) const
{
    UIndirectList<symmTensor>& sigma = response.stress();

    const UIndirectList<tensor>& gradD = kin.gradD();
    const UIndirectList<tensor>& gradD0 = kin.gradD0();

    // The previous step's stress variation, which the trial stress is built
    // on, and this step's, which is written here
    const mechanicalConstitutiveLawState& cState = state;
    const Field<symmTensor>& deltaSigma0 =
        cState.symmTensorField0("deltaSigma");
    Field<symmTensor>& deltaSigma = state.symmTensorField("deltaSigma");

    Field<scalar>& activeYield = state.scalarField("activeYield");

    // Read at old time and through a const reference, so that a tangent query
    // evaluated into a shadow state sees the value rather than a silently
    // zero field. See linearElastic for the same reasoning
    const Field<symmTensor>& sigma0 = cState.symmTensorField0("sigma0");

    const scalar muVal = mu_.value();
    const scalar lambdaVal = lambda_.value();

    forAll(sigma, i)
    {
        // Strain increment over the step
        const symmTensor DEpsilon(symm(gradD[i]) - symm(gradD0[i]));

        // Trial elastic stress, from Hooke's law
        const symmTensor DSigmaTrial
        (
            2.0*muVal*DEpsilon + lambdaVal*I*tr(DEpsilon)
        );

        // The trial stress is built on the previous step's variation, not on
        // whatever the caller's storage held. That is what lets this law sit
        // under poroMechanicalLaw, where the incoming stress is effective and
        // the pore pressure is subtracted afterwards
        sigma[i] = deltaSigma0[i] + DSigmaTrial + sigma0[i];

        calculateStress(sigma[i], activeYield[i]);

        deltaSigma[i] = sigma[i] - sigma0[i];
    }

    // Scalar tangent: only if explicitly requested
    if (response.wantsScalarTangent())
    {
        UIndirectList<scalar>& K = response.scalarTangent();

        const scalar Keff = 2.0*muVal + lambdaVal;

        forAll(K, i)
        {
            K[i] = Keff;
        }
    }
}


// ************************************************************************* //
