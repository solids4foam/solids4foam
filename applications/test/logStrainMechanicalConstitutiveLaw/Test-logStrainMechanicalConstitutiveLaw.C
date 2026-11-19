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

Application
    Test-logStrainMechanicalConstitutiveLaw

Description
    Material-point benchmarks of the logStrain law, which drives a
    small-strain law with the elastic logarithmic strain.

    Each law is driven through the manager at the one integration point of a
    one-cell mesh, over a deformation path of many steps, so that its history
    is carried from step to step as in a solid model. The following are
    checked:
      1. The tensor logarithm and exponential: closed-form values, also for
         repeated eigenvalues, and exp(log(A)) = A.
      2. With a linearElastic sub-law the law is Hencky's hyperelastic law,
         tau = 2 mu dev(ln V) + K tr(ln V) I, at the end of a path that
         stretches and rotates: elastic, so path independent.
      3. On a path whose principal directions stay fixed, the Kirchhoff
         stress of logStrain over linearElasticMisesPlastic equals the
         stress of linearElasticMisesPlastic itself driven by the
         logarithmic strain, step by step, past yield and through a change
         of loading direction. This is exact for the multiplicative split.
      4. A rigid rotation superposed after yielding rotates the stress and
         does nothing else.
      5. Simple shear to a shear of 5, with perfect plasticity: the shear
         stress approaches the yield stress in shear, sigmaY/sqrt(3), and
         agrees with neoHookeanElasticMisesPlastic, an independent
         multiplicative J2 law that differs only in its elastic response,
         which is negligible at metal yield strains. An additive split of
         the Lagrangian logarithmic strain falls well below both. The curve
         is written to simpleShear.csv.
      6. A sub-law that cannot be driven by a trial elastic strain is
         refused.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalConstitutiveLawManager.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "logStrainMechanicalConstitutiveLaw.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

label nPass = 0;
label nFail = 0;

void check
(
    const bool ok,
    const string& what,
    const scalar value,
    const scalar tol
)
{
    if (ok)
    {
        ++nPass;
        Info<< "PASS: " << what.c_str() << " (" << value << " <= " << tol
            << ")" << endl;
    }
    else
    {
        ++nFail;
        Info<< "FAIL: " << what.c_str() << " (" << value << " > " << tol
            << ")" << endl;
    }
}


//- Rotation by angle about axis (Rodrigues)
tensor rotation(const vector& axis, const scalar angle)
{
    const vector n(axis/mag(axis));
    const tensor K
    (
        0,     -n.z(),  n.y(),
        n.z(),  0,     -n.x(),
       -n.y(),  n.x(),  0
    );

    return I + Foam::sin(angle)*K + (1 - Foam::cos(angle))*(K & K);
}


//- Q & diag(d) & Q^T, written out rather than taken from the code under test
symmTensor rotatedDiag(const tensor& Q, const vector& d)
{
    const tensor D(d.x(), 0, 0, 0, d.y(), 0, 0, 0, d.z());

    return symm(Q & D & Q.T());
}


//- A law, driven at the one integration point of the cell-centred topology
class pointDriver
{
    const dictionary mechanicalProperties_;
    mechanicalConstitutiveLawManager manager_;
    const integrationPointTopology& topo_;
    const scalar dt_;
    tensor F0_;
    tensor gradD0_;

    static dictionary makeDict(const word& name, const dictionary& lawDict)
    {
        OStringStream os;
        os  << "planeStress no;" << nl
            << "mechanical" << nl << "(" << nl
            << name << nl << lawDict << nl
            << ");" << nl;

        IStringStream is(os.str());

        return dictionary(is);
    }

public:

    pointDriver
    (
        const fvMesh& mesh,
        const word& name,
        const dictionary& lawDict
    )
    :
        mechanicalProperties_(makeDict(name, lawDict)),
        manager_(mesh, mechanicalProperties_),
        topo_
        (
            manager_.topologyFor(cellCentredIntegrationPointTopology::typeName)
        ),
        dt_(mesh.time().deltaTValue()),
        F0_(I),
        gradD0_(tensor::zero)
    {}

    //- The Cauchy stress at F, at the end of a step from the F of the last
    //  call. The caller advances the time first, so that the history of
    //  the last step is committed
    symmTensor finiteStrainStep(const tensor& F)
    {
        const label n = topo_.nIntegrationPoints();

        const tensorField Ff(n, F);
        const tensorField F0f(n, F0_);
        const tensorField Finv(n, inv(F));
        const tensorField Finv0(n, inv(F0_));
        const scalarField J(n, det(F));
        const scalarField J0(n, det(F0_));
        symmTensorField sigma(n, symmTensor::zero);

        manager_.updateStressFiniteStrain
        (
            topo_, Ff, F0f, Finv, Finv0, J, J0, dt_, sigma,
            nullptr, nullptr, tangentRequest::none
        );

        F0_ = F;

        return sigma[0];
    }

    //- As above, for a small-strain law and a displacement gradient
    symmTensor smallStrainStep(const tensor& gradD)
    {
        const label n = topo_.nIntegrationPoints();

        const tensorField g(n, gradD);
        const tensorField g0(n, gradD0_);
        symmTensorField sigma(n, symmTensor::zero);

        manager_.updateStressSmallStrain
        (
            topo_, g, g0, dt_, sigma, nullptr, nullptr, tangentRequest::none
        );

        gradD0_ = gradD;

        return sigma[0];
    }
};


//- The dictionary of a law, from a string
dictionary lawDict(const string& s)
{
    IStringStream is(s);

    return dictionary(is);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    typedef logStrainMechanicalConstitutiveLaw logLaw;

    // Elastic constants shared by every law here
    const scalar E = 200e9;
    const scalar nu = 0.3;
    const scalar mu = E/(2*(1 + nu));
    const scalar K = E/(3*(1 - 2*nu));
    const scalar sigmaY = 400e6;

    const string elastic
    (
        "rho rho [1 -3 0 0 0 0 0] 7800;"
        "E E [1 -1 -2 0 0 0 0] 200e9;"
        "nu nu [0 0 0 0 0 0 0] 0.3;"
    );

    // Linear hardening, H = 1 GPa, which the laws return in closed form
    const string hardening
    (
        elastic
      + "file \"$FOAM_CASE/constant/linearHardening\"; outOfBounds clamp;"
    );

    // Perfect plasticity at sigmaY
    const string perfect
    (
        elastic
      + "file \"$FOAM_CASE/constant/perfectPlasticity\"; outOfBounds clamp;"
    );

    // 1. The tensor logarithm and exponential
    {
        Info<< nl << "1. Tensor logarithm and exponential" << endl;

        const tensor Q(rotation(vector(1, -2, 0.5), 0.7));

        scalar err = mag(logLaw::logm(symmTensor(I)));

        check(err <= 1e-14, "log(I) = 0", err, 1e-14);

        // Distinct, repeated, and nearly repeated eigenvalues
        List<vector> eigs(3);
        eigs[0] = vector(2.5, 0.4, 1.3);
        eigs[1] = vector(2.0, 2.0, 0.5);
        eigs[2] = vector(1.0, 1.0 + 1e-10, 1.0 - 1e-10);

        forAll(eigs, i)
        {
            const vector& d = eigs[i];

            const symmTensor A(rotatedDiag(Q, d));
            const symmTensor logA
            (
                rotatedDiag
                (
                    Q,
                    vector(Foam::log(d.x()), Foam::log(d.y()), Foam::log(d.z()))
                )
            );

            err = mag(logLaw::logm(A) - logA)/max(mag(logA), 1.0);

            check
            (
                err <= 1e-12,
                "log of Q diag(d) Q^T, eigenvalue set " + Foam::name(i),
                err, 1e-12
            );

            err = mag(logLaw::expm(logLaw::logm(A)) - A)/mag(A);

            check
            (
                err <= 1e-12,
                "exp(log(A)) = A, eigenvalue set " + Foam::name(i),
                err, 1e-12
            );
        }
    }

    // 2. Hencky's law, at the end of a path that stretches and rotates
    {
        Info<< nl << "2. Hencky elasticity" << endl;

        pointDriver law
        (
            mesh,
            "hencky",
            lawDict
            (
                "type logStrain; mechanicalLaw { type linearElastic; "
              + elastic + " }"
            )
        );

        const tensor G(0.3, 0.2, -0.1, 0.05, -0.2, 0.15, -0.1, 0.1, 0.25);
        const label nSteps = 20;

        tensor F(I);
        symmTensor sigma(symmTensor::zero);

        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            ++runTime;

            const scalar t = scalar(stepI)/nSteps;
            F = rotation(vector(0, 1, 1), 1.2*t) & (I + t*G);
            sigma = law.finiteStrainStep(F);
        }

        const scalar J = det(F);
        const symmTensor lnV(0.5*logLaw::logm(symm(F & F.T())));
        const symmTensor expected
        (
            (2*mu*dev(lnV) + K*tr(lnV)*symmTensor(I))/J
        );

        const scalar err = mag(sigma - expected)/mag(expected);

        check
        (
            err <= 1e-10,
            "logStrain over linearElastic is Hencky's law, after "
          + Foam::name(nSteps) + " steps with rotation",
            err, 1e-10
        );
    }

    // 3. Fixed principal directions: the small-strain law in log strain
    {
        Info<< nl << "3. Coaxial path, against the small-strain law" << endl;

        pointDriver logStrain
        (
            mesh,
            "logStrainPlastic",
            lawDict
            (
                "type logStrain; mechanicalLaw { type "
                "linearElasticMisesPlastic; " + hardening + " }"
            )
        );

        pointDriver smallStrain
        (
            mesh,
            "smallStrainPlastic",
            lawDict("type linearElasticMisesPlastic; " + hardening)
        );

        // Uniaxial stretch to 30%, then a lateral stretch and compression
        // with the axial stretch held: the second leg turns the loading
        // direction, so the return mapping is not proportional
        scalar maxErr = 0;
        scalar maxSigmaEq = 0;

        for (label stepI = 1; stepI <= 60; ++stepI)
        {
            ++runTime;

            vector lambda(1 + 0.01*stepI, 1, 1);

            if (stepI > 30)
            {
                lambda = vector
                (
                    1.3, 1 + 0.005*(stepI - 30), 1 - 0.003*(stepI - 30)
                );
            }

            const tensor F
            (
                lambda.x(), 0, 0, 0, lambda.y(), 0, 0, 0, lambda.z()
            );

            const tensor lnF
            (
                Foam::log(lambda.x()), 0, 0,
                0, Foam::log(lambda.y()), 0,
                0, 0, Foam::log(lambda.z())
            );

            const symmTensor tau(det(F)*logStrain.finiteStrainStep(F));
            const symmTensor sigmaSmall(smallStrain.smallStrainStep(lnF));

            maxErr = max(maxErr, mag(tau - sigmaSmall)/mag(sigmaSmall));
            maxSigmaEq = max(maxSigmaEq, Foam::sqrt(1.5)*mag(dev(tau)));
        }

        check
        (
            maxSigmaEq > sigmaY,
            "the coaxial path yields, sigmaY / max equivalent stress",
            sigmaY/maxSigmaEq, 1
        );

        check
        (
            maxErr <= 1e-9,
            "Kirchhoff stress equals the small-strain law's at the "
            "logarithmic strain, at every step",
            maxErr, 1e-9
        );
    }

    // 4. A rigid rotation after yielding rotates the stress
    {
        Info<< nl << "4. Rigid rotation after yield" << endl;

        pointDriver law
        (
            mesh,
            "rotatedPlastic",
            lawDict
            (
                "type logStrain; mechanicalLaw { type "
                "linearElasticMisesPlastic; " + hardening + " }"
            )
        );

        const tensor F1(1.2, 0.05, 0, 0, 0.95, 0, 0, 0, 0.95);

        symmTensor sigma1(symmTensor::zero);

        for (label stepI = 1; stepI <= 10; ++stepI)
        {
            ++runTime;

            sigma1 = law.finiteStrainStep(I + 0.1*stepI*(F1 - I));
        }

        scalar maxErr = 0;

        for (label stepI = 1; stepI <= 18; ++stepI)
        {
            ++runTime;

            const tensor R(rotation(vector(1, 1, 0), 0.1*stepI));
            const symmTensor sigma(law.finiteStrainStep(R & F1));

            maxErr =
                max(maxErr, mag(sigma - transform(R, sigma1))/mag(sigma1));
        }

        const scalar tauEq1 = Foam::sqrt(1.5)*mag(dev(det(F1)*sigma1));

        check
        (
            tauEq1 > sigmaY,
            "the stretch yields, sigmaY / equivalent stress",
            sigmaY/tauEq1, 1
        );

        check
        (
            maxErr <= 1e-9,
            "a superposed rigid rotation only rotates the stress",
            maxErr, 1e-9
        );
    }

    // 5. Simple shear with perfect plasticity
    {
        Info<< nl << "5. Simple shear" << endl;

        pointDriver logStrain
        (
            mesh,
            "logStrainPerfect",
            lawDict
            (
                "type logStrain; mechanicalLaw { type "
                "linearElasticMisesPlastic; " + perfect + " }"
            )
        );

        pointDriver neoHookean
        (
            mesh,
            "neoHookeanPerfect",
            lawDict("type neoHookeanElasticMisesPlastic; " + perfect)
        );

        const scalar gammaMax = 5;
        const label nSteps = 500;
        const scalar tauY = sigmaY/Foam::sqrt(3.0);

        OFstream csv(runTime.path()/"simpleShear.csv");
        csv << "gamma,sigmaXY_logStrain,sigmaXX_logStrain,"
            << "sigmaXY_neoHookeanMisesPlastic,"
            << "sigmaXX_neoHookeanMisesPlastic" << endl;

        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            ++runTime;

            const scalar gamma = gammaMax*stepI/nSteps;
            const tensor F(1, gamma, 0, 0, 1, 0, 0, 0, 1);

            const symmTensor sLog(logStrain.finiteStrainStep(F));
            const symmTensor sNH(neoHookean.finiteStrainStep(F));

            if (stepI % 10 == 0)
            {
                csv << gamma << "," << sLog.xy() << "," << sLog.xx() << ","
                    << sNH.xy() << "," << sNH.xx() << endl;
            }

            // At shears of 1, 2 and 5
            if (stepI % 100 == 0 && stepI != 300 && stepI != 400)
            {
                const scalar diff = mag(sLog.xy() - sNH.xy())/mag(sNH.xy());
                const scalar fromTauY = mag(sLog.xy() - tauY)/tauY;

                Info<< "    gamma = " << gamma
                    << ": sigma_xy logStrain = " << sLog.xy()
                    << ", neoHookeanElasticMisesPlastic = " << sNH.xy()
                    << ", sigmaY/sqrt(3) = " << tauY << endl;

                check
                (
                    diff <= 1e-3,
                    "simple shear stress agrees with "
                    "neoHookeanElasticMisesPlastic at gamma = "
                  + Foam::name(gamma),
                    diff, 1e-3
                );

                check
                (
                    fromTauY <= 1e-3,
                    "simple shear stress is sigmaY/sqrt(3) at gamma = "
                  + Foam::name(gamma),
                    fromTauY, 1e-3
                );
            }
        }
    }

    // 6. A sub-law that cannot be driven by a trial elastic strain
    {
        Info<< nl << "6. Refused sub-law" << endl;

        FatalError.throwExceptions();
        FatalIOError.throwExceptions();

        bool threw = false;

        try
        {
            autoPtr<mechanicalConstitutiveLaw> lawPtr
            (
                mechanicalConstitutiveLaw::New
                (
                    lawDict
                    (
                        "type logStrain; mechanicalLaw { type "
                        "neoHookeanElastic; " + elastic + " }"
                    )
                )
            );
        }
        catch (const Foam::error&)
        {
            threw = true;
        }

        check
        (
            threw,
            "a sub-law that provides no elastic predictor is refused",
            threw ? 0 : 1, 0
        );
    }

    Info<< nl << "logStrain checks: " << nPass << " passed, " << nFail
        << " failed" << endl;

    return (nFail == 0 ? 0 : 1);
}


// ************************************************************************* //
