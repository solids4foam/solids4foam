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
    Test-LemaitreDamageMechanicalConstitutiveLaw

Description
    Material-point checks of linearElasticMisesPlasticLemaitreDamage.

    The law is driven directly, at one point, over strain paths of many
    steps, so that its damage and plastic history can be read back. The
    following are checked:
      1. With the damage threshold out of reach, the law is
         linearElasticMisesPlastic, step for step.
      2. On a damaging path, for each lemaitreType, the state the law
         returns satisfies the discrete equations, recomputed here from the
         stress: the yield condition sigmaEq/(1 - D) = sigmaY(r), and the
         damage update D - D0 = dr/(1 - D) (-Y/S)^b, with Y written in the
         triaxiality form of Whelan et al. (2025), eq. 52, rather than the
         law's.
      3. Crack closure does nothing when every principal stress is tensile,
         and reduces the damage when one is compressive.
      4. Under uniaxial compression of the strain the triaxiality falls
         below -1/3, and the cut-off keeps the damage at zero.
      5. A point whose damage would pass Dc fails: its damage is 0.99 from
         then on, and its stress falls accordingly.
      6. The damage converges, at first order, as the step is refined.
      7. Under logStrain, on a path whose principal directions stay fixed,
         the Kirchhoff stress equals the small-strain law's at the
         logarithmic strain.
      8. Inconsistent entries are refused.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawManager.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "IStringStream.H"
#include "OStringStream.H"

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
        Info<< "PASS: " << what.c_str() << " (" << value << " vs " << tol
            << ")" << endl;
    }
    else
    {
        ++nFail;
        Info<< "FAIL: " << what.c_str() << " (" << value << " vs " << tol
            << ")" << endl;
    }
}


dictionary lawDict(const string& s)
{
    IStringStream is(s);

    return dictionary(is);
}


//- A small-strain law at one point, driven directly so that its state can
//  be read back
class directDriver
{
    autoPtr<mechanicalConstitutiveLaw> lawPtr_;
    mechanicalConstitutiveLawState state_;
    tensor gradD0_;

public:

    explicit directDriver(const dictionary& dict)
    :
        lawPtr_(mechanicalConstitutiveLaw::New(dict)),
        state_(1),
        gradD0_(tensor::zero)
    {
        // The declared defaults, at both times, as the manager applies them
        mechanicalConstitutiveLawStateSpec spec;
        lawPtr_->declareState(spec);

        forAll(spec.entries(), i)
        {
            const mechanicalConstitutiveLawStateSpec::entry& e =
                spec.entries()[i];

            if (e.typeName == pTraits<scalar>::typeName)
            {
                state_.scalarField(e.name) = e.scalarDefault;
                state_.scalarField0(e.name) = e.scalarDefault;
            }
            else if (e.typeName == pTraits<symmTensor>::typeName)
            {
                state_.symmTensorField(e.name) = e.symmTensorDefault;
                state_.symmTensorField0(e.name) = e.symmTensorDefault;
            }
            else
            {
                FatalErrorInFunction
                    << "Unexpected state type " << e.typeName
                    << exit(FatalError);
            }
        }
    }

    //- The stress at the end of a step to gradD, committing the last step
    symmTensor step(const tensor& gradD)
    {
        state_.storeOldTime();

        return evaluate(gradD);
    }

    //- The stress at the end of the current step, at gradD, without
    //  committing anything: may be called repeatedly within a step
    symmTensor evaluate(const tensor& gradD)
    {
        const List<tensor> g(1, gradD);
        const List<tensor> g0(1, gradD0_);
        const labelList pts(1, label(0));

        const UIndirectList<tensor> gView(g, pts);
        const UIndirectList<tensor> g0View(g0, pts);

        const smallStrainMechanicalConstitutiveLawKinematics kin(gView, g0View);

        symmTensorField sigma(1, symmTensor::zero);
        UIndirectList<symmTensor> sigmaView(sigma, pts);
        mechanicalConstitutiveLawResponse response(sigmaView);

        mechanicalConstitutiveLawInputs inputs(1.0, 0.0);

        lawPtr_->evaluate(kin, inputs, state_, response);

        return sigma[0];
    }

    //- Start a step: commit the last one
    void beginStep()
    {
        state_.storeOldTime();
    }

    scalar value(const word& name) const
    {
        return state_.getScalarField(name)[0];
    }

    scalar value0(const word& name) const
    {
        return state_.getScalarField0(name)[0];
    }
};


//- A law at the one integration point of a one-cell mesh, through the
//  manager, for the finite-strain composite
class managerDriver
{
    const dictionary mechanicalProperties_;
    mechanicalConstitutiveLawManager manager_;
    const integrationPointTopology& topo_;
    const scalar dt_;
    tensor F0_;
    tensor gradD0_;

    static dictionary makeDict(const word& name, const dictionary& d)
    {
        OStringStream os;
        os  << "planeStress no;" << nl << "mechanical" << nl << "(" << nl
            << name << nl << d << nl << ");" << nl;

        IStringStream is(os.str());

        return dictionary(is);
    }

public:

    managerDriver(const fvMesh& mesh, const word& name, const dictionary& d)
    :
        mechanicalProperties_(makeDict(name, d)),
        manager_(mesh, mechanicalProperties_),
        topo_
        (
            manager_.topologyFor(cellCentredIntegrationPointTopology::typeName)
        ),
        dt_(mesh.time().deltaTValue()),
        F0_(I),
        gradD0_(tensor::zero)
    {}

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


//- A diagonal strain
tensor diagStrain(const scalar a, const scalar b, const scalar c)
{
    return tensor(a, 0, 0, 0, b, 0, 0, 0, c);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const scalar E = 200e9;
    const scalar nu = 0.3;
    const scalar sigmaY0 = 400e6;
    const scalar H = 1e9;

    const string elastic
    (
        "rho rho [1 -3 0 0 0 0 0] 7800;"
        "E E [1 -1 -2 0 0 0 0] 200e9;"
        "nu nu [0 0 0 0 0 0 0] 0.3;"
    );

    // Linear hardening, sigmaY = 400 MPa + 1 GPa r, returned in closed form
    const string hardening
    (
        elastic
      + "file \"$FOAM_CASE/constant/linearHardening\"; outOfBounds clamp;"
    );

    // A nonlinear hardening table, which runs the Newton loop
    const string table
    (
        elastic
      + "file \"$FOAM_CASE/constant/nonlinearHardening\"; outOfBounds clamp;"
    );

    const string lemaitre("type linearElasticMisesPlasticLemaitreDamage; ");

    // Damage parameters that give a damage of order 0.1 over the paths below
    const string damage("S0 S0 [1 -1 -2 0 0 0 0] 200e6; b 1; ");

    // Uniaxial extension of the strain, past yield: every principal stress
    // is tensile and the triaxiality is above 1/3
    auto extension = [](const scalar t)
    {
        return diagStrain(0.2*t, 0, 0);
    };

    // 1. No damage: the plastic law itself
    {
        Info<< nl << "1. Without damage, linearElasticMisesPlastic" << endl;

        directDriver law
        (
            lawDict(lemaitre + table + damage + "epsilonD 1e6;")
        );
        directDriver plastic
        (
            lawDict("type linearElasticMisesPlastic; " + table)
        );

        scalar maxErr = 0;
        const label nSteps = 50;

        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            const scalar t = scalar(stepI)/nSteps;
            const tensor g
            (
                diagStrain(0.05*t, -0.01*t, 0.02*t)
              + tensor(0, 0.02*t, 0, 0.02*t, 0, 0, 0, 0, 0)
            );

            const symmTensor s1(law.step(g));
            const symmTensor s2(plastic.step(g));

            maxErr = max(maxErr, mag(s1 - s2)/mag(s2));
        }

        check
        (
            maxErr <= 1e-12 && law.value("damage") == 0,
            "with the threshold out of reach the stress is "
            "linearElasticMisesPlastic's and the damage zero",
            maxErr, 1e-12
        );
    }

    // 2. The discrete equations, for each form of S
    {
        Info<< nl << "2. The discrete yield and damage equations" << endl;

        wordList types(3);
        types[0] = "classic";
        types[1] = "Malcher";
        types[2] = "triaxiality";

        forAll(types, typeI)
        {
            const word& type = types[typeI];

            const scalar S0 = 2e6;
            string params("S0 S0 [1 -1 -2 0 0 0 0] 2e6; b 1; ");
            if (type == "Malcher")
            {
                params =
                    "S0Tension S0Tension [1 -1 -2 0 0 0 0] 2e6; "
                    "S0Shear S0Shear [1 -1 -2 0 0 0 0] 1e6; b 1; ";
            }

            directDriver law
            (
                lawDict
                (
                    lemaitre + hardening + params + "lemaitreType " + type
                  + ";"
                )
            );

            scalar maxYield = 0;
            scalar maxDamage = 0;
            scalar etaEnd = 0;
            const label nSteps = 100;

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                // Nearly isochoric extension with some shear, so that the
                // triaxiality stays between 0 and 1, where the three forms
                // of S differ, and the Lode parameter is not that of
                // uniaxial extension
                const scalar t = scalar(stepI)/nSteps;
                const tensor g
                (
                    diagStrain(0.2*t, -0.099*t, -0.099*t)
                  + tensor(0, 0.03*t, 0, 0.03*t, 0, 0, 0, 0, 0)
                );

                const symmTensor sigma(law.step(g));

                const scalar D = law.value("damage");
                const scalar D0 = law.value0("damage");
                const scalar r = law.value("epsilonPEq");
                const scalar dr = r - law.value0("epsilonPEq");

                // A failed point is set to 0.99, not solved for
                if (dr <= 0 || D >= 0.99)
                {
                    continue;
                }

                const scalar sigmaEq = Foam::sqrt(1.5)*mag(dev(sigma));
                const scalar sigmaYr = sigmaY0 + H*r;

                maxYield =
                    max(maxYield, mag(sigmaEq/(1 - D) - sigmaYr)/sigmaYr);

                // Equation 52, on the effective stress
                const scalar eta = tr(sigma)/(3*sigmaEq);
                etaEnd = eta;
                const scalar sv = sigmaEq/(1 - D);
                const scalar Y =
                   -sqr(sv)/(2*E)
                   *((2.0/3.0)*(1 + nu) + 3*(1 - 2*nu)*sqr(eta));

                scalar S = S0;
                if (type == "Malcher")
                {
                    const symmTensor s(dev(sigma));
                    const scalar xi = (27.0/2.0)*det(s)/pow3(sigmaEq);
                    S = S0/(3*mag(eta) + 2*(1 - sqr(xi)));
                }
                else if (type == "triaxiality")
                {
                    S = 2*S0/(1 + 3*eta);
                }

                const scalar dD = dr/(1 - D)*(-Y/S);

                maxDamage = max(maxDamage, mag(D - D0 - dD));
            }

            Info<< "    " << type << ": damage at the end "
                << law.value("damage") << ", triaxiality " << etaEnd << endl;

            check
            (
                law.value("damage") > 1e-2 && law.value("damage") < 0.99,
                type + ": the path damages the material without failing it, D",
                law.value("damage"), 1e-2
            );

            check
            (
                maxYield <= 1e-10,
                type + ": the returned stress is on the damaged yield "
                "surface",
                maxYield, 1e-10
            );

            check
            (
                maxDamage <= 1e-10,
                type + ": the damage satisfies the implicit update with "
                "Y from equation 52, |D - D0 - dD|",
                maxDamage, 1e-10
            );
        }
    }

    // 2b. Uniaxial stress, against the analytical solution
    {
        Info<< nl << "2b. Uniaxial stress, against Doghri's solution" << endl;

        // With b = 1 and linear hardening, under uniaxial stress, where the
        // triaxiality is 1/3 and Y = -sigmaEq^2/(2 E (1 - D)^2), the damage
        // integrates to
        //     (1 - D)^2 = 1 - sigmaY0^3/(3 E S0 H) [(1 + H r/sigmaY0)^3 - 1]
        const scalar S0 = 2e6;
        const scalar rMaxExpected = 0.3;

        scalarList errs(2);
        labelList nSteps(2);
        nSteps[0] = 100;
        nSteps[1] = 200;

        forAll(nSteps, k)
        {
            directDriver law
            (
                lawDict
                (
                    lemaitre + hardening + "S0 S0 [1 -1 -2 0 0 0 0] 2e6; b 1;"
                )
            );

            scalar lateral = 0;
            scalar maxLateralStress = 0;

            for (label stepI = 1; stepI <= nSteps[k]; ++stepI)
            {
                const scalar axial = 0.32*scalar(stepI)/nSteps[k];

                law.beginStep();

                // The lateral strain that leaves the lateral stresses zero,
                // by the secant method from the last step's
                scalar l0 = lateral;
                scalar l1 = lateral - 0.5*0.32/nSteps[k];
                scalar s0 =
                    law.evaluate(diagStrain(axial, l0, l0)).yy();
                scalar s1 =
                    law.evaluate(diagStrain(axial, l1, l1)).yy();

                for (label iter = 0; iter < 50 && mag(s1) > 1e-3; ++iter)
                {
                    const scalar l2 = l1 - s1*(l1 - l0)/(s1 - s0);
                    l0 = l1;
                    s0 = s1;
                    l1 = l2;
                    s1 = law.evaluate(diagStrain(axial, l1, l1)).yy();
                }

                lateral = l1;
                maxLateralStress = max(maxLateralStress, mag(s1));
            }

            const scalar r = law.value("epsilonPEq");
            const scalar D = law.value("damage");
            const scalar DExact =
                1
              - Foam::sqrt
                (
                    1
                  - pow3(sigmaY0)/(3*E*S0*H)
                   *(pow3(1 + H*r/sigmaY0) - 1)
                );

            errs[k] = mag(D - DExact)/DExact;

            Info<< "    " << nSteps[k] << " steps: r = " << r << ", D = "
                << D << ", exact " << DExact << ", max lateral stress "
                << maxLateralStress << endl;

            if (k == 0)
            {
                check
                (
                    r > 0.8*rMaxExpected && maxLateralStress < 1,
                    "uniaxial stress is held, to a plastic strain r",
                    r, rMaxExpected
                );
            }
        }

        check
        (
            errs[1] <= 5e-3,
            "the damage agrees with the analytical solution, relative "
            "error with 200 steps",
            errs[1], 5e-3
        );

        check
        (
            errs[0]/errs[1] > 1.6 && errs[0]/errs[1] < 2.5,
            "the error halves as the step is halved, ratio",
            errs[0]/errs[1], 2
        );
    }

    // 3. Crack closure
    {
        Info<< nl << "3. Crack closure" << endl;

        directDriver classic(lawDict(lemaitre + hardening + damage));
        directDriver closure
        (
            lawDict(lemaitre + hardening + damage + "crackClosure 0.2;")
        );

        for (label stepI = 1; stepI <= 100; ++stepI)
        {
            const tensor g(extension(scalar(stepI)/100));
            classic.step(g);
            closure.step(g);
        }

        const scalar diff =
            mag(classic.value("damage") - closure.value("damage"))
           /classic.value("damage");

        check
        (
            diff <= 1e-12,
            "with every principal stress tensile, crack closure changes "
            "nothing",
            diff, 1e-12
        );

        // Shear with lateral compression: one principal stress is
        // compressive and the triaxiality stays above the cut-off
        directDriver classicShear(lawDict(lemaitre + hardening + damage));
        directDriver closureShear
        (
            lawDict(lemaitre + hardening + damage + "crackClosure 0.2;")
        );

        for (label stepI = 1; stepI <= 100; ++stepI)
        {
            const scalar t = scalar(stepI)/100;
            const tensor g
            (
                tensor(0, 0.2*t, 0, 0.2*t, 0, 0, 0, 0, 0)
              + diagStrain(0.01*t, 0, 0)
            );
            classicShear.step(g);
            closureShear.step(g);
        }

        const scalar Dc = classicShear.value("damage");
        const scalar Dh = closureShear.value("damage");

        Info<< "    shear damage: classic " << Dc << ", crackClosure 0.2 "
            << Dh << endl;

        check
        (
            Dh < Dc && Dh > 0,
            "with a compressive principal stress, crack closure reduces the "
            "damage, D(h = 0.2)/D(h = 1)",
            Dh/Dc, 1
        );
    }

    // 4. The triaxiality cut-off
    {
        Info<< nl << "4. Triaxiality cut-off" << endl;

        directDriver law(lawDict(lemaitre + hardening + damage));
        directDriver plastic
        (
            lawDict("type linearElasticMisesPlastic; " + hardening)
        );

        scalar maxErr = 0;
        scalar eta = 0;

        for (label stepI = 1; stepI <= 100; ++stepI)
        {
            const tensor g(extension(-scalar(stepI)/100));
            const symmTensor s1(law.step(g));
            const symmTensor s2(plastic.step(g));

            maxErr = max(maxErr, mag(s1 - s2)/mag(s2));
            eta = tr(s1)/(3*Foam::sqrt(1.5)*mag(dev(s1)));
        }

        check
        (
            eta < -1.0/3.0,
            "uniaxial compression of the strain reaches eta < -1/3, eta",
            eta, -1.0/3.0
        );

        check
        (
            law.value("damage") == 0 && maxErr <= 1e-12,
            "below the cut-off there is no damage, and the stress is the "
            "undamaged law's",
            maxErr, 1e-12
        );
    }

    // 5. Failure
    {
        Info<< nl << "5. Failure" << endl;

        directDriver law
        (
            lawDict
            (
                lemaitre + hardening
              + "S0 S0 [1 -1 -2 0 0 0 0] 20e6; b 1; Dc 0.5;"
            )
        );

        label failedStep = -1;
        scalar sigmaBefore = 0;
        scalar sigmaAfter = 0;
        bool stays = true;

        for (label stepI = 1; stepI <= 100; ++stepI)
        {
            const symmTensor sigma(law.step(extension(scalar(stepI)/100)));

            if (failedStep < 0)
            {
                if (law.value("damage") >= 0.99)
                {
                    failedStep = stepI;
                    sigmaAfter = mag(sigma);
                }
                else
                {
                    sigmaBefore = mag(sigma);
                    if (law.value("damage") > 0.5)
                    {
                        stays = false;
                    }
                }
            }
            else if (law.value("damage") != 0.99)
            {
                stays = false;
            }
        }

        check
        (
            failedStep > 0 && stays,
            "a point whose damage passes Dc = 0.5 fails at 0.99 and stays "
            "failed, at step",
            failedStep, 100
        );

        check
        (
            sigmaAfter < 0.05*sigmaBefore,
            "the stress of a failed point collapses, |sigma| after/before",
            sigmaAfter/max(sigmaBefore, VSMALL), 0.05
        );
    }

    // 6. Convergence with the step
    {
        Info<< nl << "6. Step refinement" << endl;

        scalarList Ds(3);
        labelList nSteps(3);
        nSteps[0] = 25;
        nSteps[1] = 50;
        nSteps[2] = 100;

        forAll(nSteps, k)
        {
            directDriver law(lawDict(lemaitre + table + damage));

            for (label stepI = 1; stepI <= nSteps[k]; ++stepI)
            {
                const scalar t = scalar(stepI)/nSteps[k];
                law.step
                (
                    extension(t)
                  + tensor(0, 0.05*t, 0, 0.05*t, 0, 0, 0, 0, 0)
                );
            }

            Ds[k] = law.value("damage");
            Info<< "    " << nSteps[k] << " steps: D = " << Ds[k] << endl;
        }

        const scalar ratio = mag(Ds[0] - Ds[1])/max(mag(Ds[1] - Ds[2]), VSMALL);

        check
        (
            ratio > 1.6 && ratio < 2.5,
            "the damage converges at first order as the step is halved, "
            "ratio of successive differences",
            ratio, 2
        );
    }

    // 7. Under logStrain, on a coaxial path
    {
        Info<< nl << "7. Under logStrain" << endl;

        managerDriver logStrain
        (
            mesh,
            "logStrainLemaitre",
            lawDict
            (
                "type logStrain; rho rho [1 -3 0 0 0 0 0] 7800; "
                "mechanicalLaw { " + lemaitre + hardening + damage + " }"
            )
        );

        managerDriver smallStrain
        (
            mesh, "smallStrainLemaitre", lawDict(lemaitre + hardening + damage)
        );

        scalar maxErr = 0;

        for (label stepI = 1; stepI <= 60; ++stepI)
        {
            ++runTime;

            vector lambda(1 + 0.004*stepI, 1 - 0.001*stepI, 1);
            if (stepI > 30)
            {
                lambda = vector(1.12, 0.97 + 0.002*(stepI - 30), 1);
            }

            const tensor F(diagStrain(lambda.x(), lambda.y(), lambda.z()));
            const tensor lnF
            (
                diagStrain
                (
                    Foam::log(lambda.x()),
                    Foam::log(lambda.y()),
                    Foam::log(lambda.z())
                )
            );

            const symmTensor tau(det(F)*logStrain.finiteStrainStep(F));
            const symmTensor sigmaSmall(smallStrain.smallStrainStep(lnF));

            maxErr = max(maxErr, mag(tau - sigmaSmall)/mag(sigmaSmall));
        }

        check
        (
            maxErr <= 1e-9,
            "the Kirchhoff stress equals the small-strain law's at the "
            "logarithmic strain",
            maxErr, 1e-9
        );
    }

    // 8. Refused entries
    {
        Info<< nl << "8. Refused entries" << endl;

        FatalError.throwExceptions();
        FatalIOError.throwExceptions();

        auto refused = [](const string& s)
        {
            bool threw = false;

            try
            {
                autoPtr<mechanicalConstitutiveLaw> lawPtr
                (
                    mechanicalConstitutiveLaw::New(lawDict(s))
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            return threw;
        };

        check
        (
            refused
            (
                lemaitre + hardening + damage
              + "lemaitreType triaxiality; triaxialityCutOff no;"
            ),
            "lemaitreType triaxiality without the cut-off is refused",
            0, 0
        );

        check
        (
            refused(lemaitre + hardening + damage + "Dc 0.995;"),
            "a critical damage above 0.99 is refused", 0, 0
        );

        check
        (
            refused(lemaitre + hardening + damage + "lemaitreType Lemaitre;"),
            "an unknown lemaitreType is refused", 0, 0
        );

        check
        (
            refused(lemaitre + hardening + "lemaitreType Malcher;" + damage),
            "Malcher without S0Tension and S0Shear is refused", 0, 0
        );
    }

    Info<< nl << "Lemaitre checks: " << nPass << " passed, " << nFail
        << " failed" << endl;

    return (nFail == 0 ? 0 : 1);
}


// ************************************************************************* //
