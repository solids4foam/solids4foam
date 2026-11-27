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
    Test-GTNPlasticMechanicalConstitutiveLaw

Description
    Material-point checks of linearElasticGTNPlastic.

    The law is driven directly, at one point, over strain paths of many
    steps, so that its void and plastic history can be read back. The
    following are checked:
      1. With no voids, the law is linearElasticMisesPlastic, step for step.
      2. On a path of tensile mean stress, the state the law returns
         satisfies the discrete equations, recomputed here from the stress:
         the yield condition, the normality of the flow, the equivalence of
         plastic work and the void growth with nucleation; with linear and
         nonlinear hardening, with the shear term on a path with shear, and
         with coalescence.
      3. Under hydrostatic tension the voids grow, and under compression
         they close, the equations of 2 holding.
      4. A point whose void fraction would pass the failed void fraction
         fails, stays failed, and keeps about 1% of its strength; and a
         step whose void growth equation has a root near the start of the
         step returns it, even where the equation is negative again at the
         failed void fraction.
      5. The void fraction converges at first order as the step is refined.
      6. Under logStrain, on a path whose principal directions stay fixed,
         the Kirchhoff stress equals the small-strain law's at the
         logarithmic strain.
      7. Inconsistent entries are refused.
      8. The non-local form, with the non-local void fraction supplied as
         the input the solid model would solve for: with it zero the stress
         is linearElasticMisesPlastic's while the local voids nucleate; with
         the input iterated to the local void fraction within each step,
         which is the converged outer iteration for a vanishing length, the
         local form is recovered; and the law declares the field.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalConstitutiveLaw.H"
#include "mechanicalConstitutiveLawManager.H"
#include "linearElasticGTNPlasticMechanicalConstitutiveLaw.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "mathematicalConstants.H"
#include "Tuple2.H"
#include "IStringStream.H"
#include "OStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

label nPass = 0;
label nFail = 0;

#ifdef OPENFOAM_NOT_EXTEND
const scalar pi = constant::mathematical::pi;
#else
const scalar pi = mathematicalConstant::pi;
#endif

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

    //- The non-local void fraction input, for the non-local form
    scalarField fNonLocal_;

public:

    explicit directDriver(const dictionary& dict)
    :
        lawPtr_(mechanicalConstitutiveLaw::New(dict)),
        state_(1),
        gradD0_(tensor::zero),
        fNonLocal_(1, 0.0)
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

        const wordList inputNames(lawPtr_->requiredScalarInputs());
        forAll(inputNames, i)
        {
            inputs.setScalar(inputNames[i], fNonLocal_);
        }

        lawPtr_->evaluate(kin, inputs, state_, response);

        return sigma[0];
    }

    //- Start a step: commit the last one
    void beginStep()
    {
        state_.storeOldTime();
    }

    //- Set the non-local void fraction input
    void setNonLocal(const scalar f)
    {
        fNonLocal_ = f;
    }

    const mechanicalConstitutiveLaw& law() const
    {
        return lawPtr_();
    }

    scalar value(const word& name) const
    {
        return state_.getScalarField(name)[0];
    }

    scalar value0(const word& name) const
    {
        return state_.getScalarField0(name)[0];
    }

    symmTensor symmValue(const word& name) const
    {
        return state_.getSymmTensorField(name)[0];
    }

    symmTensor symmValue0(const word& name) const
    {
        return state_.getSymmTensorField0(name)[0];
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



//- The GTN equations, recomputed here from the returned state
class gtnEquations
{
    const scalar q1_;
    const scalar q2_;
    const scalar q3_;
    const scalar fN_;
    const scalar sN_;
    const scalar epsilonN_;
    const scalar kw_;
    const scalar fC_;
    const scalar fF_;

public:

    gtnEquations
    (
        const scalar q1, const scalar q2, const scalar fN, const scalar sN,
        const scalar epsilonN, const scalar kw, const scalar fC = 1,
        const scalar fF = 2
    )
    :
        q1_(q1), q2_(q2), q3_(sqr(q1)), fN_(fN), sN_(sN),
        epsilonN_(epsilonN), kw_(kw), fC_(fC), fF_(fF)
    {}

    scalar fStar(const scalar f) const
    {
        if (f <= fC_)
        {
            return f;
        }

        return fC_ + (1.0/q1_ - fC_)*(f - fC_)/(fF_ - fC_);
    }

    //- The relative residuals of the yield condition, the normality of the
    //  flow, the work equivalence and the void growth, the largest returned
    scalar residual
    (
        const symmTensor& sigma,
        const symmTensor& DEpsilonP,
        const scalar DEpsilonPEq,
        const scalar epsilonPEq,
        const scalar sigmaY,
        const scalar f,
        const scalar f0
    ) const
    {
        const scalar p = tr(sigma)/3.0;
        const symmTensor s(dev(sigma));
        const scalar q = Foam::sqrt(3.0/2.0)*mag(s);
        const scalar fs = fStar(f);
        const scalar a = 1.5*q2_*p/sigmaY;

        const scalar Phi =
            sqr(q/sigmaY) + 2.0*q1_*fs*Foam::cosh(a) - 1.0 - q3_*sqr(fs);

        const scalar dPhidp = 3.0*q1_*q2_*fs*Foam::sinh(a)/sigmaY;
        const scalar dPhidq = 2.0*q/sqr(sigmaY);

        const scalar DEpsilonV = tr(DEpsilonP);
        const scalar DEpsilonQ = Foam::sqrt(2.0/3.0)*mag(dev(DEpsilonP));

        // Normality: the deviatoric part along s, and the ratio of the parts
        const scalar dirErr =
            mag(dev(DEpsilonP) - DEpsilonQ*(3.0/2.0)*s/max(q, VSMALL))
           /max(mag(DEpsilonP), VSMALL);

        // Under a hydrostatic stress the deviatoric parts are round-off,
        // and normality is that the deviatoric flow vanishes
        scalar ratioErr = DEpsilonQ/max(mag(DEpsilonV), VSMALL);

        if (q > 1e-9*mag(p))
        {
            ratioErr =
                mag(DEpsilonV*dPhidq - DEpsilonQ*dPhidp)
               /max(mag(DEpsilonV*dPhidq) + mag(DEpsilonQ*dPhidp), VSMALL);
        }

        const scalar work = sigma && DEpsilonP;
        const scalar workErr =
            mag((1.0 - f)*sigmaY*DEpsilonPEq - work)/max(mag(work), VSMALL);

        scalar A = 0;
        if (fN_ > 0 && p > 0)
        {
            A =
                fN_/(sN_*Foam::sqrt(2.0*pi))
               *Foam::exp(-0.5*sqr((epsilonPEq - epsilonN_)/sN_));
        }

        scalar omega = 0;
        if (q > SMALL)
        {
            omega = 1.0 - sqr(27.0*det(s)/(2.0*pow3(q)));
        }

        const scalar growth =
            (1.0 - f)*DEpsilonV + A*DEpsilonPEq + kw_*f*omega*DEpsilonQ;

        const scalar growthErr =
            mag(f - f0 - growth)/max(mag(f - f0), VSMALL);

        return
            max
            (
                max(mag(Phi), dirErr),
                max(ratioErr, max(workErr, growthErr))
            );
    }
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const string elastic
    (
        "rho rho [1 -3 0 0 0 0 0] 7800;"
        "E E [1 -1 -2 0 0 0 0] 200e9;"
        "nu nu [0 0 0 0 0 0 0] 0.3;"
    );

    // Linear hardening, sigmaY = 400 MPa + 1 GPa r
    const string hardening
    (
        elastic
      + "file \"$FOAM_CASE/constant/linearHardening\"; outOfBounds clamp;"
    );

    // A nonlinear hardening table
    const string table
    (
        elastic
      + "file \"$FOAM_CASE/constant/nonlinearHardening\"; outOfBounds clamp;"
    );

    const string gtn("type linearElasticGTNPlastic; q1 1.5; q2 1; ");

    // Voids that grow and nucleate over the paths below
    const string voids
    (
        "f0 0.01; fN 0.04; sN 0.1; epsilonN 0.05; "
    );

    const gtnEquations equations(1.5, 1.0, 0.04, 0.1, 0.05, 0.0);

    // Extension under lateral strain: a tensile mean stress, so that the
    // voids grow
    auto triaxial = [](const scalar t)
    {
        return diagStrain(0.1*t, 0.03*t, 0.03*t);
    };

    // 1. No voids: the plastic law itself
    {
        Info<< nl << "1. Without voids, linearElasticMisesPlastic" << endl;

        directDriver law(lawDict(gtn + table + "f0 0;"));
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
            maxErr <= 1e-7 && law.value("voidFraction") == 0,
            "with no voids the stress is linearElasticMisesPlastic's and "
            "the void fraction stays zero",
            maxErr, 1e-7
        );
    }

    // 2. The discrete equations
    {
        Info<< nl << "2. The discrete GTN equations" << endl;

        List<Tuple2<string, string>> cases(2);
        cases[0] = Tuple2<string, string>("linear hardening", hardening);
        cases[1] = Tuple2<string, string>("nonlinear hardening", table);

        forAll(cases, caseI)
        {
            directDriver law(lawDict(gtn + cases[caseI].second() + voids));

            scalar maxErr = 0;
            const label nSteps = 100;

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                const scalar t = scalar(stepI)/nSteps;
                const symmTensor sigma(law.step(triaxial(t)));

                const scalar DEpsilonPEq =
                    law.value("epsilonPEq") - law.value0("epsilonPEq");

                if (DEpsilonPEq <= 0)
                {
                    continue;
                }

                maxErr =
                    max
                    (
                        maxErr,
                        equations.residual
                        (
                            sigma,
                            law.symmValue("epsilonP")
                          - law.symmValue0("epsilonP"),
                            DEpsilonPEq,
                            law.value("epsilonPEq"),
                            law.value("sigmaY"),
                            law.value("voidFraction"),
                            law.value0("voidFraction")
                        )
                    );
            }

            check
            (
                maxErr <= 1e-8 && law.value("voidFraction") > 0.02,
                "with " + cases[caseI].first() + ", the yield, flow, work "
                "and void growth equations hold, and the voids grow",
                maxErr, 1e-8
            );
        }

        // With the shear term, on a path with shear
        {
            directDriver law(lawDict(gtn + hardening + voids + "kw 2;"));
            const gtnEquations shearEquations(1.5, 1.0, 0.04, 0.1, 0.05, 2.0);
            directDriver noShear(lawDict(gtn + hardening + voids));

            scalar maxErr = 0;
            const label nSteps = 100;

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                const scalar t = scalar(stepI)/nSteps;
                const tensor g
                (
                    diagStrain(0.02*t, 0, 0)
                  + tensor(0, 0.1*t, 0, 0.1*t, 0, 0, 0, 0, 0)
                );

                const symmTensor sigma(law.step(g));
                noShear.step(g);

                const scalar DEpsilonPEq =
                    law.value("epsilonPEq") - law.value0("epsilonPEq");

                if (DEpsilonPEq <= 0)
                {
                    continue;
                }

                maxErr =
                    max
                    (
                        maxErr,
                        shearEquations.residual
                        (
                            sigma,
                            law.symmValue("epsilonP")
                          - law.symmValue0("epsilonP"),
                            DEpsilonPEq,
                            law.value("epsilonPEq"),
                            law.value("sigmaY"),
                            law.value("voidFraction"),
                            law.value0("voidFraction")
                        )
                    );
            }

            check
            (
                maxErr <= 1e-8
             && law.value("voidFraction") > noShear.value("voidFraction"),
                "with the shear term, the equations hold and the voids grow "
                "faster in shear",
                maxErr, 1e-8
            );
        }

        // With coalescence
        {
            directDriver law
            (
                lawDict
                (
                    gtn + hardening + voids
                  + "coalescence yes; fC 0.03; fF 0.2;"
                )
            );
            const gtnEquations coalescenceEquations
            (
                1.5, 1.0, 0.04, 0.1, 0.05, 0.0, 0.03, 0.2
            );

            scalar maxErr = 0;
            const label nSteps = 100;

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                const scalar t = scalar(stepI)/nSteps;
                const symmTensor sigma(law.step(triaxial(t)));

                const scalar DEpsilonPEq =
                    law.value("epsilonPEq") - law.value0("epsilonPEq");

                if (DEpsilonPEq <= 0)
                {
                    continue;
                }

                maxErr =
                    max
                    (
                        maxErr,
                        coalescenceEquations.residual
                        (
                            sigma,
                            law.symmValue("epsilonP")
                          - law.symmValue0("epsilonP"),
                            DEpsilonPEq,
                            law.value("epsilonPEq"),
                            law.value("sigmaY"),
                            law.value("voidFraction"),
                            law.value0("voidFraction")
                        )
                    );
            }

            check
            (
                maxErr <= 1e-8 && law.value("voidFraction") > 0.03,
                "with coalescence, past fC, the equations hold with the "
                "effective void fraction",
                maxErr, 1e-8
            );
        }
    }

    // 3. Void closure
    {
        Info<< nl << "3. Void closure" << endl;

        directDriver law(lawDict(gtn + hardening + "f0 0.05;"));

        const label nSteps = 50;
        scalar fMax = 0;

        // Dilatation, then compression past where it started
        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            const scalar t = scalar(stepI)/nSteps;
            law.step(diagStrain(0.02*t, 0.02*t, 0.02*t));
        }

        fMax = law.value("voidFraction");

        // With the equations checked while the voids close
        const gtnEquations closureEquations(1.5, 1.0, 0, 1, 0, 0);
        scalar maxErr = 0;

        for (label stepI = 1; stepI <= 2*nSteps; ++stepI)
        {
            const scalar t = scalar(stepI)/nSteps;
            const symmTensor sigma
            (
                law.step
                (
                    diagStrain(0.02 - 0.02*t, 0.02 - 0.02*t, 0.02 - 0.02*t)
                )
            );

            const scalar DEpsilonPEq =
                law.value("epsilonPEq") - law.value0("epsilonPEq");

            if (DEpsilonPEq <= 0)
            {
                continue;
            }

            maxErr =
                max
                (
                    maxErr,
                    closureEquations.residual
                    (
                        sigma,
                        law.symmValue("epsilonP") - law.symmValue0("epsilonP"),
                        DEpsilonPEq,
                        law.value("epsilonPEq"),
                        law.value("sigmaY"),
                        law.value("voidFraction"),
                        law.value0("voidFraction")
                    )
                );
        }

        const scalar fEnd = law.value("voidFraction");

        check
        (
            fMax > 0.05 && fEnd < fMax && fEnd >= 0 && maxErr <= 1e-8,
            "under hydrostatic tension the voids grow, and under "
            "compression they close, the equations holding",
            maxErr, 1e-8
        );
    }

    // 4. Failure
    {
        Info<< nl << "4. Failure" << endl;

        directDriver law
        (
            lawDict
            (
                gtn + hardening + voids + "coalescence yes; fC 0.03; fF 0.2;"
            )
        );

        const scalar fFailed =
            refCast<const linearElasticGTNPlasticMechanicalConstitutiveLaw>
            (
                law.law()
            ).fFailed();

        const label nSteps = 200;
        label failedStep = -1;
        scalar maxStressAfter = 0;

        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            const scalar t = scalar(stepI)/nSteps;
            const symmTensor sigma
            (
                law.step(diagStrain(0.2*t, 0.15*t, 0.15*t))
            );

            if (failedStep < 0 && law.value("voidFraction") >= fFailed)
            {
                failedStep = stepI;
            }
            else if (failedStep > 0)
            {
                maxStressAfter =
                    max(maxStressAfter, mag(sigma)/law.value("sigmaY"));
            }
        }

        check
        (
            failedStep > 0 && failedStep < nSteps - 10
         && law.value("voidFraction") == fFailed
         && maxStressAfter < 0.02,
            "a point fails at the failed void fraction, stays there, and "
            "keeps about 1% of its strength",
            maxStressAfter, 0.02
        );
    }

    // 4b. The root nearest the start of the step. Under hydrostatic
    // tension, with perfect plasticity, K = 100 sigmaY and coalescence from
    // 0.03, the void growth equation has a root near f0 and is negative
    // again at the failed void fraction: the step must not fail
    {
        Info<< nl << "4b. The nearest root" << endl;

        directDriver law
        (
            lawDict
            (
                "type linearElasticGTNPlastic; q1 1.5; q2 1; "
                "rho rho [1 -3 0 0 0 0 0] 7800; "
                "E E [1 -1 -2 0 0 0 0] 4.8e10; "
                "nu nu [0 0 0 0 0 0 0] 0.3; "
                "file \"$FOAM_CASE/constant/perfectPlasticity\"; "
                "outOfBounds clamp; "
                "f0 0.01; coalescence yes; fC 0.03; fF 0.031;"
            )
        );

        // tr(epsilon) = 2.800803385 sigmaY/K
        const scalar e = 0.02800803385/3.0;
        law.step(diagStrain(e, e, e));

        const scalar f = law.value("voidFraction");

        check
        (
            mag(f - 0.01002903349) < 1e-9,
            "a step whose void growth has a root near f0 and none just "
            "short of failure returns the near root rather than failing",
            f, 0.01002903349
        );
    }

    // 5. Convergence with the step
    {
        Info<< nl << "5. Convergence with the step" << endl;

        auto endVoidFraction = [&](const label nSteps)
        {
            directDriver law(lawDict(gtn + table + voids));

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                law.step(triaxial(scalar(stepI)/nSteps));
            }

            return law.value("voidFraction");
        };

        const scalar fRef = endVoidFraction(3200);
        const scalar e1 = mag(endVoidFraction(100) - fRef);
        const scalar e2 = mag(endVoidFraction(200) - fRef);
        const scalar e3 = mag(endVoidFraction(400) - fRef);

        Info<< "    errors " << e1 << " " << e2 << " " << e3 << endl;

        check
        (
            e1/e2 > 1.7 && e1/e2 < 2.4 && e2/e3 > 1.7 && e2/e3 < 2.4,
            "the void fraction converges at first order in the step",
            e2/e3, 2
        );
    }

    // 6. Under logStrain, on a coaxial path
    {
        Info<< nl << "6. Under logStrain" << endl;

        managerDriver logStrain
        (
            mesh,
            "logStrainGTN",
            lawDict
            (
                "type logStrain; rho rho [1 -3 0 0 0 0 0] 7800; "
                "mechanicalLaw { " + gtn + hardening + voids + " }"
            )
        );

        managerDriver smallStrain
        (
            mesh, "smallStrainGTN", lawDict(gtn + hardening + voids)
        );

        scalar maxErr = 0;

        for (label stepI = 1; stepI <= 60; ++stepI)
        {
            ++runTime;

            vector lambda(1 + 0.004*stepI, 1 + 0.001*stepI, 1);
            if (stepI > 30)
            {
                lambda = vector(1.12, 1.03 - 0.002*(stepI - 30), 1);
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

    // 7. Refused entries
    {
        Info<< nl << "7. Refused entries" << endl;

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
            refused(gtn + hardening + "f0 0.01; q3 3;"),
            "q3 above q1^2 is refused", 0, 0
        );

        check
        (
            refused(gtn + hardening + "f0 0.7;"),
            "an initial void fraction past failure is refused", 0, 0
        );

        check
        (
            refused(gtn + hardening + "f0 0.01; coalescence yes; fC 0.1;"),
            "coalescence without fF is refused", 0, 0
        );

        check
        (
            refused(gtn + hardening + "f0 0.01; fN 0.04;"),
            "nucleation without sN and epsilonN is refused", 0, 0
        );

        check
        (
            refused(gtn + hardening + voids + "nonLocal classic;"),
            "the non-local form without a characteristic length is refused",
            0, 0
        );

        check
        (
            refused
            (
                gtn + hardening + voids + "nonLocal effective; "
                "characteristicLength characteristicLength [0 1 0 0 0 0 0] "
                "1e-3;"
            ),
            "an unknown non-local form is refused", 0, 0
        );
    }

    // 8. The non-local form
    {
        Info<< nl << "8. The non-local form" << endl;

        const string classic
        (
            "nonLocal classic; "
            "characteristicLength characteristicLength [0 1 0 0 0 0 0] 1e-3; "
        );

        // 8a. With the input zero, the yield surface is von Mises', and the
        // local voids grow by nucleation alone
        {
            directDriver law(lawDict(gtn + table + voids + classic));
            directDriver plastic
            (
                lawDict("type linearElasticMisesPlastic; " + table)
            );

            scalar maxErr = 0;
            scalar maxGrowthErr = 0;
            const label nSteps = 50;

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                const scalar t = scalar(stepI)/nSteps;

                const symmTensor s1(law.step(triaxial(t)));
                const symmTensor s2(plastic.step(triaxial(t)));

                maxErr = max(maxErr, mag(s1 - s2)/mag(s2));

                const scalar epsilonPEq = law.value("epsilonPEq");
                const scalar DEpsilonPEq =
                    epsilonPEq - law.value0("epsilonPEq");
                const scalar A =
                    0.04/(0.1*Foam::sqrt(2.0*pi))
                   *Foam::exp(-0.5*sqr((epsilonPEq - 0.05)/0.1));

                maxGrowthErr =
                    max
                    (
                        maxGrowthErr,
                        mag
                        (
                            law.value("voidFraction")
                          - law.value0("voidFraction")
                          - A*DEpsilonPEq
                        )
                    );
            }

            check
            (
                maxErr <= 1e-7 && maxGrowthErr < 1e-12
             && law.value("voidFraction") > 0.01,
                "with the input zero the stress is linearElasticMisesPlastic's "
                "and the local voids nucleate",
                max(maxErr, maxGrowthErr), 1e-7
            );
        }

        // 8b. The input iterated to the local void fraction within each
        // step, the converged outer iteration of a vanishing length, gives
        // the local form
        {
            directDriver law(lawDict(gtn + table + voids + classic));
            directDriver local(lawDict(gtn + table + voids));

            scalar maxErr = 0;
            const label nSteps = 100;

            for (label stepI = 1; stepI <= nSteps; ++stepI)
            {
                const scalar t = scalar(stepI)/nSteps;
                const tensor g(triaxial(t));

                const symmTensor sLocal(local.step(g));

                law.beginStep();
                scalar fIn = law.value0("voidFraction");
                symmTensor s(symmTensor::zero);

                for (label iter = 0; iter < 200; ++iter)
                {
                    law.setNonLocal(fIn);
                    s = law.evaluate(g);

                    const scalar fOut = law.value("voidFraction");

                    if (mag(fOut - fIn) < 1e-15)
                    {
                        break;
                    }

                    fIn = fOut;
                }

                maxErr =
                    max
                    (
                        maxErr,
                        max
                        (
                            mag(s - sLocal)/mag(sLocal),
                            mag
                            (
                                law.value("voidFraction")
                              - local.value("voidFraction")
                            )
                        )
                    );
            }

            check
            (
                maxErr <= 1e-9,
                "iterated to the local void fraction, the non-local form "
                "gives the local form",
                maxErr, 1e-9
            );
        }

        // 8c. The declared field
        {
            directDriver law(lawDict(gtn + hardening + voids + classic));

            const List<mechanicalConstitutiveLaw::implicitGradientField>
                fields(law.law().implicitGradientFields());
            const wordList inputs(law.law().requiredScalarInputs());

            check
            (
                fields.size() == 1
             && fields[0].name == "voidFractionNonLocal"
             && fields[0].source == "voidFraction"
             && mag(fields[0].length - 1e-3) < SMALL
             && inputs.size() == 1
             && inputs[0] == "voidFractionNonLocal",
                "the law declares voidFractionNonLocal, of the local void "
                "fraction over the characteristic length, and reads it",
                fields.size(), 1
            );

            directDriver localLaw(lawDict(gtn + hardening + voids));

            check
            (
                localLaw.law().implicitGradientFields().empty()
             && localLaw.law().requiredScalarInputs().empty(),
                "the local form declares and reads no field",
                0, 0
            );
        }
    }

    Info<< nl << "GTN checks: " << nPass << " passed, " << nFail
        << " failed" << endl;

    return (nFail == 0 ? 0 : 1);
}


// ************************************************************************* //
