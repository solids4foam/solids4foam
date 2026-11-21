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
    Test-phaseFieldFractureMechanicalConstitutiveLaw

Description
    Material-point checks of linearElasticMisesPlasticPhaseFieldFracture.

    The law is driven directly, at one point, with the phase field supplied
    as the input the solid model would solve for. The following are
    checked:
      1. With the phase field zero, the law is linearElasticMisesPlastic,
         step for step.
      2. The stress is degraded by g(d), plastic flow included, except for
         its volumetric part in compression; and a crack keeps the residual
         stiffness.
      3. On a path that yields and unloads, the history of the tensile
         elastic energy, the plastic work, and the source and coefficient of
         the phase field equation are those recomputed here; and a residual
         stiffness k scales the source by 1 - k.
      4. For the brittle model, with the phase field iterated to the
         solution of its equation for a uniform field, the uniaxial strain
         response is the analytical one, peaking at 9/16 of the intact
         stress.
      5. Under logStrain, with a phase field of 0.3, on a path whose
         principal directions stay fixed, the Kirchhoff stress equals the
         small-strain law's at the logarithmic strain.
      6. Inconsistent entries are refused.
      7. The law declares the phase field, with its coefficient.

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

    //- The phase field input
    scalarField input_;

public:

    explicit directDriver(const dictionary& dict)
    :
        lawPtr_(mechanicalConstitutiveLaw::New(dict)),
        state_(1),
        gradD0_(tensor::zero),
        input_(1, 0.0)
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
            inputs.setScalar(inputNames[i], input_);
        }

        lawPtr_->evaluate(kin, inputs, state_, response);

        return sigma[0];
    }

    //- Start a step: commit the last one
    void beginStep()
    {
        state_.storeOldTime();
    }

    //- Set the phase field input
    void setInput(const scalar d)
    {
        input_ = d;
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


} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const scalar E = 200e9;
    const scalar nu = 0.3;
    const scalar mu = E/(2*(1 + nu));
    const scalar K = E/(3*(1 - 2*nu));
    const scalar Gc = 2.7e3;
    const scalar l = 1e-3;

    const string elastic
    (
        "rho rho [1 -3 0 0 0 0 0] 7800;"
        "E E [1 -1 -2 0 0 0 0] 200e9;"
        "nu nu [0 0 0 0 0 0 0] 0.3;"
    );

    const string fracture
    (
        "type linearElasticMisesPlasticPhaseFieldFracture; "
        "Gc Gc [1 0 -2 0 0 0 0] 2.7e3; "
        "characteristicLength characteristicLength [0 1 0 0 0 0 0] 1e-3; "
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

    // A yield stress out of reach: the brittle model
    const string brittle
    (
        elastic
      + "file \"$FOAM_CASE/constant/elastic\"; outOfBounds clamp;"
    );

    // 1. Intact: the plastic law itself
    {
        Info<< nl << "1. Intact, linearElasticMisesPlastic" << endl;

        directDriver law(lawDict(fracture + table));
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
            maxErr <= 1e-12,
            "with the phase field zero the stress is "
            "linearElasticMisesPlastic's",
            maxErr, 1e-12
        );
    }

    // 2. The degraded stress
    {
        Info<< nl << "2. The degraded stress" << endl;

        const scalar d = 0.4;
        const scalar g = sqr(1 - d);

        // On a plastic path in tension: the whole stress is degraded
        {
            directDriver law(lawDict(fracture + table));
            directDriver plastic
            (
                lawDict("type linearElasticMisesPlastic; " + table)
            );
            law.setInput(d);

            scalar maxErr = 0;

            for (label stepI = 1; stepI <= 50; ++stepI)
            {
                const scalar t = scalar(stepI)/50;
                const tensor gradD(diagStrain(0.02*t, 0.005*t, 0.005*t));

                const symmTensor s1(law.step(gradD));
                const symmTensor s2(plastic.step(gradD));

                maxErr = max(maxErr, mag(s1 - g*s2)/mag(g*s2));
            }

            check
            (
                maxErr <= 1e-12,
                "under a tensile volumetric strain the stress is g(d) times "
                "the intact stress, plastic flow included",
                maxErr, 1e-12
            );
        }

        // Under compression, elastic: the volumetric part is not degraded
        {
            directDriver law(lawDict(fracture + brittle));
            law.setInput(d);

            const tensor gradD(diagStrain(-1e-3, -5e-4, 2e-4));
            const symmTensor eps(symm(gradD));
            const symmTensor sigma(law.step(gradD));
            const symmTensor expected
            (
                g*2*mu*dev(eps) + K*tr(eps)*symmTensor(I)
            );

            const scalar err = mag(sigma - expected)/mag(expected);

            check
            (
                err <= 1e-12,
                "under a compressive volumetric strain only the deviatoric "
                "stress is degraded",
                err, 1e-12
            );
        }

        // With a residual stiffness
        {
            directDriver law
            (
                lawDict(fracture + brittle + "residualStiffness 1e-3;")
            );
            law.setInput(1.0);

            const tensor gradD(diagStrain(1e-3, 0, 0));
            const symmTensor eps(symm(gradD));
            const symmTensor sigma(law.step(gradD));
            const symmTensor expected
            (
                1e-3*(2*mu*dev(eps) + K*tr(eps)*symmTensor(I))
            );

            const scalar err = mag(sigma - expected)/mag(expected);

            check
            (
                err <= 1e-12,
                "a crack keeps the residual stiffness",
                err, 1e-12
            );
        }
    }

    // 3. The crack driving force
    {
        Info<< nl << "3. The crack driving force" << endl;

        directDriver law
        (
            lawDict(fracture + hardening + "w0 w0 [1 -1 -2 0 0 0 0] 1e6;")
        );

        scalar maxErr = 0;
        scalar HMax = 0;
        scalar Wp = 0;
        const label nSteps = 60;

        // Load in tension, past yield, then unload
        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            const scalar t =
                stepI <= 40 ? scalar(stepI)/40 : 1 - scalar(stepI - 40)/40;
            const tensor gradD(diagStrain(0.01*t, 0.002*t, 0));

            law.step(gradD);

            // The tensile elastic energy, from the elastic strain
            const symmTensor epsE
            (
                symm(gradD) - law.symmValue("epsilonP")
            );
            const scalar psiPos =
                mu*(dev(epsE) && dev(epsE))
              + 0.5*K*sqr(max(tr(epsE), 0.0));

            HMax = max(HMax, psiPos);

            Wp +=
                law.value("sigmaY")
               *(law.value("epsilonPEq") - law.value0("epsilonPEq"));

            const scalar A = 4*l*(HMax + max(Wp - 1e6, 0.0))/Gc;

            maxErr =
                max
                (
                    maxErr,
                    max
                    (
                        mag(law.value("elasticEnergyHistory") - HMax)
                       /max(HMax, SMALL),
                        max
                        (
                            mag(law.value("plasticWork") - Wp)
                           /max(Wp, SMALL),
                            max
                            (
                                mag(law.value("phaseFieldSource") - A)/A,
                                mag
                                (
                                    law.value("phaseFieldCoefficient")
                                  - 1 - A
                                )/(1 + A)
                            )
                        )
                    )
                );
        }

        check
        (
            maxErr <= 1e-12 && Wp > 1e6
         && law.value("elasticEnergyHistory") > 0,
            "the history of the tensile elastic energy, the plastic work and "
            "the source and coefficient of the phase field equation hold "
            "through unloading",
            maxErr, 1e-12
        );
    }

    // 3b. With a residual stiffness, the driving force is scaled by 1 - k
    {
        directDriver law(lawDict(fracture + brittle));
        directDriver lawK
        (
            lawDict(fracture + brittle + "residualStiffness 0.1;")
        );

        const tensor gradD(diagStrain(1e-3, 2e-4, 0));
        law.step(gradD);
        lawK.step(gradD);

        const scalar err =
            mag
            (
                lawK.value("phaseFieldSource")
              - 0.9*law.value("phaseFieldSource")
            )/law.value("phaseFieldSource");

        check
        (
            err <= 1e-12,
            "with a residual stiffness k the driving force is scaled by "
            "1 - k, the slope of g",
            err, 1e-12
        );
    }

    // 4. The homogeneous brittle response. Under a uniform field the phase
    // field equation is d = A/(1 + A); with it iterated to that within each
    // step, uniaxial strain gives sigma = M eps/(1 + A)^2, with
    // A = 2 l M eps^2/Gc, whose peak is at A = 1/3
    {
        Info<< nl << "4. The homogeneous brittle response" << endl;

        directDriver law(lawDict(fracture + brittle));

        const scalar M = K + 4.0*mu/3.0;
        const scalar c = 2*l*M/Gc;
        const scalar epsPeak = Foam::sqrt(1.0/(3.0*c));

        scalar maxErr = 0;
        scalar peak = 0;
        scalar epsAtPeak = 0;
        const label nSteps = 400;

        for (label stepI = 1; stepI <= nSteps; ++stepI)
        {
            const scalar eps = 3*epsPeak*scalar(stepI)/nSteps;
            const tensor gradD(diagStrain(eps, 0, 0));

            law.beginStep();
            scalar dIn = 0;
            symmTensor sigma(symmTensor::zero);

            for (label iter = 0; iter < 500; ++iter)
            {
                law.setInput(dIn);
                sigma = law.evaluate(gradD);

                const scalar A = law.value("phaseFieldSource");
                const scalar dOut = A/(1 + A);

                if (mag(dOut - dIn) < 1e-15)
                {
                    break;
                }

                dIn = dOut;
            }

            const scalar A = c*sqr(eps);
            const scalar expected = M*eps/sqr(1 + A);

            maxErr = max(maxErr, mag(sigma.xx() - expected)/expected);

            if (sigma.xx() > peak)
            {
                peak = sigma.xx();
                epsAtPeak = eps;
            }
        }

        check
        (
            maxErr <= 1e-10,
            "with the phase field iterated to its uniform solution, the "
            "uniaxial strain response is the analytical one",
            maxErr, 1e-10
        );

        check
        (
            mag(epsAtPeak - epsPeak)/epsPeak < 0.01
         && mag(peak - (9.0/16.0)*M*epsPeak)/peak < 1e-4,
            "the stress peaks at A = 1/3, at 9/16 of the intact stress",
            epsAtPeak/epsPeak, 1
        );
    }

    // 5. Under logStrain, on a coaxial path
    {
        Info<< nl << "5. Under logStrain" << endl;

        // The phase field input, cracked, which the manager reads from the
        // registered field
        volScalarField phaseField
        (
            IOobject
            (
                "phaseField",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("d", dimless, 0.3)
        );

        managerDriver logStrain
        (
            mesh,
            "logStrainPhaseField",
            lawDict
            (
                "type logStrain; rho rho [1 -3 0 0 0 0 0] 7800; "
                "mechanicalLaw { " + fracture + hardening + " }"
            )
        );

        managerDriver smallStrain
        (
            mesh, "smallStrainPhaseField", lawDict(fracture + hardening)
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

    // 6. Refused entries
    {
        Info<< nl << "6. Refused entries" << endl;

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
                "type linearElasticMisesPlasticPhaseFieldFracture; "
                "characteristicLength characteristicLength [0 1 0 0 0 0 0] "
                "1e-3; " + hardening
            ),
            "a law without Gc is refused", 0, 0
        );

        check
        (
            refused
            (
                "type linearElasticMisesPlasticPhaseFieldFracture; "
                "Gc Gc [1 -1 -2 0 0 0 0] 2.7e3; "
                "characteristicLength characteristicLength [0 1 0 0 0 0 0] "
                "1e-3; " + hardening
            ),
            "Gc with the dimensions of a stress is refused", 0, 0
        );

        check
        (
            refused(fracture + hardening + "residualStiffness 1;"),
            "a residual stiffness of 1 is refused", 0, 0
        );

        check
        (
            refused(fracture + hardening + "w0 w0 [1 -1 -2 0 0 0 0] -1;"),
            "a negative plastic work threshold is refused", 0, 0
        );
    }

    // 7. The declared field
    {
        Info<< nl << "7. The declared field" << endl;

        directDriver law(lawDict(fracture + hardening));

        const List<mechanicalConstitutiveLaw::implicitGradientField>
            fields(law.law().implicitGradientFields());
        const wordList inputs(law.law().requiredScalarInputs());

        check
        (
            fields.size() == 1
         && fields[0].name == "phaseField"
         && fields[0].source == "phaseFieldSource"
         && fields[0].coefficient == "phaseFieldCoefficient"
         && mag(fields[0].length - 2*l) < SMALL
         && inputs.size() == 1
         && inputs[0] == "phaseField",
            "the law declares phaseField, of source A and coefficient 1 + A "
            "over the length 2 l, and reads it",
            fields.size(), 1
        );
    }

    Info<< nl << "Phase field checks: " << nPass << " passed, " << nFail
        << " failed" << endl;

    return (nFail == 0 ? 0 : 1);
}


// ************************************************************************* //
