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
    Test-mechanicalConstitutiveLawChecks

Description
    Checks each mechanicalConstitutiveLaw on its own, one integration point at
    a time, against properties every law of its kind must have.

    Test-mechanicalConstitutiveLaw exercises the manager on the laws a
    tutorial happens to use; this covers every law, whether or not a tutorial
    selects it. It runs on a one-cell case whose constant/lawChecks lists the
    laws, each as

        name
        {
            finiteStrain    yes;    // which update the law implements
            isotropic       yes;    // optional, default no
            linear          yes;    // optional, default no; small strain
            referenceChecks no;     // optional, default yes: checks 1-3
            yieldBounds     (243e6 600e6);  // optional: check 6
            fields          { T T [0 0 0 1 0 0 0] 0; }  // optional
            law             { type ...; ... }
        }

    and builds a fresh manager for each, from a mechanicalProperties with that
    one law. The fields are registered for the laws that read coupling inputs,
    at the values that make the reference state stress free.

    The following are checked, for each law:
      1. The reference state is stress free: F = I, or no strain.
      2. The tangent at the reference state, by central differences of the
         stress, has major symmetry. Every law here is elastic about its
         reference state, hyperelastic or linear, and so derives from an
         energy there.
      3. For an isotropic law, the law's scalar tangent is the normal
         stiffness C_xxxx of that tangent, lambda + 2 mu, and its deviatoric
         scalar tangent is (4/3) C_xyxy, (4/3) mu. The scalar tangents steer
         the segregated solvers; a wrong one costs iterations and does not
         show in any answer, which is why they are checked here. The second
         is the sharper test of the shear modulus: for a nearly
         incompressible law the bulk modulus dominates the first.
      4. For a finite-strain law, the Cauchy stress is objective: a rotation
         Q superposed on a 10% deformation gives Q sigma Q^T. This is at a
         strain large enough to take the plastic laws past yield.
      5. For a linear small-strain law, doubling the strain doubles the
         stress.
      6. For a plastic law given yieldBounds, the equivalent stress past yield
         - at 1% strain for a small-strain law, and at the 10% of check 4 for
         a finite-strain one, as a Kirchhoff stress - lies between the initial
         yield stress and the last of its hardening table. This is what runs
         the return mapping.

    Every law is checked for a stiffness that is positive and finite and a
    stress that is not zero under load, so that a law that is never reached
    cannot pass the checks above by returning nothing. And every law the
    runtime selection table holds must appear in constant/lawChecks, so a
    law added without an entry fails rather than going unchecked.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalConstitutiveLawManager.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "mat66.H"
#include "compatibilityFunctions.H"
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
    const word& law,
    const string& what,
    const scalar value,
    const scalar tol
)
{
    if (ok)
    {
        ++nPass;
        Info<< "PASS: " << law << ": " << what.c_str() << " ("
            << value << " <= " << tol << ")" << endl;
    }
    else
    {
        ++nFail;
        Info<< "FAIL: " << law << ": " << what.c_str() << " ("
            << value << " > " << tol << ")" << endl;
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


//- The displacement gradient of a unit engineering strain in Voigt component
//  i, in OpenFOAM's symmTensor order: XX, XY, XZ, YY, YZ, ZZ. A shear is
//  split between its two off-diagonal entries, so that it is a pure strain
tensor unitStrain(const label i)
{
    tensor g(tensor::zero);

    switch (i)
    {
        case symmTensor::XX: g.xx() = 1; break;
        case symmTensor::YY: g.yy() = 1; break;
        case symmTensor::ZZ: g.zz() = 1; break;
        case symmTensor::XY: g.xy() = 0.5; g.yx() = 0.5; break;
        case symmTensor::YZ: g.yz() = 0.5; g.zy() = 0.5; break;
        case symmTensor::XZ: g.xz() = 0.5; g.zx() = 0.5; break;
    }

    return g;
}


//- One law, evaluated at the first integration point of the cell-centred
//  topology of a one-cell mesh
class lawProbe
{
    mechanicalConstitutiveLawManager& manager_;
    const integrationPointTopology& topo_;
    const bool finiteStrain_;
    const scalar dt_;
    const label n_;

public:

    lawProbe
    (
        mechanicalConstitutiveLawManager& manager,
        const bool finiteStrain,
        const scalar dt
    )
    :
        manager_(manager),
        topo_
        (
            manager.topologyFor(cellCentredIntegrationPointTopology::typeName)
        ),
        finiteStrain_(finiteStrain),
        dt_(dt),
        n_(topo_.nIntegrationPoints())
    {}

    //- The stress for a displacement gradient (F = I + gradD at finite
    //  strain), from the reference state, with the scalar tangent if asked
    symmTensor stress
    (
        const tensor& gradD,
        scalar* scalarTangentPtr = nullptr,
        const tangentRequest scalarReq = tangentRequest::scalar
    )
    {
        symmTensorField sigma(n_, symmTensor::zero);
        scalarField K(n_, -GREAT);

        const tangentRequest req =
            scalarTangentPtr ? scalarReq : tangentRequest::none;

        if (finiteStrain_)
        {
            const tensor Fi(I + gradD);
            const tensorField F(n_, Fi);
            const tensorField Finv(n_, inv(Fi));
            const scalarField J(n_, det(Fi));
            const tensorField F0(n_, I);
            const tensorField Finv0(n_, I);
            const scalarField J0(n_, 1.0);

            manager_.updateStressFiniteStrain
            (
                topo_, F, F0, Finv, Finv0, J, J0, dt_, sigma,
                scalarTangentPtr ? &K : nullptr, nullptr, req
            );
        }
        else
        {
            const tensorField g(n_, gradD);
            const tensorField g0(n_, tensor::zero);

            manager_.updateStressSmallStrain
            (
                topo_, g, g0, dt_, sigma,
                scalarTangentPtr ? &K : nullptr, nullptr, req
            );
        }

        if (scalarTangentPtr)
        {
            *scalarTangentPtr = K[0];
        }

        return sigma[0];
    }

    //- The tangent at the reference state by central differences, in Voigt
    //  form against the engineering strain
    mat66 tangent(const scalar h)
    {
        mat66 C;

        for (label j = 0; j < 6; ++j)
        {
            const symmTensor sp(stress(h*unitStrain(j)));
            const symmTensor sm(stress(-h*unitStrain(j)));

            for (label i = 0; i < 6; ++i)
            {
                // A shear stress component pairs with the engineering shear
                // strain, which is what makes the Voigt matrix symmetric
                C(i, j) = (sp[i] - sm[i])/(2*h);
            }
        }

        return C;
    }
};

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const IOdictionary lawChecks
    (
        IOobject
        (
            "lawChecks",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar dt = runTime.deltaTValue();

    const wordList names(lawChecks.toc());

    HashSet<word> testedTypes;

    forAll(names, nameI)
    {
        const word& name = names[nameI];

        if (!lawChecks.isDict(name))
        {
            continue;
        }

        const dictionary& spec = lawChecks.subDict(name);

        const bool finiteStrain = Switch(spec.lookup("finiteStrain"));
        const bool isotropic =
            spec.lookupOrDefault<Switch>("isotropic", false);
        const bool linear = spec.lookupOrDefault<Switch>("linear", false);
        const bool referenceChecks =
            spec.lookupOrDefault<Switch>("referenceChecks", true);

        testedTypes.insert(word(spec.subDict("law").lookup("type")));

        Info<< nl << "Law checks: " << name << endl;

        // The coupling inputs, registered as another model would register
        // them, at the values of the law's reference state
        PtrList<volScalarField> inputs;

        if (spec.isDict("fields"))
        {
            const dictionary& fieldsDict = spec.subDict("fields");
            const wordList fieldNames(fieldsDict.toc());

            inputs.setSize(fieldNames.size());

            forAll(fieldNames, fieldI)
            {
                inputs.set
                (
                    fieldI,
                    new volScalarField
                    (
                        IOobject
                        (
                            fieldNames[fieldI],
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar(fieldsDict.lookup(fieldNames[fieldI]))
                    )
                );
            }
        }

        // A mechanicalProperties holding this one law
        OStringStream os;
        os  << "planeStress no;" << nl
            << "mechanical" << nl << "(" << nl
            << name << nl << spec.subDict("law") << nl
            << ");" << nl;

        IStringStream is(os.str());
        const dictionary mechanicalProperties(is);

        mechanicalConstitutiveLawManager manager(mesh, mechanicalProperties);

        lawProbe probe(manager, finiteStrain, dt);

        // Scales: the normal stiffness, from the tangent below, turns a
        // stress error into a strain-sized one
        const mat66 C(probe.tangent(1e-7));

        scalar Cmax = 0;
        for (label i = 0; i < 6; ++i)
        {
            for (label j = 0; j < 6; ++j)
            {
                Cmax = max(Cmax, mag(C(i, j)));
            }
        }

        // The law is reached and responds: a finite, positive stiffness, and
        // a stress under load that is not zero
        const tensor loadGradD
        (
            finiteStrain
          ? tensor(0.08, 0.03, -0.02, 0.01, -0.05, 0.04, 0.02, -0.01, 0.06)
          : tensor
            (
                1e-4, 0.3e-4, 0.0, 0.3e-4, -0.5e-4, 0.2e-4, 0.0, 0.2e-4, 0.7e-4
            )
        );

        {
            const scalar C11 = C(symmTensor::XX, symmTensor::XX);
            const scalar sigmaLoad = mag(probe.stress(loadGradD));

            const bool ok =
                Cmax > 0 && Cmax < GREAT && C11 > 0 && sigmaLoad > 0
             && sigmaLoad < GREAT;

            check
            (
                ok,
                name, "responds: positive, finite stiffness and nonzero stress",
                ok ? 0 : 1, 0
            );
        }

        // 1. The reference state is stress free
        if (referenceChecks)
        {
            const symmTensor sigma0(probe.stress(tensor::zero));

            check
            (
                mag(sigma0) <= 1e-12*Cmax,
                name, "reference state is stress free, |sigma|/C",
                mag(sigma0)/max(Cmax, VSMALL), 1e-12
            );
        }

        // 2. Major symmetry of the reference tangent
        if (referenceChecks)
        {
            scalar asym = 0;
            for (label i = 0; i < 6; ++i)
            {
                for (label j = 0; j < 6; ++j)
                {
                    asym = max(asym, mag(C(i, j) - C(j, i)));
                }
            }

            check
            (
                asym <= 1e-5*Cmax,
                name, "reference tangent has major symmetry, |C - C^T|/|C|",
                asym/max(Cmax, VSMALL), 1e-5
            );
        }

        // 3. The scalar tangent is lambda + 2 mu for an isotropic law
        if (isotropic && referenceChecks)
        {
            scalar K = -GREAT;
            probe.stress(tensor::zero, &K);

            const scalar C11 = C(symmTensor::XX, symmTensor::XX);
            const scalar err = mag(K - C11)/max(mag(C11), VSMALL);

            check
            (
                err <= 1e-4,
                name, "scalar tangent is C_xxxx, |K - C11|/C11", err, 1e-4
            );

            scalar Kdev = -GREAT;
            probe.stress(tensor::zero, &Kdev, tangentRequest::scalarDeviatoric);

            const scalar mu = C(symmTensor::XY, symmTensor::XY);
            const scalar errDev = mag(Kdev - 4*mu/3)/max(mag(4*mu/3), VSMALL);

            check
            (
                errDev <= 1e-4,
                name, "deviatoric scalar tangent is (4/3) C_xyxy", errDev, 1e-4
            );
        }

        // 4. Objectivity at finite strain
        if (finiteStrain)
        {
            const tensor& A = loadGradD;

            // 40 degrees, written without the fork-specific pi constant
            const tensor Q
            (
                rotation(vector(1, 2, 3), 40*Foam::atan(1.0)/45)
            );

            const symmTensor sigma(probe.stress(A));
            const symmTensor sigmaQ(probe.stress((Q & (I + A)) - I));
            const symmTensor expected(transform(Q, sigma));

            const scalar err =
                mag(sigmaQ - expected)/max(mag(sigma), VSMALL);

            check
            (
                err <= 1e-6,
                name, "objective under a superposed rotation", err, 1e-6
            );
        }

        // 5. Linearity of a linear small-strain law
        if (linear && !finiteStrain)
        {
            const tensor e
            (
                1e-4,  0.3e-4, 0.0,
                0.3e-4, -0.5e-4, 0.2e-4,
                0.0,    0.2e-4, 0.7e-4
            );

            const symmTensor s1(probe.stress(e));
            const symmTensor s2(probe.stress(2*e));

            const scalar err = mag(s2 - 2*s1)/max(mag(s2), VSMALL);

            check
            (
                err <= 1e-10,
                name, "doubling the strain doubles the stress", err, 1e-10
            );
        }

        // 6. Past yield, the equivalent stress is on the hardening curve
        if (spec.found("yieldBounds"))
        {
            const scalarList bounds(spec.lookup("yieldBounds"));

            // 1% at small strain, well past yield; at finite strain the 10%
            // of the objectivity check, measured as a Kirchhoff stress, which
            // is what the finite-strain plastic laws' yield function bounds
            const tensor g
            (
                finiteStrain
              ? loadGradD
              : tensor(1e-2, 0.3e-2, 0.0, 0.3e-2, -0.5e-2, 0.2e-2, 0.0, 0.2e-2,
                    -0.5e-2)
            );

            const scalar J = finiteStrain ? det(I + g) : 1.0;
            const symmTensor sigma(probe.stress(g));
            const scalar sigmaEq = J*Foam::sqrt(1.5)*mag(dev(sigma));

            const bool ok =
                sigmaEq >= (1 - 1e-6)*bounds[0]
             && sigmaEq <= (1 + 1e-6)*bounds[1];

            check
            (
                ok,
                name,
                "past yield, the equivalent stress lies within yieldBounds",
                sigmaEq, bounds[1]
            );
        }
    }

    // Every law that can be selected is checked
    {
        const wordList registered
        (
            mechanicalConstitutiveLaw::
            mechanicalConstitutiveLawConstructorTablePtr_->sortedToc()
        );

        forAll(registered, i)
        {
            if (!testedTypes.found(registered[i]))
            {
                ++nFail;
                Info<< "FAIL: " << registered[i] << ": no entry in "
                    << "constant/lawChecks, so it is not checked" << endl;
            }
        }

        if (testedTypes.empty())
        {
            ++nFail;
            Info<< "FAIL: constant/lawChecks lists no laws" << endl;
        }
    }

    Info<< nl << "Law checks: " << nPass << " passed, " << nFail
        << " failed" << endl;

    return (nFail == 0 ? 0 : 1);
}


// ************************************************************************* //
