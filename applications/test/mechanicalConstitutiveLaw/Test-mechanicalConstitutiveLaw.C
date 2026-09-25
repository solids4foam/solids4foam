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
    Test-mechanicalConstitutiveLaw

Description
    Exercises the mechanicalConstitutiveLawManager evaluation paths on the
    case mesh, using the case constant/mechanicalProperties dictionary.

    The reference stress and tangent are recomputed here from E, nu and
    planeStress, rather than being taken from the constitutive laws, so that
    the closed-form checks are independent of the code under test.

    Run this on a case with more than one material (e.g. the layeredPipe
    tutorial, after setSet and setsToZones) to cover the multi-material
    addressing as well.

    The following are checked:
      1. linearElastic reproduces the closed-form isotropic stress and the
         scalar tangent, per material.
      2. The flat-list primitive agrees with the volTensorField overload.
      3. The CompactListList overload agrees with the volTensorField overload.
         Note that 2 and 3 check that the overloads agree, not that they are
         right: the overloads share the primitive, so a fault in the primitive
         moves them together and only check 1 sees it.
      4. updateTangentSmallStrain returns that same tangent and leaves the
         caller's stress storage untouched.
      5. dualFaceIntegrationPointTopology inverts the dual-face-to-cell map
         correctly, registerTopology is idempotent, and a stress evaluated on
         the dual faces uses the material of the owning primary cell - also
         on a second topology built from the boundary part of the map, the
         way vertexCentredLinGeomSolid evaluates its boundary dual faces.
      6. A fourth-order tangent on the dual faces matches the closed-form
         isotropic stiffness, including with more than one material.
      7. The finite-difference fourth-order tangent reproduces the analytical
         one, which also checks its Voigt shear convention.
      8. A tangent query leaves the constitutive state untouched, so a stress
         evaluated before and after intervening queries at wildly different
         kinematics is identical.
      9. updateScalarTangent agrees with the tangent from a stress update.
     9a. The finite-strain finite-difference tangent of a hyperelastic law
         reproduces the analytical small-strain isotropic tangent near F = I.
     10. The misuse guards fire: a fourth-order tangent on a topology that
         cannot carry one, a flat-list update on a topology whose integration
         points are shared between cells and where more than one material
         could claim them, a tangent request with no storage, a duplicate
         registerTopology key, and - where the case has more than one
         material - a face shared by two materials with no collapse rule to
         combine them, which must be refused where a rule is accepted.
     11. A field given per cell reaches a point on a face as the interpolated
         face value, for a face-centred and a compact face topology alike,
         including on a processor face, and a point in a cell as the cell
         value.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalConstitutiveLawManager.H"
#include "mechanicalConstitutiveLawInputs.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "mat66.H"
#include "Switch.H"
#include "IOmanip.H"
#include "compatibilityFunctions.H"
#include "mechanicalConstitutiveLaw.H"
#include "finiteStrainMechanicalConstitutiveLawKinematics.H"
#include "OFstream.H"
#include "mechanicalConstitutiveLawGather.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Number of failed checks
static label nFailed_ = 0;

//- Report the outcome of a check
void report(const word& name, const bool ok, const string& detail = "")
{
    if (ok)
    {
        Info<< "    PASS: " << name;
    }
    else
    {
        Info<< "    FAIL: " << name;
        nFailed_++;
    }

    if (!detail.empty())
    {
        Info<< " (" << detail.c_str() << ")";
    }

    Info<< endl;
}


//- Report a check on a maximum error against a tolerance
void reportError
(
    const word& name,
    const scalar error,
    const scalar tolerance
)
{
    OStringStream detail;
    detail << "max error " << error << ", tolerance " << tolerance;

    report(name, error <= tolerance, detail.str());
}


//- Largest difference between two lists, relative to the largest magnitude in
//  either of them.
//  Both are used to set the scale so that a result which wrongly collapses to
//  zero still reports a relative error of order one rather than of
//  order 1/SMALL
template<class Type>
scalar relativeDifference(const UList<Type>& a, const UList<Type>& b)
{
    if (a.size() != b.size())
    {
        FatalErrorInFunction
            << "Comparing lists of different size: " << a.size() << " and "
            << b.size() << exit(FatalError);
    }

    scalar scale = SMALL;
    scalar maxDiff = 0.0;

    forAll(a, i)
    {
        scale = max(scale, max(mag(a[i]), mag(b[i])));
        maxDiff = max(maxDiff, mag(a[i] - b[i]));
    }

    return maxDiff/scale;
}


//- A smooth, non-symmetric displacement gradient at the given position.
//  Polynomial rather than trigonometric, so that it is bounded and cheap and
//  needs no transcendental functions
tensor testGradD(const vector& p)
{
    const scalar s = 1e-3;

    const scalar x = p.x();
    const scalar y = p.y();
    const scalar z = p.z();

    return
        s*tensor
        (
            0.70 + 0.30*x*x,  0.20 - 0.15*y,    0.10 + 0.05*z,
            0.40 - 0.20*x,   -0.50 + 0.25*y*y,  0.30 - 0.10*z,
            0.15 + 0.10*x*y,  0.25 - 0.05*y*z,  0.60 + 0.20*z*z
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< nl << "Reading mechanicalProperties" << endl;

    const IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // ---------------------------------------------------------------------
    // Reference material properties per cell.
    // These are derived here from the dictionary rather than from the
    // constitutive laws, so that the closed-form checks below are independent
    // of the code under test
    // ---------------------------------------------------------------------

    const Switch planeStress
    (
        mechanicalProperties.lookupOrDefault<Switch>("planeStress", false)
    );

    const PtrList<entry> lawEntries(mechanicalProperties.lookup("mechanical"));

    scalarField refMu(mesh.nCells(), 0.0);
    scalarField refLambda(mesh.nCells(), 0.0);

    // The closed-form checks need a law whose stress and tangent are known
    // here. Everything else - path agreement, state preservation, the
    // topologies and the guards - applies to any law, and a history-dependent
    // law is the only thing that exercises the shadow state properly
    bool allLinearElastic = true;
    bool allNeoHookean = true;
    bool allStVenantKirchhoff = true;
    bool allMooneyRivlin = true;
    bool allNeoHookeanPlastic = true;
    bool allViscoelastic = true;
    bool allHGO = true;
    forAll(lawEntries, lawI)
    {
        const word type(lawEntries[lawI].dict().lookup("type"));

        if (type != "linearElastic")
        {
            allLinearElastic = false;
        }

        if (type != "neoHookeanElastic")
        {
            allNeoHookean = false;
        }

        if (type != "StVenantKirchhoffElastic")
        {
            allStVenantKirchhoff = false;
        }

        if (type != "MooneyRivlinElastic")
        {
            allMooneyRivlin = false;
        }

        if (type != "HolzapfelGasserOgdenElastic")
        {
            allHGO = false;
        }

        if (type != "neoHookeanElasticMisesPlastic")
        {
            allNeoHookeanPlastic = false;
        }

        if (type != "viscousHookeanElastic")
        {
            allViscoelastic = false;
        }
    }

    // Both are finite-strain-only laws: they implement no small-strain
    // evaluation, and both linearise to isotropic elasticity near F = I
    const bool allFiniteStrainOnly =
        allNeoHookean
     || allStVenantKirchhoff
     || allMooneyRivlin
     || allNeoHookeanPlastic;

    forAll(lawEntries, lawI)
    {
        const dictionary& lawDict = lawEntries[lawI].dict();

        if (!allLinearElastic && !allFiniteStrainOnly)
        {
            continue;
        }

        // A law may be given as E and nu or as mu and K, so accept either
        // here too rather than assuming the first form
        scalar E = 0.0;
        scalar nu = 0.0;

        if (allMooneyRivlin)
        {
            // Mooney-Rivlin is given as c10, c01, c11 and either K or nu.
            // Its small-strain limit has mu = 2*(c10 + c01), and where nu is
            // given the bulk modulus follows from E = 6*(c10 + c01), exactly
            // as the law itself derives them
            const scalar c10 =
                dimensionedScalar(lawDict.lookup("c10")).value();
            const scalar c01 =
                dimensionedScalar(lawDict.lookup("c01")).value();

            const scalar muIn = 2.0*(c10 + c01);

            scalar KIn = 0.0;
            if (lawDict.found("K"))
            {
                KIn = dimensionedScalar(lawDict.lookup("K")).value();
            }
            else
            {
                const scalar nuIn =
                    dimensionedScalar(lawDict.lookup("nu")).value();

                KIn = 6.0*(c10 + c01)/(3.0*(1.0 - 2.0*nuIn));
            }

            E = 9.0*KIn*muIn/(3.0*KIn + muIn);
            nu = (3.0*KIn - 2.0*muIn)/(2.0*(3.0*KIn + muIn));
        }
        else if (lawDict.found("mu") && lawDict.found("K"))
        {
            const scalar muIn = dimensionedScalar(lawDict.lookup("mu")).value();
            const scalar KIn = dimensionedScalar(lawDict.lookup("K")).value();

            E = 9.0*KIn*muIn/(3.0*KIn + muIn);
            nu = (3.0*KIn - 2.0*muIn)/(2.0*(3.0*KIn + muIn));
        }
        else
        {
            E = dimensionedScalar(lawDict.lookup("E")).value();
            nu = dimensionedScalar(lawDict.lookup("nu")).value();
        }

        const scalar mu = E/(2.0*(1.0 + nu));

        const scalar lambda =
            planeStress
          ? E*nu/((1.0 + nu)*(1.0 - nu))
          : E*nu/((1.0 + nu)*(1.0 - 2.0*nu));

        // The manager gives a single law the whole domain, and otherwise
        // matches each law to the cellZone of the same name
        if (lawEntries.size() == 1)
        {
            refMu = mu;
            refLambda = lambda;
        }
        else
        {
            const word& lawName = lawEntries[lawI].keyword();
            const label zoneID = mesh.cellZones().findZoneID(lawName);

            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "CellZone " << lawName << " not found. If this case has "
                    << "more than one material, run setSet and setsToZones "
                    << "before this test." << exit(FatalError);
            }

            const labelList& cells = mesh.cellZones()[zoneID];
            forAll(cells, i)
            {
                refMu[cells[i]] = mu;
                refLambda[cells[i]] = lambda;
            }
        }
    }

    Info<< "    materials: " << lawEntries.size()
        << ", planeStress: " << planeStress
        << ", cells: " << mesh.nCells()
        << ", closed-form checks: " << Switch(allLinearElastic) << endl;

    // ---------------------------------------------------------------------
    // Construct the manager and the test kinematics
    // ---------------------------------------------------------------------

    mechanicalConstitutiveLawManager manager(mesh, mechanicalProperties);

    volTensorField gradD
    (
        IOobject
        (
            "gradD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    );

    forAll(gradD, cellI)
    {
        gradD[cellI] = testGradD(mesh.C()[cellI]);
    }

    forAll(gradD.boundaryField(), patchI)
    {
        const vectorField& Cf = mesh.boundary()[patchI].Cf();
        tensorField& pGradD = Foam::boundaryFieldRef(gradD)[patchI];

        forAll(pGradD, faceI)
        {
            pGradD[faceI] = testGradD(Cf[faceI]);
        }
    }

    // The old-time gradient is not used by a linear elastic law, but it must
    // be supplied
    const volTensorField gradD0
    (
        IOobject
        (
            "gradD0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimless, tensor::zero)
    );

    const scalar dt = runTime.deltaTValue();

    // Some laws are finite strain only, and asking them for a small-strain
    // stress is a fatal error rather than a wrong answer. Ask once, quietly,
    // so that the checks which need small strain can be skipped for such a law
    // instead of taking the whole run down with them
    bool smallStrainCapable = true;
    {
        symmTensorField probeSigma(mesh.nCells(), symmTensor::zero);

        FatalError.throwExceptions();

        try
        {
            manager.updateStressSmallStrain
            (
                manager.topologyFor
                (
                    cellCentredIntegrationPointTopology::typeName
                ),
                Foam::primitiveField(gradD),
                Foam::primitiveField(gradD0),
                dt,
                probeSigma,
                nullptr,
                nullptr,
                tangentRequest::none
            );
        }
        catch (const Foam::error&)
        {
            smallStrainCapable = false;
        }

        FatalError.dontThrowExceptions();
    }

    // ------------------------------------------------------------------
    Info<< nl << "A declared isochoric split is an honest one" << endl;

    // A law that declares it can separate its isochoric stress from its
    // volumetric response is taken at its word by every mixed formulation, and
    // the declaration is a claim the framework cannot otherwise check: a law
    // written on the full deformation can return dev() of its total stress and
    // look exactly like one written on Cbar.
    //
    // What tells them apart is a superposed dilation. Under F -> c*F the
    // isochoric deformation Fbar = J^(-1/3)*F is unchanged, so a law whose
    // energy depends on Fbar alone returns the same Kirchhoff isochoric
    // stress, J*sigma_iso. A law whose energy sees the whole deformation does
    // not. The Cauchy stress itself does change, by 1/c^3, because it is per
    // current area - so the comparison is on J*sigma_iso and not on sigma_iso
    // Finite strain only: the check superposes a dilation on F, so it needs
    // laws that evaluate a finite-strain kinematics. A small-strain law may
    // separate its volumetric response perfectly well - linearElastic does -
    // and still have no finite-strain evaluation to call
    if
    (
        !smallStrainCapable
     && manager.allLawsProvideVolumetricSplit()
     && manager.allLawsHaveDilationInvariantIsochoricStress()
    )
    {
        const label n = mesh.nCells();

        volTensorField Fd
        (
            IOobject
            (
                "Fd",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("I", dimless, I)
        );
        volTensorField Fd0(Fd), Finvd(Fd), Finvd0(Fd);
        volScalarField Jd
        (
            IOobject
            (
                "Jd",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless, 1.0)
        );
        volScalarField Jd0(Jd);

        volSymmTensorField isoStress
        (
            IOobject
            (
                "isoStress",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
        );
        volScalarField volResponse
        (
            IOobject
            (
                "volResponse",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimPressure, 0.0)
        );

        // A deformation with shear and stretch, so the isochoric part is not
        // trivially zero, and the same one scaled by a pure dilation
        const tensor gradDbase
        (
            0.03, 0.012, 0.0,
            0.008, -0.02, 0.005,
            0.0, 0.004, 0.017
        );

        symmTensorField kirchhoffA(n, symmTensor::zero);
        symmTensorField kirchhoffB(n, symmTensor::zero);

        for (label pass = 0; pass < 2; ++pass)
        {
            const scalar c = (pass == 0 ? 1.0 : 1.19);

            forAll(Fd, cellI)
            {
                const tensor Fi = c*(I + gradDbase);
                Foam::primitiveFieldRef(Fd)[cellI] = Fi;
                Foam::primitiveFieldRef(Finvd)[cellI] = inv(Fi);
                Foam::primitiveFieldRef(Jd)[cellI] = det(Fi);
            }

            manager.updateStressFiniteStrainSplit
            (
                Fd, Fd0, Finvd, Finvd0, Jd, Jd0, dt, isoStress, volResponse
            );

            forAll(isoStress, cellI)
            {
                const symmTensor tau =
                    Foam::primitiveField(Jd)[cellI]
                   *Foam::primitiveField(isoStress)[cellI];

                if (pass == 0)
                {
                    kirchhoffA[cellI] = tau;
                }
                else
                {
                    kirchhoffB[cellI] = tau;
                }
            }
        }

        scalar maxDiff = 0.0;
        scalar scale = SMALL;

        forAll(kirchhoffA, cellI)
        {
            maxDiff = max(maxDiff, mag(kirchhoffA[cellI] - kirchhoffB[cellI]));
            scale = max(scale, mag(kirchhoffA[cellI]));
        }

        const scalar relDiff = maxDiff/scale;

        report
        (
            "the isochoric stress ignores a superposed dilation",
            relDiff < 1e-10,
            "relative change " + Foam::name(relDiff)
        );

        // And the isochoric stress must be trace-free for a law with no
        // spherical stress of its own. A law that adds one - an active tension
        // or a pore pressure - legitimately fails this, so it is only checked
        // where the total and the split differ by the volumetric response
        // alone
        scalar maxTrace = 0.0;

        forAll(isoStress, cellI)
        {
            maxTrace =
                max(maxTrace, mag(tr(Foam::primitiveField(isoStress)[cellI])));
        }

        // Asserted, not merely reported. This is the condition under which
        // a deviatoric projection of the total stress and the law's own
        // isochoric stress are the same thing - which is what the solid
        // models did before they could ask. A law that declares a dilation
        // invariant split and then returns a stress with a trace has quietly
        // made that substitution wrong wherever it is still used
        scalar maxStress = 0.0;

        forAll(isoStress, cellI)
        {
            maxStress =
                max
                (
                    maxStress,
                    mag(Foam::primitiveField(isoStress)[cellI])
                );
        }

        reportError
        (
            "the isochoric stress is trace free",
            maxTrace/max(maxStress, SMALL),
            1e-10
        );
    }
    else
    {
        Info<< "    SKIP: this check does not apply here - "
            << (
                   smallStrainCapable
                 ? "these are small-strain laws, and this check superposes a "
                   "dilation on the deformation gradient"
                 : manager.allLawsProvideVolumetricSplit()
                 ? "a law's split is not dilation invariant - it adds a "
                   "stress that is not derived from a potential, or its "
                   "yield surface scales with J"
                 : "no law here separates its isochoric and volumetric "
                   "responses"
               )
            << endl;
    }

    // ------------------------------------------------------------------
    Info<< nl << "A declared split recomposes the total stress" << endl;

    // The check above needs a law whose isochoric stress ignores a superposed
    // dilation, and skips the rest. This one needs nothing but the
    // declaration: whatever a law hands back as its isochoric stress and its
    // volumetric response, the two together must be the total it returns when
    // asked for the total, or a mixed formulation that replaces the second
    // with a solved pressure is solving for a different material. It is the
    // only check a law with history gets - neoHookeanElasticMisesPlastic,
    // whose yield surface scales with J - so the deformation is large enough
    // to take such a law well past yield
    if (!smallStrainCapable && manager.allLawsProvideVolumetricSplit())
    {
        volTensorField Fr
        (
            IOobject
            (
                "Fr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("I", dimless, I)
        );
        volTensorField Fr0(Fr), Finvr(Fr), Finvr0(Fr);
        volScalarField Jr
        (
            IOobject
            (
                "Jr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless, 1.0)
        );
        volScalarField Jr0(Jr);

        // Shear, stretch and a volume change, so that no part of the split is
        // trivially zero
        const tensor Fi
        (
            1.12, 0.07, 0.0,
            0.03, 0.93, 0.02,
            0.0, 0.01, 1.04
        );

        forAll(Fr, cellI)
        {
            Foam::primitiveFieldRef(Fr)[cellI] = Fi;
            Foam::primitiveFieldRef(Finvr)[cellI] = inv(Fi);
            Foam::primitiveFieldRef(Jr)[cellI] = det(Fi);
        }

        volSymmTensorField total
        (
            IOobject
            (
                "totalStress",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
        );
        volSymmTensorField iso(total);
        volScalarField vol
        (
            IOobject
            (
                "volResponseR",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimPressure, 0.0)
        );

        // A fully incompressible law has no total stress to compare with, and
        // the manager refuses to form one. Asked before evaluating, so that
        // any other failure of the evaluation fails this check rather than
        // being taken for that one
        const bool haveTotal = !manager.anyLawIncompressible();
        const string whyNoTotal
        (
            "a law is fully incompressible, so it has no total stress"
        );

        if (haveTotal)
        {
            manager.updateStressFiniteStrain
            (
                Fr, Fr0, Jr, Jr0, Finvr, Finvr0, dt, total
            );

            manager.updateStressFiniteStrainSplit
            (
                Fr, Fr0, Finvr, Finvr0, Jr, Jr0, dt, iso, vol
            );

            scalar maxErr = 0.0;
            scalar scale = SMALL;

            forAll(total, cellI)
            {
                const symmTensor& t = Foam::primitiveField(total)[cellI];
                const symmTensor recomposed
                (
                    Foam::primitiveField(iso)[cellI]
                  + Foam::primitiveField(vol)[cellI]*symmTensor(I)
                );

                maxErr = max(maxErr, mag(t - recomposed));
                scale = max(scale, mag(t));
            }

            reportError
            (
                "the isochoric stress plus the volumetric response is the "
                "total",
                maxErr/scale,
                1e-10
            );
        }
        else
        {
            Info<< "    SKIP: the total stress could not be formed here, so "
                << "there is nothing to compare with: " << whyNoTotal.c_str()
                << endl;
        }
    }
    else
    {
        Info<< "    SKIP: this check does not apply here - "
            << (
                   smallStrainCapable
                 ? "these are small-strain laws, and this check evaluates a "
                   "finite-strain split"
                 : "no law here separates its isochoric and volumetric "
                   "responses"
               )
            << endl;
    }

    // ------------------------------------------------------------------
    // The fibre term, against a closed form
    //
    // The two checks above are necessary and not sufficient: deleting the
    // fibre term entirely leaves a law that is still dilation invariant and
    // still trace free, so both would pass a law that had lost half its
    // physics. This pins the fibre contribution to a number.
    //
    // Under a uniaxial isochoric stretch F = diag(l, 1/sqrt(l), 1/sqrt(l))
    // with the fibres along x - fibreAngle zero, so both families coincide
    // with the stretch direction - the deformation is already isochoric, so
    // Fbar = F and J = 1. Then I4 = I6 = l^2, both families pull along x, and
    // eliminating the pressure by requiring zero lateral stress leaves
    //
    //     sigma_xx - sigma_yy = mu*(l^2 - 1/l)
    //                         + 4*k1*l^2*(l^2 - 1)*exp(k2*(l^2 - 1)^2)
    //
    // which is what the difference of the returned isochoric stresses must
    // be, since the volumetric response is spherical and cancels from it
    if (allHGO)
    {
        Info<< nl << "The fibre term against a closed form" << endl;

        const dictionary& hgoDict = lawEntries[0].dict();

        const scalar muVal =
            dimensionedScalar(hgoDict.lookup("mu")).value();
        const scalar k1Val =
            dimensionedScalar(hgoDict.lookup("k1")).value();
        // Dimensioned, as the law reads them and as case dictionaries write
        // them
        const scalar k2Val =
            dimensionedScalar(hgoDict.lookup("k2")).value();
        const scalar angle =
            dimensionedScalar(hgoDict.lookup("fibreAngle")).value();

        if (mag(angle) > SMALL)
        {
            Info<< "    SKIP: this check needs fibreAngle 0, and this case "
                << "sets " << angle << endl;
        }
        else
        {
            const label n = mesh.nCells();

            volTensorField Fd
            (
                IOobject
                (
                    "Fu",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedTensor("I", dimless, I)
            );
            volTensorField Fd0(Fd), Finvd(Fd), Finvd0(Fd);
            volScalarField Jd
            (
                IOobject
                (
                    "Ju",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("one", dimless, 1.0)
            );
            volScalarField Jd0(Jd);

            volSymmTensorField isoStress
            (
                IOobject
                (
                    "isoU",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
            );
            volScalarField volResponse
            (
                IOobject
                (
                    "volU",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimPressure, 0.0)
            );

            // Well past the exponential's knee, so that a wrong coefficient
            // or a missing push-forward shows up as a large error rather than
            // a small one
            const scalar lambda = 1.35;
            const scalar s = 1.0/Foam::sqrt(lambda);

            const tensor Fu(lambda, 0, 0, 0, s, 0, 0, 0, s);

            forAll(Fd, cellI)
            {
                Foam::primitiveFieldRef(Fd)[cellI] = Fu;
                Foam::primitiveFieldRef(Finvd)[cellI] = inv(Fu);
                Foam::primitiveFieldRef(Jd)[cellI] = det(Fu);
            }

            manager.updateStressFiniteStrainSplit
            (
                Fd, Fd0, Finvd, Finvd0, Jd, Jd0, dt, isoStress, volResponse
            );

            const scalar l2 = sqr(lambda);

            const scalar expected =
                muVal*(l2 - 1.0/lambda)
              + 4.0*k1Val*l2*(l2 - 1.0)*Foam::exp(k2Val*sqr(l2 - 1.0));

            scalar maxErr = 0.0;

            forAll(isoStress, cellI)
            {
                const symmTensor& sig =
                    Foam::primitiveField(isoStress)[cellI];

                const scalar got =
                    sig[symmTensor::XX] - sig[symmTensor::YY];

                maxErr = max(maxErr, mag(got - expected));
            }

            Info<< "        (uniaxial stretch " << lambda
                << ", expected sigma_xx - sigma_yy = " << expected << ')'
                << endl;

            reportError
            (
                "the fibre stress matches the closed form",
                maxErr/max(mag(expected), SMALL),
                1e-10
            );
        }
    }

    // ---------------------------------------------------------------------
    // Finite-strain finite-difference tangent
    //
    // Run first, and on its own, because a finite-strain law such as
    // neoHookeanElastic or StVenantKirchhoffElastic implements no small-strain
    // evaluation: every check below would fatal on it
    // ---------------------------------------------------------------------

    if (allFiniteStrainOnly)
    {
        Info<< nl << "Finite-strain finite-difference tangent" << endl;

        // A hyperelastic law linearises to isotropic elasticity as F
        // approaches the identity, so at a small deformation its
        // finite-difference spatial tangent must reproduce the analytical
        // small-strain tangent built from the same constants. That checks the
        // perturbation, the recomputed inverse and determinant, and the Voigt
        // convention at once
        // Face-centred: a cell-centred topology deliberately cannot carry a
        // fourth-order tangent, since the operators that consume one are
        // assembled from face fluxes
        const integrationPointTopology& topo =
            manager.topologyFor
            (
                faceCentredIntegrationPointTopology::typeName
            );

        const label n = topo.nIntegrationPoints();

        tensorField F(n, tensor::zero);
        tensorField Finv(n, tensor::zero);
        scalarField J(n, 0.0);
        const tensorField F0(n, I);
        const tensorField Finv0(n, I);
        const scalarField J0(n, 1.0);

        // A uniform small deformation is enough: the check is on the tangent,
        // not on any particular strain state.
        //
        // It must also stay below yield for an elasto-plastic law, since the
        // target below is the elastic tangent. That sets the scale: the
        // deviatoric trial stress is about 2*mu*strain, and cylinderExpansion
        // yields at 0.5 MPa with mu = 3.8 GPa, so a 1e-4 strain would already
        // be plastic and the elastic tangent would be the wrong target. 1e-6
        // is elastic for any realistic material.
        //
        // It costs no accuracy in the difference: the perturbation is
        // max(1e-8, 1e-6*mag(F - I)), which is at its 1e-8 floor for both
        // strains, so the stress difference being measured is the same size
        // either way
        const tensor gradDSmall
        (
            1e-6,  0.5e-6, 0.0,
            0.5e-6, -0.7e-6, 0.0,
            0.0,    0.0,   0.3e-6
        );

        forAll(F, ipI)
        {
            F[ipI] = I + gradDSmall;
            Finv[ipI] = inv(F[ipI]);
            J[ipI] = det(F[ipI]);
        }

        // Poisoned, so that an integration point the manager fails to
        // reach fails the comparison deterministically. mat66 is a POD, so
        // left alone its contents would be whatever memory held, which is a
        // test that passes or fails by luck
        List<mat66> fdC(n);
        forAll(fdC, ipI)
        {
            for (label i = 0; i < 6; ++i)
            {
                for (label j = 0; j < 6; ++j)
                {
                    fdC[ipI](i, j) = GREAT;
                }
            }
        }

        manager.updateTangentFiniteStrain
        (
            topo, F, F0, Finv, Finv0, J, J0, dt,
            nullptr, &fdC, tangentRequest::fourthOrderFiniteDifference
        );

        const label XX = symmTensor::XX;
        const label YY = symmTensor::YY;
        const label ZZ = symmTensor::ZZ;
        const label XY = symmTensor::XY;

        // Every integration point of the topology, boundary faces included.
        // The deformation is uniform, so the tangent is the same everywhere
        // and the boundary points are held to exactly the same standard as
        // the internal ones.
        //
        // This is deliberate, and is the regression guard for the defect of
        // section 8.14: the face-centred topology reports nFaces integration
        // points, and until the manager evaluated the boundary ones, every
        // entry from nInternalFaces upwards was left unwritten. Restricting
        // this loop to internal faces would hide a repeat of that defect
        const label nInternal = mesh.nInternalFaces();

        // Single material here, so the constants are uniform
        const scalar mu = refMu[0];
        const scalar lambda = refLambda[0];
        const scalar scale = lambda + 2.0*mu;

        scalar maxRelErrorInternal = 0.0;
        scalar maxRelErrorBoundary = 0.0;

        for (label ipI = 0; ipI < n; ++ipI)
        {
            const mat66& C = fdC[ipI];

            scalar e = 0.0;
            e = max(e, mag(C(XX, XX) - (lambda + 2.0*mu))/scale);
            e = max(e, mag(C(ZZ, ZZ) - (lambda + 2.0*mu))/scale);
            e = max(e, mag(C(XX, YY) - lambda)/scale);
            e = max(e, mag(C(XY, XY) - mu)/scale);
            e = max(e, mag(C(XX, XY))/scale);

            if (ipI < nInternal)
            {
                maxRelErrorInternal = max(maxRelErrorInternal, e);
            }
            else
            {
                maxRelErrorBoundary = max(maxRelErrorBoundary, e);
            }
        }

        reportError
        (
            "reproduces the small-strain isotropic tangent near F = I",
            maxRelErrorInternal,
            1e-3
        );

        reportError
        (
            "reproduces that tangent on boundary faces too",
            maxRelErrorBoundary,
            1e-3
        );

        // -----------------------------------------------------------------
        // Plasticity: the return mapping must actually return
        //
        // The check above stays deliberately below yield, so it exercises the
        // elastic predictor and nothing else. This one drives the material
        // well past yield and asserts the property that distinguishes a
        // working return map from a broken one: the deviatoric stress
        // saturates. Doubling the strain in the elastic range doubles the
        // deviatoric stress; once yielding, it must grow far more slowly,
        // and for a perfectly plastic curve hardly at all.
        //
        // This is deliberately independent of the hardening curve, so it does
        // not need to read the yield stress table the case supplies
        // -----------------------------------------------------------------
        if (allNeoHookeanPlastic)
        {
            Info<< nl << "Plastic return mapping" << endl;

            const scalar strainA = 1e-2;
            const scalar strainB = 2e-2;

            scalar magDevA = 0.0;
            scalar magDevB = 0.0;

            for (label pass = 0; pass < 2; ++pass)
            {
                const scalar e = (pass == 0 ? strainA : strainB);

                const tensor gradDLarge
                (
                    e,      0.5*e, 0.0,
                    0.5*e, -0.7*e, 0.0,
                    0.0,    0.0,   0.3*e
                );

                forAll(F, ipI)
                {
                    F[ipI] = I + gradDLarge;
                    Finv[ipI] = inv(F[ipI]);
                    J[ipI] = det(F[ipI]);
                }

                symmTensorField sigmaLarge(n, symmTensor::zero);

                // A tangent-free stress update. Both passes start from the
                // same old-time state, so this compares two trial states from
                // one history rather than a load path
                manager.updateStressFiniteStrain
                (
                    topo, F, F0, Finv, Finv0, J, J0, dt,
                    sigmaLarge, nullptr, nullptr,
                    tangentRequest::none
                );

                scalar acc = 0.0;
                for (label ipI = 0; ipI < n; ++ipI)
                {
                    acc = max(acc, mag(dev(sigmaLarge[ipI])));
                }

                if (pass == 0)
                {
                    magDevA = acc;
                }
                else
                {
                    magDevB = acc;
                }
            }

            // Elastic would give 2.0; a working return map gives close to 1
            const scalar growth = magDevB/max(magDevA, SMALL);

            reportError
            (
                "deviatoric stress saturates once yielding",
                mag(growth - 1.0),
                0.5
            );

            // And it must genuinely have yielded, or the check above is
            // vacuous: the stress must be far below the elastic prediction
            const scalar elasticPrediction = 2.0*refMu[0]*strainB;

            reportError
            (
                "the large deformation is well past yield",
                magDevB/elasticPrediction,
                0.5
            );
        }

        Info<< nl
            << "========================================================="
            << nl;

        if (nFailed_ == 0)
        {
            Info<< "All mechanicalConstitutiveLaw checks passed" << nl
                << "========================================================="
                << nl << nl << "End" << nl << endl;
            return 0;
        }

        Info<< nFailed_ << " mechanicalConstitutiveLaw check(s) FAILED" << nl
            << "========================================================="
            << nl << endl;

        return 1;
    }

    // ---------------------------------------------------------------------
    // 0b. Viscoelastic relaxation, and the time increment reaching the law
    //
    // This is the first law whose response depends on the time increment, so
    // it is the first end-to-end check that dt travels through the inputs
    // object. Two evaluations are made from the same rest state:
    //
    //   dt -> 0   no relaxation, so every Maxwell arm carries the full
    //             deviatoric stress and the response is the instantaneous
    //             elastic one
    //   dt -> inf every arm has relaxed to nothing and only the equilibrium
    //             branch remains
    //
    // The ratio of the two deviatoric stresses is therefore exactly
    // gammaInf = EInfinity/(EInfinity + sum(E)), which the case dictionary
    // gives, so this is an exact target rather than a bound.
    //
    // This must run before any other section, because both evaluations have
    // to start from the same rest state, and a later section commits a time
    // step after which the old-time state is no longer rest
    // ---------------------------------------------------------------------

    if (allViscoelastic)
    {
        Info<< nl << "0b. Viscoelastic relaxation" << endl;

        const dictionary& lawDict = lawEntries[0].dict();

        const scalar EInf =
            dimensionedScalar(lawDict.lookup("EInfinity")).value();
        const scalarList EArms(lawDict.lookup("E"));

        scalar E0 = EInf;
        forAll(EArms, i)
        {
            E0 += EArms[i];
        }

        const scalar gammaInf = EInf/E0;

        // A uniform deviatoric deformation
        forAll(gradD, cellI)
        {
            gradD[cellI] = tensor(1e-5, 0, 0, 0, -1e-5, 0, 0, 0, 0);
        }
        gradD.correctBoundaryConditions();

        volSymmTensorField sigmaInst
        (
            IOobject
            (
                "sigmaInst",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        );

        volSymmTensorField sigmaLong("sigmaLong", sigmaInst);

        const scalar tauMin =
            min(scalarList(lawDict.lookup("relaxationTimes")));

        manager.updateStressSmallStrain
        (
            gradD, gradD0, 1e-8*tauMin, sigmaInst
        );

        manager.updateStressSmallStrain
        (
            gradD, gradD0, 1e8*tauMin, sigmaLong
        );

        scalar maxErr = 0.0;
        scalar maxInst = 0.0;
        forAll(sigmaInst, cellI)
        {
            const scalar mInst = mag(dev(sigmaInst[cellI]));
            const scalar mLong = mag(dev(sigmaLong[cellI]));

            maxInst = max(maxInst, mInst);

            if (mInst > SMALL)
            {
                maxErr = max(maxErr, mag(mLong/mInst - gammaInf));
            }
        }

        reportError
        (
            "relaxes from the instantaneous to the long-term modulus",
            maxErr,
            1e-6
        );

        // And it must actually have relaxed, or the ratio check is vacuous
        report
        (
            "the instantaneous and long-term responses differ",
            maxInst > SMALL && gammaInf < 0.99
        );
    }



    if (!smallStrainCapable)
    {
        Info<< nl << "The remaining checks are small strain, and no law here "
            << "evaluates a small-strain" << nl << "kinematics, so they are "
            << "skipped." << nl;

        Info<< nl << "========================================================="
            << nl;

        if (nFailed_ == 0)
        {
            Info<< "All mechanicalConstitutiveLaw checks passed" << nl
                << "========================================================="
                << nl << endl;

            Info<< "End\n" << endl;

            return 0;
        }

        Info<< nFailed_ << " mechanicalConstitutiveLaw check(s) FAILED" << nl
            << "========================================================="
            << nl << endl;

        return 1;
    }

    // 1. Closed-form stress and scalar tangent through the volField overload
    // ---------------------------------------------------------------------

    Info<< nl << "1. Closed-form stress and scalar tangent" << endl;

    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
    );

    volScalarField impK
    (
        IOobject
        (
            "impK",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    );

    manager.updateStressSmallStrain
    (
        gradD, gradD0, dt, sigma, &impK, tangentRequest::scalar
    );

    {
        symmTensorField refSigma(mesh.nCells(), symmTensor::zero);
        scalarField refImpK(mesh.nCells(), 0.0);

        forAll(refSigma, cellI)
        {
            refSigma[cellI] =
                refMu[cellI]*twoSymm(gradD[cellI])
              + refLambda[cellI]*tr(gradD[cellI])*I;

            refImpK[cellI] = 2.0*refMu[cellI] + refLambda[cellI];
        }

        if (allLinearElastic)
        {
            reportError
            (
                "stress matches the closed form",
                relativeDifference(Foam::primitiveField(sigma), refSigma),
                1e-12
            );

            reportError
            (
                "scalar tangent matches 2*mu + lambda",
                relativeDifference(Foam::primitiveField(impK), refImpK),
                1e-12
            );
        }
        else
        {
            Info<< "    SKIP: closed-form stress and tangent "
                << "(not all materials are linearElastic)" << endl;
        }
    }

    // ---------------------------------------------------------------------
    // 2. The flat-list primitive agrees with the volField overload
    // ---------------------------------------------------------------------

    Info<< nl << "2. Flat-list primitive on a cell-centred topology" << endl;

    const integrationPointTopology& cellTopo =
        manager.topologyFor(cellCentredIntegrationPointTopology::typeName);

    {
        symmTensorField flatSigma(mesh.nCells(), symmTensor::zero);
        scalarField flatImpK(mesh.nCells(), 0.0);

        manager.updateStressSmallStrain
        (
            cellTopo,
            Foam::primitiveField(gradD),
            Foam::primitiveField(gradD0),
            dt,
            flatSigma,
            &flatImpK,
            nullptr,
            tangentRequest::scalar
        );

        reportError
        (
            "stress matches the volField overload",
            relativeDifference(Foam::primitiveField(sigma), flatSigma),
            1e-15
        );

        reportError
        (
            "scalar tangent matches the volField overload",
            relativeDifference(Foam::primitiveField(impK), flatImpK),
            1e-15
        );
    }

    // ---------------------------------------------------------------------
    // 3. The CompactListList overload agrees with the volField overload
    // ---------------------------------------------------------------------

    Info<< nl << "3. CompactListList overload with one point per cell" << endl;

    {
        const labelList sizes(mesh.nCells(), 1);

        CompactListList<tensor> compactGradD(sizes, tensor::zero);
        CompactListList<tensor> compactGradD0(sizes, tensor::zero);
        CompactListList<symmTensor> compactSigma(sizes, symmTensor::zero);
        scalarField compactImpK(mesh.nCells(), 0.0);

        forAll(gradD, cellI)
        {
            compactGradD(cellI, 0) = gradD[cellI];
            compactGradD0(cellI, 0) = gradD0[cellI];
        }

        manager.updateStressSmallStrain
        (
            compactGradD,
            compactGradD0,
            dt,
            compactSigma,
            &compactImpK,
            tangentRequest::scalar
        );

        // A layout with the same number of values and different rows. The
        // flat index of a point then differs between the two lists, so one
        // cell's strain would be read as another's, and every list is still
        // exactly the length it should be - nothing downstream can notice
        if (mesh.nCells() >= 2)
        {
            labelList crossed(sizes);
            crossed[0] = sizes[0] + 1;
            crossed[1] = sizes[1] - 1;

            if (crossed[1] >= 0)
            {
                const CompactListList<tensor> oddGradD0
                (
                    crossed, tensor::zero
                );

                bool threw = false;
                FatalError.throwExceptions();
                try
                {
                    manager.updateStressSmallStrain
                    (
                        compactGradD,
                        oddGradD0,
                        dt,
                        compactSigma,
                        &compactImpK,
                        tangentRequest::scalar
                    );
                }
                catch (const Foam::error&)
                {
                    threw = true;
                }
                FatalError.dontThrowExceptions();

                report
                (
                    "a compact layout with the same length and different "
                    "rows is rejected",
                    threw
                );
            }
        }

        reportError
        (
            "stress matches the volField overload",
            relativeDifference(Foam::primitiveField(sigma), compactSigma.m()),
            1e-15
        );

        reportError
        (
            "scalar tangent matches the volField overload",
            relativeDifference(Foam::primitiveField(impK), compactImpK),
            1e-15
        );
    }

    // ---------------------------------------------------------------------
    // 4. updateTangentSmallStrain leaves the caller's stress alone
    // ---------------------------------------------------------------------

    Info<< nl << "4. Tangent-only update" << endl;

    {
        const symmTensor sentinel(1, 2, 3, 4, 5, 6);

        symmTensorField untouched(mesh.nCells(), sentinel);
        scalarField tangentOnly(mesh.nCells(), 0.0);

        manager.updateTangentSmallStrain
        (
            cellTopo,
            Foam::primitiveField(gradD),
            Foam::primitiveField(gradD0),
            dt,
            &tangentOnly,
            nullptr,
            tangentRequest::scalar
        );

        bool stressUntouched = true;
        forAll(untouched, cellI)
        {
            if (untouched[cellI] != sentinel)
            {
                stressUntouched = false;
                break;
            }
        }

        report("the caller's stress storage is untouched", stressUntouched);

        reportError
        (
            "tangent matches the volField overload",
            relativeDifference(Foam::primitiveField(impK), tangentOnly),
            1e-15
        );
    }

    // ---------------------------------------------------------------------
    // 5. dualFaceIntegrationPointTopology
    // ---------------------------------------------------------------------

    Info<< nl << "5. Dual-face topology" << endl;

    // A synthetic dual-face-to-cell map. The dual faces are interleaved across
    // the cells so that the inversion is not trivially ordered, and a handful
    // of trailing entries stand in for boundary dual faces, which are not
    // integration points of the internal topology. They are spread over the
    // whole cell range, so that they do not all share one material
    const label nPerCell = 3;
    const label nInternalDualFaces = nPerCell*mesh.nCells();
    const label nBoundaryDualFaces = 5;

    labelList dualFaceToCell(nInternalDualFaces + nBoundaryDualFaces, 0);
    for (label i = 0; i < nInternalDualFaces; ++i)
    {
        dualFaceToCell[i] = i % mesh.nCells();
    }
    for (label i = 0; i < nBoundaryDualFaces; ++i)
    {
        dualFaceToCell[nInternalDualFaces + i] =
            mesh.nCells() - 1 - (i*mesh.nCells())/nBoundaryDualFaces;
    }

    const integrationPointTopology& dualTopo =
        manager.registerTopology
        (
            "testDualFaces",
            autoPtr<integrationPointTopology>
            (
                new dualFaceIntegrationPointTopology
                (
                    mesh, dualFaceToCell, nInternalDualFaces
                )
            )
        );

    report
    (
        "only internal dual faces are integration points",
        dualTopo.nIntegrationPoints() == nInternalDualFaces
    );

    report
    (
        "a fourth-order tangent is supported",
        dualTopo.supportsFourthOrderTangent()
    );

    report
    (
        "integration points are not shared between cells",
        !dualTopo.requiresUniqueIntegrationPointsPerMaterial()
    );

    {
        // Every internal dual face must appear exactly once, under the cell
        // the map assigns it to
        labelList timesSeen(nInternalDualFaces, 0);
        bool addressingCorrect = true;

        for (label cellI = 0; cellI < mesh.nCells(); ++cellI)
        {
            const labelUList ips = dualTopo.cellIntegrationPointIDs(cellI);

            forAll(ips, i)
            {
                const label ip = ips[i];

                if (ip < 0 || ip >= nInternalDualFaces)
                {
                    addressingCorrect = false;
                    break;
                }

                if (dualFaceToCell[ip] != cellI)
                {
                    addressingCorrect = false;
                    break;
                }

                timesSeen[ip]++;
            }
        }

        forAll(timesSeen, ip)
        {
            if (timesSeen[ip] != 1)
            {
                addressingCorrect = false;
                break;
            }
        }

        report
        (
            "the dual-face-to-cell map is inverted correctly",
            addressingCorrect
        );
    }

    report
    (
        "registerTopology is idempotent",
        &manager.registerTopology
        (
            "testDualFaces",
            autoPtr<integrationPointTopology>
            (
                new dualFaceIntegrationPointTopology
                (
                    mesh, dualFaceToCell, nInternalDualFaces
                )
            )
        ) == &dualTopo
    );

    // vertexCentredLinGeomSolid evaluates its boundary dual faces as a second
    // dual-face topology, built from the boundary part of the same map, so
    // that they keep their own constitutive state
    const labelList boundaryDualFaceToCell
    (
        SubList<label>(dualFaceToCell, nBoundaryDualFaces, nInternalDualFaces)
    );

    const integrationPointTopology& dualBoundaryTopo =
        manager.registerTopology
        (
            "testDualBoundaryFaces",
            autoPtr<integrationPointTopology>
            (
                new dualFaceIntegrationPointTopology
                (
                    mesh, boundaryDualFaceToCell, nBoundaryDualFaces
                )
            )
        );

    report
    (
        "the boundary dual faces form a separate topology",
        &dualBoundaryTopo != &dualTopo
     && dualBoundaryTopo.nIntegrationPoints() == nBoundaryDualFaces
    );

    // Kinematics on the dual faces, then a stress check per dual face against
    // the material of the primary cell that owns it
    tensorField dualGradD(nInternalDualFaces, tensor::zero);
    const tensorField dualGradD0(nInternalDualFaces, tensor::zero);

    forAll(dualGradD, dualFaceI)
    {
        // Offset the position so that neighbouring dual faces of one cell do
        // not share a value
        dualGradD[dualFaceI] =
            testGradD
            (
                mesh.C()[dualFaceToCell[dualFaceI]]
              + vector(0.001*dualFaceI, 0, 0)
            );
    }

    {
        symmTensorField dualSigma(nInternalDualFaces, symmTensor::zero);

        manager.updateStressSmallStrain
        (
            dualTopo,
            dualGradD,
            dualGradD0,
            dt,
            dualSigma,
            nullptr,
            nullptr,
            tangentRequest::none
        );

        symmTensorField refDualSigma(nInternalDualFaces, symmTensor::zero);
        forAll(refDualSigma, dualFaceI)
        {
            const label cellI = dualFaceToCell[dualFaceI];

            refDualSigma[dualFaceI] =
                refMu[cellI]*twoSymm(dualGradD[dualFaceI])
              + refLambda[cellI]*tr(dualGradD[dualFaceI])*I;
        }

        if (allLinearElastic)
        {
            reportError
            (
                "dual-face stress uses the owning cell's material",
                relativeDifference(dualSigma, refDualSigma),
                1e-12
            );
        }
        else
        {
            Info<< "    SKIP: dual-face closed-form stress "
                << "(not all materials are linearElastic)" << endl;
        }

        // The same on the boundary topology, whose integration point i is
        // boundary dual face i
        tensorField bGradD(nBoundaryDualFaces, tensor::zero);
        const tensorField bGradD0(nBoundaryDualFaces, tensor::zero);
        forAll(bGradD, i)
        {
            bGradD[i] =
                testGradD
                (
                    mesh.C()[boundaryDualFaceToCell[i]]
                  + vector(0, 0.001*i, 0)
                );
        }

        symmTensorField bSigma(nBoundaryDualFaces, symmTensor::zero);

        manager.updateStressSmallStrain
        (
            dualBoundaryTopo, bGradD, bGradD0, dt, bSigma
        );

        symmTensorField refBSigma(nBoundaryDualFaces, symmTensor::zero);
        forAll(refBSigma, i)
        {
            const label cellI = boundaryDualFaceToCell[i];

            refBSigma[i] =
                refMu[cellI]*twoSymm(bGradD[i])
              + refLambda[cellI]*tr(bGradD[i])*I;
        }

        if (allLinearElastic)
        {
            reportError
            (
                "boundary dual-face stress uses the owning cell's material",
                relativeDifference(bSigma, refBSigma),
                1e-12
            );
        }
        else
        {
            Info<< "    SKIP: boundary dual-face closed-form stress "
                << "(not all materials are linearElastic)" << endl;
        }
    }

    // ---------------------------------------------------------------------
    // 6. A fourth-order tangent on the dual faces
    // ---------------------------------------------------------------------

    Info<< nl << "6. Fourth-order tangent on the dual faces" << endl;

    if (!allLinearElastic)
    {
        Info<< "    SKIP: no analytical fourth-order tangent for these "
            << "materials" << endl;
    }
    else
    {
        List<mat66> dualC(nInternalDualFaces);

        manager.updateTangentSmallStrain
        (
            dualTopo,
            dualGradD,
            dualGradD0,
            dt,
            nullptr,
            &dualC,
            tangentRequest::fourthOrder
        );

        const label XX = symmTensor::XX;
        const label YY = symmTensor::YY;
        const label ZZ = symmTensor::ZZ;
        const label XY = symmTensor::XY;

        scalar maxRelError = 0.0;

        forAll(dualC, dualFaceI)
        {
            const label cellI = dualFaceToCell[dualFaceI];

            const scalar mu = refMu[cellI];
            const scalar lambda = refLambda[cellI];
            const scalar scale = lambda + 2.0*mu;

            const mat66& C = dualC[dualFaceI];

            maxRelError =
                max(maxRelError, mag(C(XX, XX) - (lambda + 2.0*mu))/scale);
            maxRelError =
                max(maxRelError, mag(C(YY, YY) - (lambda + 2.0*mu))/scale);
            maxRelError =
                max(maxRelError, mag(C(ZZ, ZZ) - (lambda + 2.0*mu))/scale);
            maxRelError = max(maxRelError, mag(C(XX, YY) - lambda)/scale);
            maxRelError = max(maxRelError, mag(C(YY, ZZ) - lambda)/scale);
            maxRelError = max(maxRelError, mag(C(XY, XY) - mu)/scale);

            // The tangent is isotropic, so there is no normal-shear coupling
            maxRelError = max(maxRelError, mag(C(XX, XY))/scale);
        }

        if (allLinearElastic)
        {
            reportError
            (
                "the tangent matches the closed-form isotropic stiffness",
                maxRelError,
                1e-12
            );
        }
        else
        {
            Info<< "    SKIP: closed-form isotropic stiffness "
                << "(not all materials are linearElastic)" << endl;
        }
    }

    // ---------------------------------------------------------------------
    // 7. Finite-difference fourth-order tangent
    // ---------------------------------------------------------------------

    Info<< nl << "7. Finite-difference fourth-order tangent" << endl;

    {
        List<mat66> fdC(nInternalDualFaces);
        List<mat66> fdCagain(nInternalDualFaces);

        manager.updateTangentSmallStrain
        (
            dualTopo, dualGradD, dualGradD0, dt,
            nullptr, &fdC, tangentRequest::fourthOrderFiniteDifference
        );

        if (allLinearElastic)
        {
            // Compare against the analytical tangent rather than a
            // hand-written closed form, so the check is on the finite
            // difference itself, including its Voigt shear convention
            List<mat66> analyticC(nInternalDualFaces);

            manager.updateTangentSmallStrain
            (
                dualTopo, dualGradD, dualGradD0, dt,
                nullptr, &analyticC, tangentRequest::fourthOrder
            );

            scalar maxRelError = 0.0;
            forAll(fdC, dualFaceI)
            {
                const label cellI = dualFaceToCell[dualFaceI];
                const scalar scale = refLambda[cellI] + 2.0*refMu[cellI];

                for (label i = 0; i < 6; ++i)
                {
                    for (label j = 0; j < 6; ++j)
                    {
                        maxRelError =
                            max
                            (
                                maxRelError,
                                mag
                                (
                                    fdC[dualFaceI](i, j)
                                  - analyticC[dualFaceI](i, j)
                                )/scale
                            );
                    }
                }
            }

            reportError("matches the analytical tangent", maxRelError, 1e-6);
        }

        // With no analytical tangent to compare against, require that the
        // finite difference is finite and that repeating it gives exactly the
        // same answer. For a history-dependent law that is a direct check that
        // the perturbed evaluations left no trace in the state
        manager.updateTangentSmallStrain
        (
            dualTopo, dualGradD, dualGradD0, dt,
            nullptr, &fdCagain, tangentRequest::fourthOrderFiniteDifference
        );

        bool finiteAndRepeatable = true;
        forAll(fdC, dualFaceI)
        {
            for (label i = 0; i < 6; ++i)
            {
                for (label j = 0; j < 6; ++j)
                {
                    const scalar a = fdC[dualFaceI](i, j);

                    if (a != a || mag(a) > GREAT)
                    {
                        finiteAndRepeatable = false;
                    }

                    if (a != fdCagain[dualFaceI](i, j))
                    {
                        finiteAndRepeatable = false;
                    }
                }
            }
        }

        report
        (
            "is finite and exactly repeatable",
            finiteAndRepeatable
        );
    }

    // ---------------------------------------------------------------------
    // 8. A tangent query leaves the constitutive state alone
    // ---------------------------------------------------------------------

    Info<< nl << "8. Tangent queries preserve constitutive state" << endl;

    {
        // A law is a function of the kinematics and the OLD-time state, so a
        // stress evaluated straight after a tangent query cannot see anything
        // the query wrote: it recomputes from the same history. The damage a
        // query can do is to the CURRENT-time fields, which endTimeStep()
        // reads and which storeOldTime() promotes to history at the next time
        // step. So the query has to be straddled by a time step to be seen.
        symmTensorField sigmaA(nInternalDualFaces, symmTensor::zero);
        symmTensorField sigmaB(nInternalDualFaces, symmTensor::zero);

        // Establish history at the working strain
        manager.updateStressSmallStrain
        (
            dualTopo, dualGradD, dualGradD0, dt, sigmaA
        );

        // A tangent query at a strain far beyond the working one. Left in the
        // current-time fields, this is what would be committed below
        tensorField wildGradD(nInternalDualFaces, tensor::zero);
        forAll(wildGradD, i)
        {
            wildGradD[i] = 50.0*dualGradD[i];
        }

        scalarField throwaway(nInternalDualFaces, 0.0);
        List<mat66> throwawayC(nInternalDualFaces);

        manager.updateTangentSmallStrain
        (
            dualTopo, wildGradD, dualGradD0, dt,
            &throwaway, nullptr, tangentRequest::scalar
        );

        if (allLinearElastic)
        {
            manager.updateTangentSmallStrain
            (
                dualTopo, wildGradD, dualGradD0, dt,
                nullptr, &throwawayC, tangentRequest::fourthOrder
            );
        }

        manager.updateTangentSmallStrain
        (
            dualTopo, wildGradD, dualGradD0, dt,
            nullptr, &throwawayC, tangentRequest::fourthOrderFiniteDifference
        );

        // Cross a time step, which commits the current-time fields to history,
        // then evaluate the same strain again.
        // setTime rather than operator++: this application is not the solver,
        // so running a time step would execute the case's function objects,
        // and those generally expect a registered solidModel. Only the time
        // index matters here, since that is what the manager keys its
        // old-time rollover on
        runTime.setTime
        (
            runTime.value() + runTime.deltaTValue(),
            runTime.timeIndex() + 1
        );

        manager.updateStressSmallStrain
        (
            dualTopo, dualGradD, dualGradD0, dt, sigmaB
        );

        // A viscoelastic law relaxes, so its stress at a given strain is a
        // function of how much time has passed. Time-independence is the
        // wrong property to demand of it
        if (!allViscoelastic)
        {
        reportError
        (
            "the same strain gives the same stress across a committed "
            "time step",
            relativeDifference(sigmaA, sigmaB),
            1e-12
        );
        }
    }

    // ---------------------------------------------------------------------
    // 9. updateScalarTangent, the cell-centred convenience form
    // ---------------------------------------------------------------------

    Info<< nl << "9. Cell-centred scalar tangent query" << endl;

    {
        volScalarField queriedImpK
        (
            IOobject
            (
                "queriedImpK",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimPressure, 0.0)
        );

        manager.updateScalarTangent
        (
            gradD, gradD0, dt, queriedImpK, tangentRequest::scalar
        );

        reportError
        (
            "agrees with the tangent from the stress update",
            relativeDifference
            (
                Foam::primitiveField(impK), Foam::primitiveField(queriedImpK)
            ),
            1e-15
        );
    }

    // ---------------------------------------------------------------------
    // 10. Misuse guards
    // ---------------------------------------------------------------------

    Info<< nl << "10. Misuse guards" << endl;

    {
        symmTensorField scratch(mesh.nCells(), symmTensor::zero);
        List<mat66> scratchC(mesh.nCells());
        scalarField scratchK(mesh.nCells(), 0.0);

        FatalError.throwExceptions();

        // A cell-centred topology cannot carry a fourth-order tangent
        {
            bool threw = false;
            try
            {
                manager.updateStressSmallStrain
                (
                    cellTopo,
                    Foam::primitiveField(gradD),
                    Foam::primitiveField(gradD0),
                    dt,
                    scratch,
                    nullptr,
                    &scratchC,
                    tangentRequest::fourthOrder
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            report("a fourth-order tangent on cells is rejected", threw);
        }

        // A face-centred topology shares integration points between cells, so
        // the flat-list update, which does not collapse, must refuse it
        {
            const integrationPointTopology& faceTopo =
                manager.topologyFor
                (
                    faceCentredIntegrationPointTopology::typeName
                );

            symmTensorField faceScratch(faceTopo.nIntegrationPoints());
            tensorField faceGradD(faceTopo.nIntegrationPoints(), tensor::zero);

            bool threw = false;
            try
            {
                manager.updateStressSmallStrain
                (
                    faceTopo,
                    faceGradD,
                    faceGradD,
                    dt,
                    faceScratch,
                    nullptr,
                    nullptr,
                    tangentRequest::none
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            // The flat-list update performs no collapse, so a topology whose
            // integration points are shared between cells is refused only
            // when there is more than one law: with a single material no
            // integration point can belong to two materials and there is
            // nothing to collapse
            if (lawEntries.size() > 1)
            {
                report
                (
                    "a flat-list update on shared points is rejected with "
                    "several materials",
                    threw
                );
            }
            else
            {
                report
                (
                    "a flat-list update on shared points is allowed with one "
                    "material",
                    !threw
                );
            }
        }

        // The inputs object's own contract. No law reads a live input yet,
        // so without this the class would ship unexercised
        {
            const scalar dtIn = 0.125;
            mechanicalConstitutiveLawInputs inputs(dtIn);

            report
            (
                "inputs carries the time increment",
                mag(inputs.dt() - dtIn) < SMALL
            );

            report
            (
                "an unsupplied scalar input is absent, not zero",
                !inputs.foundScalar("T") && inputs.findScalar("T") == nullptr
            );

            const scalarField T(3, 300.0);
            inputs.setScalar("T", T);

            report
            (
                "a supplied scalar input is found and readable",
                inputs.foundScalar("T")
             && inputs.findScalar("T") != nullptr
             && mag(inputs.getScalar("T")[1] - 300.0) < SMALL
            );

            // A required input that was never supplied must fail rather than
            // read as zero, which would be a plausible wrong answer
            bool threw = false;
            try
            {
                inputs.getScalar("thisWasNeverSupplied");
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            report
            (
                "a missing required input is rejected, not defaulted",
                threw
            );
        }

        // A tangent request with no storage to put it in
        {
            bool threw = false;
            try
            {
                manager.updateStressSmallStrain
                (
                    cellTopo,
                    Foam::primitiveField(gradD),
                    Foam::primitiveField(gradD0),
                    dt,
                    scratch,
                    nullptr,
                    nullptr,
                    tangentRequest::scalar
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            report("a tangent request with no storage is rejected", threw);
        }

        // A face on a material interface is reached by two laws, so it has
        // two stresses. stressCollapseRule::none gives no rule for combining
        // them and must refuse; a rule that does must be accepted
        if (lawEntries.size() > 1)
        {
            const surfaceTensorField faceGradD(fvc::interpolate(gradD));

            surfaceSymmTensorField faceSigma
            (
                IOobject
                (
                    "faceSigma",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedSymmTensor("0", dimPressure, symmTensor::zero)
            );

            bool threw = false;
            try
            {
                manager.updateStressSmallStrain
                (
                    faceGradD,
                    faceGradD,
                    dt,
                    faceSigma,
                    stressCollapseRule::none
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            report("a shared face with no collapse rule is rejected", threw);

            bool threwWithRule = false;
            try
            {
                manager.updateStressSmallStrain
                (
                    faceGradD,
                    faceGradD,
                    dt,
                    faceSigma,
                    stressCollapseRule::average
                );
            }
            catch (const Foam::error&)
            {
                threwWithRule = true;
            }

            report
            (
                "a shared face with a collapse rule is accepted",
                !threwWithRule
            );
        }

        // A key already in use by a topology of a different type
        {
            bool threw = false;
            try
            {
                manager.registerTopology
                (
                    "testDualFaces",
                    autoPtr<integrationPointTopology>
                    (
                        new faceCentredIntegrationPointTopology(mesh)
                    )
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            report("a clashing registerTopology key is rejected", threw);
        }

        // A topology the manager never took ownership of. Its state is held
        // against the key it was registered under, so one that was never
        // registered has no state to find; it must be refused rather than
        // quietly given a fresh set of history variables
        {
            bool threw = false;

            const faceCentredIntegrationPointTopology foreign(mesh);

            List<tensor> foreignGradD
            (
                foreign.nIntegrationPoints(), tensor::zero
            );
            List<tensor> foreignGradD0(foreignGradD);
            List<symmTensor> foreignSigma
            (
                foreignGradD.size(), symmTensor::zero
            );

            try
            {
                manager.updateStressSmallStrain
                (
                    foreign,
                    foreignGradD,
                    foreignGradD0,
                    dt,
                    foreignSigma
                );
            }
            catch (const Foam::error&)
            {
                threw = true;
            }

            report
            (
                "a topology the manager does not own is rejected", threw
            );
        }

        FatalError.dontThrowExceptions();
    }

    // ---------------------------------------------------------------------
    // The convergence scale a plastic law reports to the manager
    //
    // A law that normalises its Newton tolerance by a scale over its points
    // cannot reduce for itself: the manager evaluates a law only where this
    // rank holds its points, and which points that is depends on the
    // decomposition, so the reductions would not pair up. The law reports a
    // local value and the manager reduces once per law.
    //
    // These check the value reported, and the property that makes reducing
    // it with max the right thing to do. The law is constructed directly
    // rather than through the manager, because no solid model drives the
    // framework at this point in the stack and the arithmetic is what is at
    // issue.
    // ---------------------------------------------------------------------
    {
        // More than two points in the hardening table is what makes the law
        // non-linearly plastic, and so ask for a scale at all
        const fileName tableName
        (
            runTime.constant()/"Test-convergenceScaleHardening"
        );
        {
            OFstream os(tableName);
            os  << "(" << nl
                << "    (0      1.0e9)" << nl
                << "    (0.001  1.1e9)" << nl
                << "    (0.01   1.3e9)" << nl
                << ")" << endl;
        }

        dictionary lawDict;
        lawDict.add("type", word("neoHookeanElasticMisesPlastic"));
        lawDict.add("rho", dimensionedScalar("rho", dimDensity, 7800.0));
        lawDict.add("E", dimensionedScalar("E", dimPressure, 200e9));
        lawDict.add("nu", dimensionedScalar("nu", dimless, 0.3));
        lawDict.add("outOfBounds", word("clamp"));

        // Spelled both ways: OpenFOAM.com's interpolationTable reads "file"
        // and foam-extend's reads "fileName"
        lawDict.add("file", tableName);
        lawDict.add("fileName", tableName);

        autoPtr<mechanicalConstitutiveLaw> plasticPtr
        (
            mechanicalConstitutiveLaw::New(lawDict)
        );

        // Two points: the first undeformed, the second stretched by two
        // along x, so the larger value comes from the second and the checks
        // see which point was taken as well as what was computed
        const label nPts = 2;

        List<tensor> Fs(nPts, tensor::I);
        Fs[1] = tensor(2, 0, 0, 0, 1, 0, 0, 0, 1);

        const List<tensor> F0s(nPts, tensor::I);
        const List<tensor> Finvs(nPts, tensor::I);
        const List<tensor> Finv0s(nPts, tensor::I);

        List<scalar> Js(nPts, 1.0);
        Js[1] = 2.0;

        const List<scalar> J0s(nPts, 1.0);

        mechanicalConstitutiveLawState plasticState(nPts);
        plasticState.symmTensorField0("bEbar") = symmTensor::I;

        labelList allPts(nPts);
        forAll(allPts, i)
        {
            allPts[i] = i;
        }

        const labelList firstPt(1, label(0));
        const labelList secondPt(1, label(1));
        const labelList noPts;

        scalarList scales(4, 0.0);
        const labelList* sets[4] =
        {
            &allPts, &firstPt, &secondPt, &noPts
        };

        for (label s = 0; s < 4; ++s)
        {
            const labelList& pts = *sets[s];

            const UIndirectList<tensor> FView(Fs, pts);
            const UIndirectList<tensor> F0View(F0s, pts);
            const UIndirectList<tensor> FinvView(Finvs, pts);
            const UIndirectList<tensor> Finv0View(Finv0s, pts);
            const UIndirectList<scalar> JView(Js, pts);
            const UIndirectList<scalar> J0View(J0s, pts);

            const finiteStrainMechanicalConstitutiveLawKinematics kin
            (
                FView, F0View, JView, J0View, FinvView, Finv0View
            );

            scales[s] = plasticPtr->localConvergenceScale(kin, plasticState);
        }

        const scalar scaleAll = scales[0];
        const scalar scaleFirst = scales[1];
        const scalar scaleSecond = scales[2];
        const scalar scaleNone = scales[3];

        // 1. Closed form.
        //
        //    bEbar0, F0 and Finv0 are all the identity, so
        //    relFbar = (J/J0)^(-1/3) F and the quantity is
        //    |symm(relFbar & relFbar^T)|.
        //
        //    Undeformed: relFbar = I, symm(I) = I, |I| = sqrt(3).
        //    Stretched:  relFbar = 2^(-1/3) diag(2,1,1), so
        //                relFbar & relFbar^T = 2^(-2/3) diag(4,1,1), whose
        //                magnitude is 2^(-2/3) sqrt(16 + 1 + 1)
        const scalar expectedFirst = Foam::sqrt(3.0);
        const scalar expectedSecond =
            Foam::pow(2.0, -2.0/3.0)*Foam::sqrt(18.0);

        reportError
        (
            "the plastic convergence scale matches the closed form",
            max
            (
                mag(scaleFirst - expectedFirst)/expectedFirst,
                mag(scaleSecond - expectedSecond)/expectedSecond
            ),
            1e-10
        );

        // The larger value is the stretched point, so this also says which
        // point the maximum was taken from
        report
        (
            "the plastic convergence scale is taken over all the points",
            mag(scaleAll - expectedSecond)/expectedSecond < 1e-10
         && expectedSecond > expectedFirst
        );

        // 2. The scale over a set of points is the largest of the scales over
        //    any partition of it. That is what makes reducing with max across
        //    ranks give the answer a serial run gives, and it is the part a
        //    serial test can still check
        report
        (
            "the plastic convergence scale reduces by max over a partition",
            mag(max(scaleFirst, scaleSecond) - scaleAll) < 1e-10*scaleAll
        );

        // 3. A rank holding none of this law's points reports nothing, rather
        //    than a value that would win the reduction
        report
        (
            "a law with no points reports no convergence scale",
            mag(scaleNone) < SMALL
        );

        // 4. A law that needs no scale reports none, which is what makes
        //    asking every law harmless
        {
            dictionary elasticDict;
            elasticDict.add("type", word("linearElastic"));
            elasticDict.add
            (
                "rho", dimensionedScalar("rho", dimDensity, 7800.0)
            );
            elasticDict.add("E", dimensionedScalar("E", dimPressure, 200e9));
            elasticDict.add("nu", dimensionedScalar("nu", dimless, 0.3));

            autoPtr<mechanicalConstitutiveLaw> elasticPtr
            (
                mechanicalConstitutiveLaw::New(elasticDict)
            );

            const UIndirectList<tensor> FView(Fs, allPts);
            const UIndirectList<tensor> F0View(F0s, allPts);
            const UIndirectList<tensor> FinvView(Finvs, allPts);
            const UIndirectList<tensor> Finv0View(Finv0s, allPts);
            const UIndirectList<scalar> JView(Js, allPts);
            const UIndirectList<scalar> J0View(J0s, allPts);

            const finiteStrainMechanicalConstitutiveLawKinematics kin
            (
                FView, F0View, JView, J0View, FinvView, Finv0View
            );

            mechanicalConstitutiveLawState elasticState(nPts);

            report
            (
                "a law that needs no convergence scale reports none",
                mag(elasticPtr->localConvergenceScale(kin, elasticState))
              < SMALL
            );
        }

        Foam::rm(tableName);
    }

    // ---------------------------------------------------------------------
    // Child states
    //
    // A composite law gives each sub-law a state of its own. These check the
    // three things that make that safe: a child is sized like its parent, the
    // old-time rollover reaches it, and a shadow of a parent presents shadows
    // of the children rather than the children themselves
    // ---------------------------------------------------------------------
    {
        Info<< nl << "Child states" << nl;

        mechanicalConstitutiveLawState parent(4);

        report
        (
            "child is absent until asked for",
            !parent.foundChild("sub")
        );

        mechanicalConstitutiveLawState& sub = parent.child("sub");

        report("child is created on first use", parent.foundChild("sub"));
        report
        (
            "child is sized like its parent",
            sub.size() == parent.size(),
            "got " + Foam::name(sub.size())
        );
        report
        (
            "the same child comes back each time",
            &parent.child("sub") == &sub
        );

        // A child's own history must roll over with its parent's. Both times
        // are created up front, as a law's own initialisation does: the
        // rollover walks the old-time table, so a field with no old-time entry
        // is deliberately not history
        sub.scalarField("h") = 1.0;
        sub.scalarField0("h") = 0.0;
        parent.storeOldTime();
        sub.scalarField("h") = 2.0;

        const mechanicalConstitutiveLawState& csub = sub;

        report
        (
            "storeOldTime reaches the child",
            mag(csub.scalarField0("h")[0] - 1.0) < SMALL
         && mag(csub.scalarField("h")[0] - 2.0) < SMALL,
            "old " + Foam::name(csub.scalarField0("h")[0])
          + ", current " + Foam::name(csub.scalarField("h")[0])
        );

        parent.setSize(6);

        report
        (
            "setSize reaches the child",
            sub.size() == 6,
            "got " + Foam::name(sub.size())
        );

        // A shadow must shadow all the way down. Writing through the shadow's
        // child must leave the real child alone, and reading history through
        // it must give the real child's history
        {
            mechanicalConstitutiveLawState shadow
            (
                parent, mechanicalConstitutiveLawState::SHADOW
            );

            mechanicalConstitutiveLawState& shadowSub = shadow.child("sub");

            report
            (
                "a shadow's child is not the parent's child",
                &shadowSub != &sub
            );

            report("a shadow's child is itself a shadow", shadowSub.isShadow());

            const mechanicalConstitutiveLawState& cShadowSub = shadowSub;

            report
            (
                "a shadow's child reads the real child's history",
                mag(cShadowSub.scalarField0("h")[0] - 1.0) < SMALL,
                "got " + Foam::name(cShadowSub.scalarField0("h")[0])
            );

            shadowSub.scalarField("h") = 99.0;

            report
            (
                "writing through a shadow's child leaves the child alone",
                mag(csub.scalarField("h")[0] - 2.0) < SMALL,
                "got " + Foam::name(csub.scalarField("h")[0])
            );
        }
    }

    // ---------------------------------------------------------------------
    // 11. Gathering a cell field onto integration points
    //
    // A prescribed field or coupling input is given per cell. A point on a
    // face must see the interpolated face value whichever topology holds it:
    // a compact face topology reaches an internal face's points from both
    // cells, and a processor face's from one, so anything built from the
    // cells around a point would differ between serial and parallel runs
    // ---------------------------------------------------------------------
    {
        Info<< nl << "11. Gathering a cell field onto integration points"
            << endl;

        // A field that differs from cell to cell, so that a face value and
        // either cell's value are told apart
        volScalarField src
        (
            IOobject
            (
                "gatherTestField",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        forAll(src, cellI)
        {
            src[cellI] = 1.0 + mesh.C()[cellI].x() + 2*mesh.C()[cellI].y();
        }

        src.correctBoundaryConditions();

        const surfaceScalarField srcf(linearInterpolate(src));

        const labelList allCells(identity(mesh.nCells()));

        // The value a point on face faceI must take: the interpolated face
        // value, or the owner cell's on an empty patch, which holds none
        const label nInternal = mesh.nInternalFaces();
        const polyBoundaryMesh& bm = mesh.boundaryMesh();

        const labelList& own = mesh.faceOwner();

        const auto faceValue = [&](const label faceI) -> scalar
        {
            if (faceI < nInternal)
            {
                return srcf[faceI];
            }

            const label patchI = bm.whichPatch(faceI);

            if (srcf.boundaryField()[patchI].empty())
            {
                return src[own[faceI]];
            }

            return srcf.boundaryField()[patchI][faceI - bm[patchI].start()];
        };

        // A compact face topology with two points on every face
        {
            const labelList sizes(mesh.nFaces(), 2);
            CompactListList<label> rows(sizes);

            for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
            {
                rows(faceI, 0) = 2*faceI;
                rows(faceI, 1) = 2*faceI + 1;
            }

            const compactFaceIntegrationPointTopology topo
            (
                mesh, std::move(rows)
            );

            const labelList ipIDs(identity(topo.nIntegrationPoints()));
            scalarField fld(ipIDs.size(), -GREAT);

            Foam::gatherToIntegrationPoints
            (
                mesh, allCells, src, topo, ipIDs, fld
            );

            scalar maxErr = 0;
            for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
            {
                const scalar expected = faceValue(faceI);
                maxErr = max(maxErr, mag(fld[2*faceI] - expected));
                maxErr = max(maxErr, mag(fld[2*faceI + 1] - expected));
            }

            reduce(maxErr, maxOp<scalar>());

            report
            (
                "compact face points take the interpolated face value",
                maxErr < SMALL,
                "max error " + Foam::name(maxErr)
            );
        }

        // A face-centred topology, whose points are the faces themselves
        {
            const faceCentredIntegrationPointTopology topo(mesh);

            const labelList ipIDs(identity(topo.nIntegrationPoints()));
            scalarField fld(ipIDs.size(), -GREAT);

            Foam::gatherToIntegrationPoints
            (
                mesh, allCells, src, topo, ipIDs, fld
            );

            scalar maxErr = 0;
            forAll(ipIDs, faceI)
            {
                maxErr = max(maxErr, mag(fld[faceI] - faceValue(faceI)));
            }

            reduce(maxErr, maxOp<scalar>());

            report
            (
                "face-centred points take the interpolated face value",
                maxErr < SMALL,
                "max error " + Foam::name(maxErr)
            );
        }

        // A cell-centred topology, whose points are the cells
        {
            const cellCentredIntegrationPointTopology topo(mesh);

            const labelList ipIDs(identity(topo.nIntegrationPoints()));
            scalarField fld(ipIDs.size(), -GREAT);

            Foam::gatherToIntegrationPoints
            (
                mesh, allCells, src, topo, ipIDs, fld
            );

            scalar maxErr = 0;
            forAll(ipIDs, cellI)
            {
                maxErr = max(maxErr, mag(fld[cellI] - src[cellI]));
            }

            reduce(maxErr, maxOp<scalar>());

            report
            (
                "cell-centred points take the cell value",
                maxErr < SMALL,
                "max error " + Foam::name(maxErr)
            );
        }
    }

    // ---------------------------------------------------------------------


    Info<< nl << "========================================================="
        << nl;

    if (nFailed_ == 0)
    {
        Info<< "All mechanicalConstitutiveLaw checks passed" << nl
            << "========================================================="
            << nl << endl;

        Info<< "End\n" << endl;

        return 0;
    }

    Info<< nFailed_ << " mechanicalConstitutiveLaw check(s) FAILED" << nl
        << "========================================================="
        << nl << endl;

    return 1;
}


// ************************************************************************* //
