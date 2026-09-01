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
         the dual faces uses the material of the owning primary cell.
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
         could claim them, a tangent request with no storage, and a duplicate
         registerTopology key.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "mechanicalConstitutiveLawManager.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "mat66.H"
#include "Switch.H"
#include "IOmanip.H"
#include "compatibilityFunctions.H"

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
//  zero still reports a relative error of order one rather than of order 1/SMALL
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
    }

    forAll(lawEntries, lawI)
    {
        const dictionary& lawDict = lawEntries[lawI].dict();

        if (!allLinearElastic && !allNeoHookean)
        {
            continue;
        }

        // A law may be given as E and nu or as mu and K, so accept either
        // here too rather than assuming the first form
        scalar E = 0.0;
        scalar nu = 0.0;

        if (lawDict.found("mu") && lawDict.found("K"))
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

    // ---------------------------------------------------------------------
    // Finite-strain finite-difference tangent
    //
    // Run first, and on its own, because a finite-strain law such as
    // neoHookeanElastic implements no small-strain evaluation: every check
    // below would fatal on it
    // ---------------------------------------------------------------------

    if (allNeoHookean)
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
        // not on any particular strain state
        const tensor gradDSmall
        (
            1e-4,  0.5e-4, 0.0,
            0.5e-4, -0.7e-4, 0.0,
            0.0,    0.0,   0.3e-4
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

        Info<< nl
            << "============================================================"
            << nl;

        if (nFailed_ == 0)
        {
            Info<< "All mechanicalConstitutiveLaw checks passed" << nl
                << "============================================================"
                << nl << nl << "End" << nl << endl;
            return 0;
        }

        Info<< nFailed_ << " mechanicalConstitutiveLaw check(s) FAILED" << nl
            << "============================================================"
            << nl << endl;

        return 1;
    }

    // ---------------------------------------------------------------------
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
    // integration points
    const label nPerCell = 3;
    const label nInternalDualFaces = nPerCell*mesh.nCells();
    const label nBoundaryDualFaces = 5;

    labelList dualFaceToCell(nInternalDualFaces + nBoundaryDualFaces, 0);
    for (label i = 0; i < nInternalDualFaces; ++i)
    {
        dualFaceToCell[i] = i % mesh.nCells();
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

        reportError
        (
            "the same strain gives the same stress across a committed "
            "time step",
            relativeDifference(sigmaA, sigmaB),
            1e-12
        );
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

        FatalError.dontThrowExceptions();
    }

    // ---------------------------------------------------------------------

    Info<< nl << "============================================================"
        << nl;

    if (nFailed_ == 0)
    {
        Info<< "All mechanicalConstitutiveLaw checks passed" << nl
            << "============================================================"
            << nl << endl;

        Info<< "End\n" << endl;

        return 0;
    }

    Info<< nFailed_ << " mechanicalConstitutiveLaw check(s) FAILED" << nl
        << "============================================================"
        << nl << endl;

    return 1;
}



// ************************************************************************* //
