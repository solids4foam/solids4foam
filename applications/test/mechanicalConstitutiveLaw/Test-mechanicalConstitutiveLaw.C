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
      7. The misuse guards fire: a fourth-order tangent on a topology that
         cannot carry one, a flat-list update on a topology whose integration
         points are shared between cells, a tangent request with no storage,
         and a duplicate registerTopology key.

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

    forAll(lawEntries, lawI)
    {
        const dictionary& lawDict = lawEntries[lawI].dict();

        if (word(lawDict.lookup("type")) != "linearElastic")
        {
            FatalErrorInFunction
                << "This test assumes every material is linearElastic, but "
                << lawEntries[lawI].keyword() << " is of type "
                << word(lawDict.lookup("type")) << exit(FatalError);
        }

        const scalar E = dimensionedScalar(lawDict.lookup("E")).value();
        const scalar nu = dimensionedScalar(lawDict.lookup("nu")).value();

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
        << ", cells: " << mesh.nCells() << endl;

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

        reportError
        (
            "dual-face stress uses the owning cell's material",
            relativeDifference(dualSigma, refDualSigma),
            1e-12
        );
    }

    // ---------------------------------------------------------------------
    // 6. A fourth-order tangent on the dual faces
    // ---------------------------------------------------------------------

    Info<< nl << "6. Fourth-order tangent on the dual faces" << endl;

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

        reportError
        (
            "the tangent matches the closed-form isotropic stiffness",
            maxRelError,
            1e-12
        );
    }

    // ---------------------------------------------------------------------
    // 7. Misuse guards
    // ---------------------------------------------------------------------

    Info<< nl << "7. Misuse guards" << endl;

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

            report("a flat-list update on shared points is rejected", threw);
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
