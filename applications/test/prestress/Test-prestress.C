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
    Test-prestress

Description
    Checks the prestress switch of the mechanicalConstitutiveLawManager on
    the case mesh, using the case constant/mechanicalProperties dictionary,
    which must hold finite-strain laws.

    A smooth, non-uniform, non-symmetric prestress deformation gradient Fpre
    is written to the start time, and two managers are built from the
    dictionary: one with prestress on, which reads it, and one with prestress
    off. The prestressed manager is given the solver's deformation gradient
    Fstep; the plain one is given Fstep & Fpre, with its inverse and
    determinant taken here from the product rather than composed, so the
    reference does not share the composition it checks.

    The following are checked:
      1. The volField stress, internal and boundary values, matches the
         plain manager's at the composed kinematics.
      2. Fpre has an effect: at Fstep = I the stress is that of
         Fpre, and not zero.
      3. The order of composition is checked: Fpre & Fstep gives a
         different stress, so check 1 would see the wrong order.
      4. The scalar tangent matches the plain manager's.
      5. The surfaceField stress matches the plain manager's at Ff composed
         with Fpre interpolated to faces.
      6. The flat-list update on the cell-centred topology matches the
         volField update.
      7. The guards fire: a small-strain update, a flat-list update on a
         topology other than the cell-centred one, a fourth-order tangent,
         a flat list of the wrong size, asking a manager without prestress for its Fpre, building a
         prestressed manager with a law that cannot be prestressed, and
         reading an Fpre whose patches set their own values.

Author
    Conor Moore, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalConstitutiveLawManager.H"
#include "integrationPointTopologies.H"
#include "mechanicalConstitutiveLawTangentRequest.H"
#include "mat66.H"
#include "Switch.H"
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


//- Largest difference between two lists, relative to the largest magnitude
//  in either of them
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

    return returnReduce(maxDiff, maxOp<scalar>())
        /returnReduce(scale, maxOp<scalar>());
}


//- As above, over the internal values and every patch of two fields
template<class Type, template<class> class PatchField, class GeoMesh>
scalar relativeDifference
(
    const GeometricField<Type, PatchField, GeoMesh>& a,
    const GeometricField<Type, PatchField, GeoMesh>& b
)
{
    scalar err =
        relativeDifference(Foam::primitiveField(a), Foam::primitiveField(b));

    forAll(a.boundaryField(), patchI)
    {
        err = max
        (
            err,
            relativeDifference
            (
                a.boundaryField()[patchI], b.boundaryField()[patchI]
            )
        );
    }

    return err;
}


//- The mesh bounds, set once the mesh is read, so that the test fields
//  below are the same few per cent on a mesh of any size
static boundBox bounds_;


//- A point mapped into the unit cube of the mesh bounds
vector unitCoordinates(const vector& p)
{
    const vector span(bounds_.span());

    return vector
    (
        (p.x() - bounds_.min().x())/max(span.x(), SMALL),
        (p.y() - bounds_.min().y())/max(span.y(), SMALL),
        (p.z() - bounds_.min().z())/max(span.z(), SMALL)
    );
}


//- The solver's deformation gradient at a point: a smooth, non-symmetric
//  stretch and shear of a few per cent
tensor testFstep(const vector& p)
{
    const scalar s = 2e-2;

    const vector q(unitCoordinates(p));

    const scalar x = q.x();
    const scalar y = q.y();
    const scalar z = q.z();

    return I + s*tensor
    (
        0.70 + 0.30*x*x,  0.20 - 0.15*y,    0.10 + 0.05*z,
        0.40 - 0.20*x,   -0.50 + 0.25*y*y,  0.30 - 0.10*z,
        0.15 + 0.10*x*y,  0.25 - 0.05*y*z,  0.60 + 0.20*z*z
    );
}


//- Fpre at a point: larger than the step and with different
//  principal directions, so that the two do not commute
tensor testFpre(const vector& p)
{
    const scalar s = 8e-2;

    const vector q(unitCoordinates(p));

    const scalar x = q.x();
    const scalar y = q.y();
    const scalar z = q.z();

    return I + s*tensor
    (
        0.90 - 0.20*y,    0.60 + 0.10*z,   -0.30 + 0.10*x,
       -0.40 + 0.10*z,    0.50 + 0.20*x*x,  0.70 - 0.20*y,
        0.20 + 0.30*y,   -0.60 + 0.10*x,    0.40 + 0.15*y*z
    );
}


//- A calculated field holding fn at the cell and boundary-face centres
volTensorField pointwiseField
(
    const fvMesh& mesh,
    const word& name,
    tensor (*fn)(const vector&),
    const bool registerObject = false
)
{
    volTensorField fld
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        ),
        mesh,
        dimensionedTensor("I", dimless, I),
        "calculated"
    );

    forAll(fld, cellI)
    {
        fld[cellI] = fn(mesh.C()[cellI]);
    }

    forAll(fld.boundaryField(), patchI)
    {
        fvPatchTensorField& pF = Foam::boundaryFieldRef(fld)[patchI];
        const vectorField& Cf = mesh.C().boundaryField()[patchI];

        forAll(pF, faceI)
        {
            pF[faceI] = fn(Cf[faceI]);
        }
    }

    return fld;
}


//- An unregistered field to hold a stress or tangent
template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> resultField
(
    const fvMesh& mesh,
    const word& name,
    const dimensionSet& dims
)
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>("zero", dims, pTraits<Type>::zero),
            "calculated"
        )
    );
}


//- True if f throws a FatalError
template<class Function>
bool throwsFatal(Function f)
{
    FatalError.throwExceptions();

    bool threw = false;

    try
    {
        f();
    }
    catch (const Foam::error&)
    {
        threw = true;
    }

    FatalError.dontThrowExceptions();

    return threw;
}


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const scalar dt = runTime.deltaTValue();

    bounds_ = boundBox(mesh.points(), true);

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

    dictionary plainDict(mechanicalProperties);
    plainDict.set("prestress", Switch(false));

    dictionary preDict(mechanicalProperties);
    preDict.set("prestress", Switch(true));

    // Written, then released, so that the prestressed manager reads it from
    // the start time as a run would, and registers it under its own name
    {
        volTensorField FpreOut(pointwiseField(mesh, "Fpre", testFpre, true));

        // At full precision, so that the read-back check below is exact
        const unsigned int precision = IOstream::defaultPrecision();
        IOstream::defaultPrecision(17);
        FpreOut.write();
        IOstream::defaultPrecision(precision);
    }

    mechanicalConstitutiveLawManager plain(mesh, plainDict);
    mechanicalConstitutiveLawManager pre(mesh, preDict);

    report("the prestressed manager says so", pre.prestressed());
    report("the plain manager says not", !plain.prestressed());

    const volTensorField& Fpre = pre.Fpre();

    {
        const volTensorField FpreRef(pointwiseField(mesh, "FpreRef", testFpre));

        // The internal values only: a constraint patch, a symmetry plane
        // say, is evaluated from them on reading rather than read
        reportError
        (
            "Fpre is read back as written",
            relativeDifference
            (
                Foam::primitiveField(Fpre), Foam::primitiveField(FpreRef)
            ),
            1e-12
        );
    }

    // The solver's kinematics, and the composed kinematics the plain
    // manager is given, its inverse and determinant taken from the product
    const volTensorField F(pointwiseField(mesh, "F", testFstep));
    const volTensorField Finv("Finv", inv(F));
    const volScalarField J("J", det(F));

    const volTensorField Ftot("Ftot", F & Fpre);
    const volTensorField FtotInv("FtotInv", inv(Ftot));
    const volScalarField Jtot("Jtot", det(Ftot));


    // -------------------------------------------------------------------------
    Info<< nl << "1. volField stress against the composed kinematics" << endl;

    tmp<volSymmTensorField> tSigmaPre
    (
        resultField<symmTensor, fvPatchField, volMesh>
        (
            mesh, "sigmaPre", dimPressure
        )
    );
    volSymmTensorField& sigmaPre = Foam::tmpRef(tSigmaPre);

    tmp<volSymmTensorField> tSigmaRef
    (
        resultField<symmTensor, fvPatchField, volMesh>
        (
            mesh, "sigmaRef", dimPressure
        )
    );
    volSymmTensorField& sigmaRef = Foam::tmpRef(tSigmaRef);

    pre.updateStressFiniteStrain(F, F, Finv, Finv, J, J, dt, sigmaPre);

    plain.updateStressFiniteStrain
    (
        Ftot, Ftot, FtotInv, FtotInv, Jtot, Jtot, dt, sigmaRef
    );

    reportError
    (
        "prestressed stress at Fstep is the plain stress at Fstep & Fpre",
        relativeDifference(sigmaPre, sigmaRef),
        1e-10
    );


    // -------------------------------------------------------------------------
    Info<< nl << "2. Fpre alone carries a stress" << endl;

    {
        const volTensorField Fid(pointwiseField(mesh, "Fid", [](const vector&)
        {
            return tensor(I);
        }));
        const volScalarField Jid("Jid", det(Fid));

        const volTensorField FpreInv("FpreInv", inv(Fpre));
        const volScalarField Jpre("Jpre", det(Fpre));

        tmp<volSymmTensorField> tSigmaIdPre
        (
            resultField<symmTensor, fvPatchField, volMesh>
            (
                mesh, "sigmaIdPre", dimPressure
            )
        );
        tmp<volSymmTensorField> tSigmaIdRef
        (
            resultField<symmTensor, fvPatchField, volMesh>
            (
                mesh, "sigmaIdRef", dimPressure
            )
        );

        pre.updateStressFiniteStrain
        (
            Fid, Fid, Fid, Fid, Jid, Jid, dt, Foam::tmpRef(tSigmaIdPre)
        );

        plain.updateStressFiniteStrain
        (
            Fpre, Fpre, FpreInv, FpreInv, Jpre, Jpre, dt,
            Foam::tmpRef(tSigmaIdRef)
        );

        reportError
        (
            "prestressed stress at Fstep = I is the plain stress at Fpre",
            relativeDifference(tSigmaIdPre(), tSigmaIdRef()),
            1e-10
        );

        const scalar magSigma =
            gMax(mag(Foam::primitiveField(tSigmaIdPre())));

        report
        (
            "and it is not zero",
            magSigma > SMALL,
            "max |sigma| " + Foam::name(magSigma)
        );
    }


    // -------------------------------------------------------------------------
    Info<< nl << "3. The order of composition is visible" << endl;

    {
        const volTensorField Fwrong("Fwrong", Fpre & F);
        const volTensorField FwrongInv("FwrongInv", inv(Fwrong));
        const volScalarField Jwrong("Jwrong", det(Fwrong));

        tmp<volSymmTensorField> tSigmaWrong
        (
            resultField<symmTensor, fvPatchField, volMesh>
            (
                mesh, "sigmaWrong", dimPressure
            )
        );

        plain.updateStressFiniteStrain
        (
            Fwrong, Fwrong, FwrongInv, FwrongInv, Jwrong, Jwrong, dt,
            Foam::tmpRef(tSigmaWrong)
        );

        const scalar diff = relativeDifference(sigmaPre, tSigmaWrong());

        report
        (
            "Fpre & Fstep gives a different stress from Fstep & Fpre",
            diff > 1e-6,
            "relative difference " + Foam::name(diff)
        );
    }


    // -------------------------------------------------------------------------
    Info<< nl << "4. Scalar tangent" << endl;

    {
        tmp<volScalarField> tKPre
        (
            resultField<scalar, fvPatchField, volMesh>
            (
                mesh, "impKPre", dimPressure
            )
        );
        tmp<volScalarField> tKRef
        (
            resultField<scalar, fvPatchField, volMesh>
            (
                mesh, "impKRef", dimPressure
            )
        );

        pre.updateScalarTangentFiniteStrain
        (
            F, F, Finv, Finv, J, J, dt, Foam::tmpRef(tKPre)
        );

        plain.updateScalarTangentFiniteStrain
        (
            Ftot, Ftot, FtotInv, FtotInv, Jtot, Jtot, dt, Foam::tmpRef(tKRef)
        );

        reportError
        (
            "prestressed scalar tangent matches the plain one",
            relativeDifference
            (
                Foam::primitiveField(tKPre()), Foam::primitiveField(tKRef())
            ),
            1e-10
        );
    }


    // -------------------------------------------------------------------------
    Info<< nl << "5. surfaceField stress" << endl;

    {
        const surfaceTensorField Ff("Ff", fvc::interpolate(F));
        const surfaceTensorField Finvf("Finvf", inv(Ff));
        const surfaceScalarField Jf("Jf", det(Ff));

        // The face Fpre the switch is documented to use: the cell
        // Fpre interpolated
        const surfaceTensorField Ftotf("Ftotf", Ff & fvc::interpolate(Fpre));
        const surfaceTensorField FtotInvf("FtotInvf", inv(Ftotf));
        const surfaceScalarField Jtotf("Jtotf", det(Ftotf));

        tmp<surfaceSymmTensorField> tSigmafPre
        (
            resultField<symmTensor, fvsPatchField, surfaceMesh>
            (
                mesh, "sigmafPre", dimPressure
            )
        );
        tmp<surfaceSymmTensorField> tSigmafRef
        (
            resultField<symmTensor, fvsPatchField, surfaceMesh>
            (
                mesh, "sigmafRef", dimPressure
            )
        );

        pre.updateStressFiniteStrain
        (
            Ff, Ff, Finvf, Finvf, Jf, Jf, dt, Foam::tmpRef(tSigmafPre)
        );

        plain.updateStressFiniteStrain
        (
            Ftotf, Ftotf, FtotInvf, FtotInvf, Jtotf, Jtotf, dt,
            Foam::tmpRef(tSigmafRef)
        );

        reportError
        (
            "prestressed face stress matches the plain one at the "
            "interpolated Fpre",
            relativeDifference(tSigmafPre(), tSigmafRef()),
            1e-10
        );
    }


    // -------------------------------------------------------------------------
    Info<< nl << "6. Flat list on the cell-centred topology" << endl;

    const integrationPointTopology& cellTopo =
        pre.topologyFor(cellCentredIntegrationPointTopology::typeName);

    {
        symmTensorField sigmaFlat(mesh.nCells(), symmTensor::zero);

        pre.updateStressFiniteStrain
        (
            cellTopo,
            Foam::primitiveField(F),
            Foam::primitiveField(F),
            Foam::primitiveField(Finv),
            Foam::primitiveField(Finv),
            Foam::primitiveField(J),
            Foam::primitiveField(J),
            dt,
            sigmaFlat
        );

        reportError
        (
            "flat-list stress matches the volField stress",
            relativeDifference(sigmaFlat, Foam::primitiveField(sigmaPre)),
            1e-12
        );
    }


    // -------------------------------------------------------------------------
    Info<< nl << "7. Guards" << endl;

    {
        const volTensorField gradD
        (
            "gradD", F - dimensionedTensor("I", dimless, I)
        );

        tmp<volSymmTensorField> tSigmaSS
        (
            resultField<symmTensor, fvPatchField, volMesh>
            (
                mesh, "sigmaSS", dimPressure
            )
        );

        report
        (
            "a small-strain update is refused",
            throwsFatal([&]()
            {
                pre.updateStressSmallStrain
                (
                    gradD, gradD, dt, Foam::tmpRef(tSigmaSS)
                );
            })
        );

        const integrationPointTopology& faceTopo =
            pre.topologyFor(faceCentredIntegrationPointTopology::typeName);

        const label nFaceIP = faceTopo.nIntegrationPoints();
        const tensorField faceF(nFaceIP, I);
        const scalarField faceJ(nFaceIP, 1.0);
        symmTensorField faceSigma(nFaceIP, symmTensor::zero);

        report
        (
            "a flat-list update off the cell-centred topology is refused",
            throwsFatal([&]()
            {
                pre.updateStressFiniteStrain
                (
                    faceTopo, faceF, faceF, faceF, faceF, faceJ, faceJ, dt,
                    faceSigma
                );
            })
        );

        List<mat66> tangent4(mesh.nCells());
        symmTensorField sigmaFlat(mesh.nCells(), symmTensor::zero);

        report
        (
            "a fourth-order tangent is refused",
            throwsFatal([&]()
            {
                pre.updateTangentFiniteStrain
                (
                    cellTopo,
                    Foam::primitiveField(F),
                    Foam::primitiveField(F),
                    Foam::primitiveField(Finv),
                    Foam::primitiveField(Finv),
                    Foam::primitiveField(J),
                    Foam::primitiveField(J),
                    dt,
                    nullptr,
                    &tangent4,
                    tangentRequest::fourthOrderFiniteDifference
                );
            })
        );

        // A list of the wrong size is refused before it is composed with
        // Fpre, which would read it at every cell
        {
            const label nShort = max(mesh.nCells() - 1, 0);
            const tensorField shortF(nShort, I);
            const scalarField shortJ(nShort, 1.0);
            symmTensorField shortSigma(mesh.nCells(), symmTensor::zero);

            report
            (
                "a flat list of the wrong size is refused",
                throwsFatal([&]()
                {
                    pre.updateStressFiniteStrain
                    (
                        cellTopo, shortF, shortF, shortF, shortF, shortJ,
                        shortJ, dt, shortSigma
                    );
                })
            );
        }

        report
        (
            "a manager without prestress has no Fpre to give",
            throwsFatal([&]() { plain.Fpre(); })
        );

        // A law that has not said it can be prestressed is refused when the
        // manager is built, before Fpre is read
        const dictionary unprestressableDict
        (
            IStringStream
            (
                "prestress yes;"
                "mechanical "
                "("
                "    steel"
                "    {"
                "        type linearElastic;"
                "        rho rho [1 -3 0 0 0 0 0] 7800;"
                "        E E [1 -1 -2 0 0 0 0] 2e11;"
                "        nu nu [0 0 0 0 0 0 0] 0.3;"
                "    }"
                ");"
            )()
        );

        report
        (
            "a law that cannot be prestressed is refused",
            throwsFatal([&]()
            {
                mechanicalConstitutiveLawManager
                    refused(mesh, unprestressableDict);
            })
        );

        // An Fpre whose patches set their own values, here zeroGradient on
        // every patch that is not a constraint one, is refused when read. It
        // replaces the Fpre in the start time, which the managers above have
        // already read
        {
            wordList patchTypes(mesh.boundary().size());

            forAll(mesh.boundary(), patchI)
            {
                const word& type = mesh.boundary()[patchI].type();

                patchTypes[patchI] =
                    polyPatch::constraintType(type) ? type : "zeroGradient";
            }

            volTensorField FpreZeroGrad
            (
                IOobject
                (
                    "Fpre",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedTensor("I", dimless, I),
                patchTypes
            );

            FpreZeroGrad.write();
        }

        report
        (
            "an Fpre with patches that set their own values is refused",
            throwsFatal([&]()
            {
                mechanicalConstitutiveLawManager refused(mesh, preDict);
            })
        );
    }


    // -------------------------------------------------------------------------

    Info<< nl << "========================================================="
        << nl;

    if (nFailed_ == 0)
    {
        Info<< "All prestress checks passed" << nl
            << "========================================================="
            << nl << endl;

        Info<< "End\n" << endl;

        return 0;
    }

    Info<< nFailed_ << " prestress check(s) FAILED" << nl
        << "========================================================="
        << nl << endl;

    return 1;
}


// ************************************************************************* //
