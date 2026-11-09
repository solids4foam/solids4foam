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
    Test-patchQuadratureValues

Description
    Tests the boundary values that the high-order least-squares schemes take
    at the face quadrature points (patchQuadratureValues), for every
    registered scalar and vector patch field type that fixes the value.

    Each type is constructed on the first non-coupled patch of type patch or
    wall with more than one face. A type must behave as follows, and anything
    else is a failure.
      - A patch field derived from quadratureValuePatchField gives values at
        the face quadrature points of the solid model, and stops with an error
        when the values are requested at different quadrature points.
      - A fixedValue patch field gives its face value at every quadrature
        point when the value is uniform, also when it is uniform only up to
        round-off, and stops with an error when it is not uniform. A patch
        field that ignores the assigned value (e.g. zeroValue) is checked
        against the value it keeps.
      - Any other patch field that fixes the value (e.g. mixed and
        directionMixed) stops with an error.

    Types that cannot be constructed on this patch without a dictionary are
    skipped and counted, as are quadratureValuePatchField types that cannot be
    evaluated without the entries of their dictionary (e.g. fixedRotation
    without its rotation axis). Each type is tested in a child process (POSIX fork),
    because some constructors stop with abort(FatalError), which cannot be
    caught. The case must be a high-order solid case, as for
    Test-highOrderGrad, and the utility runs in serial. It exits with zero when
    every type behaves as expected and with one otherwise.

Author
    Amirhossein Taran, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "compatibilityFunctions.H"
#include "fvMeshQuadrature.H"
#include "patchQuadratureValues.H"
#include "wallPolyPatch.H"
#include "solidModel.H"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

//- Return true if patchQuadratureValues gives values for pf at quad, and false
//  if it stops with an error
//  correctValue is set to whether every value equals the face value of pf
template<class Type>
bool givesValues
(
    const fvPatchField<Type>& pf,
    const fvMeshQuadrature& quad,
    bool& correctValue
)
{
    correctValue = false;

    try
    {
        autoPtr<CompactListList<Type>> valuesPtr =
            patchQuadratureValues(pf, quad);
        auto& values = compactListListCRef(autoPtrRef(valuesPtr));

        correctValue = true;

        forAll(values, faceI)
        {
            forAll(values[faceI], pointI)
            {
                if
                (
                    mag(values[faceI][pointI] - pf[faceI])
                  > 1e-12*mag(pf[faceI])
                )
                {
                    correctValue = false;
                }
            }
        }

        return true;
    }
    catch (const Foam::error&)
    {
        return false;
    }
}


//- Set value on fixedValue patch field pf and check that patchQuadratureValues
//  gives the face values when pf is uniform, and an error otherwise
//  The outcome is appended to description.
template<class Type>
bool checkFixedValue
(
    fvPatchField<Type>& pf,
    const Field<Type>& value,
    const bool valueIsUniform,
    const fvMeshQuadrature& quad,
    const std::string& name,
    std::string& description
)
{
    // Some patch fields (e.g. zeroValue) ignore the assigned value
    pf == value;
    const bool uniform = uniformPatchValue(pf);

    bool correctValue = false;
    const bool given = givesValues(pf, quad, correctValue);

    if (!description.empty())
    {
        description += ", ";
    }

    description += name + " (";

    if (uniform && !valueIsUniform)
    {
        description += "kept uniform by the patch field, ";
    }

    description +=
        std::string
        (
            given ? (correctValue ? "exact values" : "wrong values") : "error"
        )
      + ")";

    return uniform ? (given && correctValue) : !given;
}


//- Construct patch field type on patch and test it, writing the outcome to
//  the file descriptor fd
//  This runs in a child process. It writes "C" once the patch field is
//  constructed. It then writes "N" if the patch field does not fix the value.
//  Otherwise it writes the kind of patch field, "I" (quadratureValuePatchField),
//  "V" (fixedValue) or "O" (other), before evaluating it, and then "P"
//  (passed) or "F" (failed) followed by a description of the outcome.
template<class Type>
void testTypeInChild
(
    const int fd,
    const word& type,
    const fvPatch& patch,
    const DimensionedField<Type, volMesh>& iF,
    const fvMeshQuadrature& quad,
    const fvMeshQuadrature& otherQuad,
    const Type& uniformValue
)
{
    typedef fvPatchField<Type> patchFieldType;

    autoPtr<patchFieldType> pfPtr;

    try
    {
        pfPtr.reset(patchFieldType::New(type, patch, iF).ptr());
    }
    catch (const Foam::error&)
    {
        return;
    }

    patchFieldType& pf = autoPtrRef(pfPtr);

    std::string report("C");

    if (!pf.fixesValue())
    {
        report += "N";
        if (::write(fd, report.data(), report.size()) < 0) {}
        return;
    }

    const bool hasInterface =
        dynamic_cast<const quadratureValuePatchField<Type>*>(&pf);
    const bool isFixedValue = isA<fixedValueFvPatchField<Type>>(pf);

    report += (hasInterface ? "I" : (isFixedValue ? "V" : "O"));
    if (::write(fd, report.data(), report.size()) < 0)
    {
        return;
    }

    const Field<Type> uniformField(patch.size(), uniformValue);

    if (hasInterface)
    {
        pf == uniformField;

        bool unused = false;
        const bool ownGivesValues = givesValues(pf, quad, unused);
        const bool otherGivesValues = givesValues(pf, otherQuad, unused);

        report =
            std::string(ownGivesValues && !otherGivesValues ? "P" : "F")
          + "own quadrature points ("
          + (ownGivesValues ? "values" : "error")
          + "), other quadrature points ("
          + (otherGivesValues ? "values" : "error") + ")";
    }
    else if (isFixedValue)
    {
        // A value that is uniform only up to round-off on every second face
        Field<Type> roundOffField(uniformField);
        for (label faceI = 1; faceI < roundOffField.size(); faceI += 2)
        {
            roundOffField[faceI] = (1.0 + 1e-15)*uniformValue;
        }

        // A value that is not uniform
        Field<Type> nonUniformField(uniformField);
        forAll(nonUniformField, faceI)
        {
            nonUniformField[faceI] = (1.0 + faceI)*uniformValue;
        }

        std::string description;

        bool passed =
            checkFixedValue
            (
                pf, uniformField, true, quad, "uniform", description
            );
        passed =
            checkFixedValue
            (
                pf, roundOffField, true, quad, "uniform to round-off",
                description
            ) && passed;
        passed =
            checkFixedValue
            (
                pf, nonUniformField, false, quad, "non-uniform", description
            ) && passed;

        report = std::string(passed ? "P" : "F") + description;
    }
    else
    {
        pf == uniformField;

        bool unused = false;
        const bool given = givesValues(pf, quad, unused);

        report =
            std::string(given ? "F" : "P")
          + "not a fixedValue patch field ("
          + (given ? "values" : "error") + ")";
    }

    if (::write(fd, report.data(), report.size()) < 0) {}
}


//- Test every registered patch field type for fields of Type that fixes the
//  value on patch, and return the number of types that failed
//  Each type is tested in a child process, because some constructors stop
//  with abort(FatalError) on a patch type they do not support (e.g. a wall
//  function on a non-wall patch), and an abort cannot be caught
template<class Type>
label testPatchFieldTypes
(
    const fvPatch& patch,
    const DimensionedField<Type, volMesh>& iF,
    const fvMeshQuadrature& quad,
    const fvMeshQuadrature& otherQuad,
    const Type& uniformValue
)
{
    const wordList types =
        fvPatchField<Type>::patchConstructorTablePtr_->sortedToc();

    label nFixesValue = 0;
    label nSkipped = 0;
    label nNotEvaluated = 0;
    label nFailed = 0;

    Info<< "Patch field types for " << pTraits<Type>::typeName
        << " fields that fix the value" << nl;

    forAll(types, typeI)
    {
        // Flush, so that the child does not repeat buffered output
        Info<< flush;
        std::cout.flush();
        std::cerr.flush();

        int fd[2];
        if (::pipe(fd) != 0)
        {
            FatalErrorInFunction
                << "Cannot create a pipe" << abort(FatalError);
        }

        const pid_t pid = ::fork();

        if (pid < 0)
        {
            FatalErrorInFunction
                << "Cannot create a child process" << abort(FatalError);
        }

        if (pid == 0)
        {
            // Child process: silence its output, test the type and leave
            // without running any exit handlers
            ::close(fd[0]);
            const int devNull = ::open("/dev/null", O_WRONLY);
            ::dup2(devNull, STDOUT_FILENO);
            ::dup2(devNull, STDERR_FILENO);

            testTypeInChild<Type>
            (
                fd[1], types[typeI], patch, iF, quad, otherQuad, uniformValue
            );

            ::close(fd[1]);
            ::_exit(0);
        }

        // Parent process: read the outcome of the child
        ::close(fd[1]);

        std::string report;
        char buffer[256];
        ssize_t n = 0;
        while ((n = ::read(fd[0], buffer, sizeof(buffer))) > 0)
        {
            report.append(buffer, n);
        }
        ::close(fd[0]);

        int status = 0;
        ::waitpid(pid, &status, 0);

        if (report.empty())
        {
            // The type cannot be constructed on this patch without a
            // dictionary
            ++nSkipped;
            continue;
        }

        if (report.size() < 2 || report[1] == 'N')
        {
            continue;
        }

        ++nFixesValue;

        const char kind = report[1];
        bool passed = false;
        std::string result;

        if (report.size() > 2)
        {
            passed = (report[2] == 'P');
            result = report.substr(3);
        }
        else if (kind == 'I')
        {
            // The patch field stopped (abort) when evaluated: it needs the
            // entries of its dictionary (e.g. the axis of fixedRotation)
            ++nNotEvaluated;
            Info<< "    SKIPPED " << types[typeI] << ": cannot be evaluated "
                << "without a dictionary (stopped with an abort)" << nl;
            continue;
        }
        else if (kind == 'O')
        {
            passed = true;
            result = "not a fixedValue patch field (stopped with an abort)";
        }
        else
        {
            result = "stopped with an abort";
        }

        if (!passed)
        {
            ++nFailed;
        }

        Info<< "    " << (passed ? "PASSED " : "FAILED ") << types[typeI]
            << ": " << result.c_str() << nl;
    }

    Info<< "    " << nFixesValue << " types fix the value and " << nFailed
        << " failed. " << nNotEvaluated << " of these cannot be evaluated "
        << "without a dictionary." << nl
        << "    " << nSkipped << " of " << types.size()
        << " registered types cannot be constructed on patch "
        << patch.name() << " without a dictionary." << nl << endl;

    return nFailed;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "This test runs in serial" << abort(FatalError);
    }

    autoPtr<solidModel> solidPtr = solidModel::New
    (
        runTime,
        dynamicFvMesh::defaultRegion
    );
    solidModel& solid = autoPtrRef(solidPtr);
    const fvMesh& mesh = solid.mesh();

    // The quadrature of the solid model, and one with different face points
    const fvMeshQuadrature& quad =
        solid.displacementLeastSquares().quadrature();
    const fvMeshQuadrature otherQuad
    (
        mesh,
        quad.cellOrder(),
        quad.faceOrder() + 1,
        true
    );

    label patchI = -1;
    forAll(mesh.boundary(), pI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[pI];

        if
        (
            (
                pp.type() == polyPatch::typeName
             || pp.type() == wallPolyPatch::typeName
            )
         && pp.size() > 1
        )
        {
            patchI = pI;
            break;
        }
    }

    if (patchI < 0)
    {
        FatalErrorInFunction
            << "No patch of type patch or wall with more than one face was "
            << "found" << abort(FatalError);
    }

    const fvPatch& patch = mesh.boundary()[patchI];

    volScalarField testScalar
    (
        IOobject
        (
            "testScalar",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Errors are thrown, so that each type can report its outcome
    FatalError.throwExceptions();

    label nFailed = 0;

    nFailed += testPatchFieldTypes<scalar>
    (
        patch, testScalar, quad, otherQuad, 2.5
    );

    nFailed += testPatchFieldTypes<vector>
    (
        patch, solid.solutionD(), quad, otherQuad, vector(1, -2, 3)
    );

    FatalError.dontThrowExceptions();

    Info<< "Overall result: " << (nFailed ? "FAILED" : "PASSED") << nl
        << endl;

    return nFailed ? 1 : 0;
}


// ************************************************************************* //
