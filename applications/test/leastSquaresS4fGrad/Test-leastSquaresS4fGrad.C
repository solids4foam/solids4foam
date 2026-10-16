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
    Test-leastSquaresS4fGrad

Description
    Tests scalar-field gradients against analytical fields on the case mesh.

    Note: Underlying fields are not symmetric, i.e. this means that this
          utility can't be used with cases with symmetry boundary conditions.

Author
    Ivan Batistic, UCD
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "boolIOList.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchFields.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedValueFvPatchFields.H"
#include "IOmanip.H"
#include "processorFvPatch.H"
#include "compatibilityFunctions.H"

#include <string>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct scalarGradTest
{
    word name;
    scalar (*value)(const vector&);
    vector (*grad)(const vector&);
};


struct errorStats
{
    scalar L2;
    scalar LInf;
    scalar boundaryL2;
    scalar boundaryLInf;
    label nBoundaryCells;
};


scalar mathematicalPi()
{
#ifdef OPENFOAM_NOT_EXTEND
    return constant::mathematical::pi;
#else
    return mathematicalConstant::pi;
#endif
}


bool optionFound(const argList& args, const word& opt)
{
#ifdef OPENFOAM_COM
    return args.found(opt);
#else
    return args.optionFound(opt);
#endif
}


string paddingString(const label width)
{
    return string
    (
        std::string
        (
            static_cast<std::string::size_type>(max(width, label(0))),
            ' '
        )
    );
}


scalar phiXValue(const vector& x)
{
    return x.x();
}


vector phiXGrad(const vector&)
{
    return vector(1.0, 0.0, 0.0);
}


scalar phiLinearValue(const vector& x)
{
    return 1.0 + x.x() + 2.0*x.y() - 0.5*x.z();
}


vector phiLinearGrad(const vector&)
{
    return vector(1.0, 2.0, -0.5);
}


scalar phiQuadraticValue(const vector& x)
{
    return sqr(x.x()) + 0.5*sqr(x.y()) + 0.25*sqr(x.z());
}


vector phiQuadraticGrad(const vector& x)
{
    return vector(2.0*x.x(), x.y(), 0.5*x.z());
}


scalar phiQuadraticCrossValue(const vector& x)
{
    return x.x()*x.y() + 0.25*x.y()*x.z();
}


vector phiQuadraticCrossGrad(const vector& x)
{
    return vector(x.y(), x.x() + 0.25*x.z(), 0.25*x.y());
}


scalar phiTrigonometricValue(const vector& x)
{
    const scalar pi(mathematicalPi());

    return
        Foam::sin(2.0*pi*x.x())*Foam::cos(3.0*pi*x.y())
      + 0.25*Foam::sin(pi*(x.x() + x.y() + x.z()))
      + 0.1*x.x()*x.y()*x.z();
}


vector phiTrigonometricGrad(const vector& x)
{
    const scalar pi(mathematicalPi());

    const scalar common =
        0.25*pi*Foam::cos(pi*(x.x() + x.y() + x.z()));

    return vector
    (
        2.0*pi*Foam::cos(2.0*pi*x.x())*Foam::cos(3.0*pi*x.y())
      + common
      + 0.1*x.y()*x.z(),
       -3.0*pi*Foam::sin(2.0*pi*x.x())*Foam::sin(3.0*pi*x.y())
      + common
      + 0.1*x.x()*x.z(),
        common
      + 0.1*x.x()*x.y()
    );
}


vector activeGradient
(
    const vector& grad,
    const Vector<label>& geometricD
)
{
    vector activeGrad(grad);

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (geometricD[cmpt] <= 0)
        {
            activeGrad[cmpt] = 0.0;
        }
    }

    return activeGrad;
}


boolList boundaryOwnerCells(const fvMesh& mesh)
{
    boolList isBoundaryCell(mesh.nCells(), false);

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if
        (
            isA<emptyFvPatch>(patch)
         || isA<processorFvPatch>(patch)
        )
        {
            continue;
        }

        const labelUList& faceCells = patch.faceCells();
        forAll(faceCells, faceI)
        {
            isBoundaryCell[faceCells[faceI]] = true;
        }
    }

    return isBoundaryCell;
}


volScalarField makeField
(
    const scalarGradTest& test,
    const fvMesh& mesh,
    const Time& runTime
)
{
    wordList patchTypes
    (
        mesh.boundary().size(),
        fixedValueFvPatchScalarField::typeName
    );

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (patch.coupled())
        {
            patchTypes[patchI] = patch.type();
        }
        else if (isA<emptyFvPatch>(patch))
        {
            patchTypes[patchI] = emptyFvPatchScalarField::typeName;
        }
    }

    volScalarField field
    (
        IOobject
        (
            test.name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        patchTypes
    );

    scalarField& fieldI = primitiveFieldRef(field);
    const volVectorField& C = mesh.C();

    forAll(fieldI, cellI)
    {
        fieldI[cellI] = test.value(C[cellI]);
    }

    forAll(field.boundaryField(), patchI)
    {
        fvPatchScalarField& patchField = boundaryFieldRef(field)[patchI];

        if
        (
            patchField.coupled()
         || isA<emptyFvPatchScalarField>(patchField)
        )
        {
            continue;
        }

        scalarField values(patchField.size(), 0.0);
        const vectorField& Cf = patchField.patch().Cf();

        forAll(values, faceI)
        {
            values[faceI] = test.value(Cf[faceI]);
        }

        patchField == values;
    }

    field.correctBoundaryConditions();

    return field;
}


volVectorField makeExactGradField
(
    const scalarGradTest& test,
    const fvMesh& mesh,
    const Time& runTime
)
{
    volVectorField exactGrad
    (
        IOobject
        (
            "gradExact_" + test.name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );

    vectorField& exactGradI = primitiveFieldRef(exactGrad);
    const volVectorField& C = mesh.C();
    const Vector<label> geometricD(mesh.geometricD());

    forAll(exactGradI, cellI)
    {
        exactGradI[cellI] = activeGradient(test.grad(C[cellI]), geometricD);
    }

    forAll(exactGrad.boundaryField(), patchI)
    {
        fvPatchVectorField& patchField =
            boundaryFieldRef(exactGrad)[patchI];

        if
        (
            patchField.coupled()
         || isA<emptyFvPatchVectorField>(patchField)
        )
        {
            continue;
        }

        vectorField values(patchField.size(), vector::zero);
        const vectorField& Cf = patchField.patch().Cf();

        forAll(values, faceI)
        {
            values[faceI] =
                activeGradient(test.grad(Cf[faceI]), geometricD);
        }

        patchField == values;
    }

    return exactGrad;
}


boolIOList makeUseBoundaryFaceValues
(
    const volScalarField& field,
    const Time& runTime,
    const volVectorField& D
)
{
    const fvMesh& mesh = field.mesh();
    boolList useBoundaryFaceValuesList(mesh.boundary().size(), false);

    forAll(mesh.boundary(), patchI)
    {
        if
        (
            isA<fixedDisplacementFvPatchVectorField>
            (
                D.boundaryField()[patchI]
            )
        )
        {
            useBoundaryFaceValuesList[patchI] = true;
        }
    }

    return boolIOList
    (
        IOobject
        (
            "useBoundaryFaceValues_" + field.name(),
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        useBoundaryFaceValuesList
    );
}


errorStats calcErrorStats
(
    const volScalarField& error,
    const boolList& isBoundaryCell
)
{
    scalar LInf = 0.0;
    scalar sumSqrError = 0.0;
    label nCells = 0;

    scalar boundaryLInf = 0.0;
    scalar boundarySumSqrError = 0.0;
    label nBoundaryCells = 0;

#ifdef OPENFOAM_NOT_EXTEND
    const scalarField& errorI = error.primitiveField();
#else
    const scalarField& errorI = error.internalField();
#endif

    forAll(errorI, cellI)
    {
        const scalar e = errorI[cellI];

        LInf = max(LInf, mag(e));
        sumSqrError += magSqr(e);
        nCells++;

        if (isBoundaryCell[cellI])
        {
            boundaryLInf = max(boundaryLInf, mag(e));
            boundarySumSqrError += magSqr(e);
            nBoundaryCells++;
        }
    }

    reduce(LInf, maxOp<scalar>());
    reduce(sumSqrError, sumOp<scalar>());
    reduce(nCells, sumOp<label>());

    reduce(boundaryLInf, maxOp<scalar>());
    reduce(boundarySumSqrError, sumOp<scalar>());
    reduce(nBoundaryCells, sumOp<label>());

    errorStats stats;
    stats.L2 =
        Foam::sqrt(sumSqrError/max(scalar(nCells), VSMALL));
    stats.LInf = LInf;
    stats.boundaryL2 =
        Foam::sqrt
        (
            boundarySumSqrError/max(scalar(nBoundaryCells), VSMALL)
        );
    stats.boundaryLInf = boundaryLInf;
    stats.nBoundaryCells = nBoundaryCells;

    return stats;
}


errorStats testGrad
(
    const scalarGradTest& test,
    const fvMesh& mesh,
    const Time& runTime,
    const boolList& isBoundaryCell,
    const volVectorField& D,
    const bool writeFields
)
{
    volScalarField field(makeField(test, mesh, runTime));
    volVectorField exactGrad(makeExactGradField(test, mesh, runTime));
    boolIOList useBoundaryFaceValuesList
    (
        makeUseBoundaryFaceValues
        (
            field,
            runTime,
            D
        )
    );

    tmp<volVectorField> tGrad(fvc::grad(field));
    volVectorField grad
    (
        IOobject
        (
            "grad_" + test.name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tGrad()
    );

    volScalarField error
    (
        IOobject
        (
            "magGradError_" + test.name,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(grad - exactGrad)
    );

    if (writeFields)
    {
        field.write();
        grad.write();
        exactGrad.write();
        error.write();
    }

    return calcErrorStats(error, isBoundaryCell);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "writeFields",
        "Write analytical fields, gradients and error fields"
    );
    #include "setRootCase.H"

    const bool writeFields = optionFound(args, "writeFields");

    #include "createTime.H"
    #include "createMesh.H"

    volVectorField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const boolList isBoundaryCell(boundaryOwnerCells(mesh));

    const scalarGradTest tests[] =
    {
        {"phi_x", phiXValue, phiXGrad},
        {"phi_linear", phiLinearValue, phiLinearGrad},
        {"phi_quadratic", phiQuadraticValue, phiQuadraticGrad},
        {"phi_quadratic_cross", phiQuadraticCrossValue, phiQuadraticCrossGrad},
        {"phi_trigonometric", phiTrigonometricValue, phiTrigonometricGrad}
    };
    const label nTests = 5;
    const label fieldWidth = 24;
    const label valueWidth = 16;
    const string headerPadding(paddingString(fieldWidth - 5));

    Info<< nl
        << "Testing fvc::grad for scalar analytical fields" << nl
        << "The active gradScheme is read from system/fvSchemes." << nl
        << nl
        << "field" << headerPadding.c_str()
        << setw(valueWidth) << "L2"
        << setw(valueWidth) << "LInf"
        << setw(valueWidth) << "boundaryL2"
        << setw(valueWidth) << "boundaryLInf"
        << nl;

    for (label testI = 0; testI < nTests; testI++)
    {
        const errorStats stats =
            testGrad
            (
                tests[testI],
                mesh,
                runTime,
                isBoundaryCell,
                D,
                writeFields
            );

        const word& fieldName = tests[testI].name;
        const string fieldPadding
        (
            paddingString(max(fieldWidth - label(fieldName.size()), label(1)))
        );

        Info<< fieldName
            << fieldPadding.c_str()
            << setw(valueWidth) << stats.L2
            << setw(valueWidth) << stats.LInf
            << setw(valueWidth) << stats.boundaryL2
            << setw(valueWidth) << stats.boundaryLInf
            << nl;
    }

    Info<< nl << "End" << nl;

    return 0;
}


// ************************************************************************* //
