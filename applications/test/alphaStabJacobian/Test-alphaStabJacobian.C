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
    Test-alphaStabJacobian

Description
    Compares the analytical high-order alpha-stabilisation Jacobian-vector
    product with a centred finite-difference derivative of the explicit
    alpha-stabilisation residual.

    Internal-face and fixed-value-boundary contributions are tested
    separately by applying two face-diffusivity masks. The test uses the
    displacement reconstruction and alphaStab settings selected in
    constant/solidProperties, and therefore applies to movingLeastSquares and
    kExactLeastSquares.

    This initial utility tests serial assembly. Parallel execution is skipped
    until processor-face assembly is implemented in
    hofvm::insertAlphaStabIntoPETScMatrix(). A fixed-value-boundary check is
    reported as skipped when the case has no fixed-value displacement patch.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifdef USE_PETSC
#include "PetscHandles.H"
#endif

#include "fvCFD.H"

#if defined(FOAMEXTEND) || !defined(USE_PETSC)

using namespace Foam;

int main(int argc, char *argv[])
{
    Info<< "Test-alphaStabJacobian: SKIPPED (requires OpenFOAM with PETSc)"
        << nl << endl;

    return 0;
}

#else

#include "compatibilityFunctions.H"
#include "foamPetscSnesHelper.H"
#include "hofvm.H"
#include "leastSquaresScheme.H"
#include "petscErrorHandling.H"
#include "solidModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

struct jacobianError
{
    scalar maximumAbsolute;
    scalar maximumReference;
    scalar relative;
    bool passed;

    jacobianError()
    :
        maximumAbsolute(0.0),
        maximumReference(0.0),
        relative(0.0),
        passed(false)
    {}
};


vectorField explicitStabilisationResidual
(
    const stabilisationModel& stabilisation,
    volVectorField& displacement,
    const volTensorField& dummyGrad,
    const surfaceScalarField& diffusivity,
    const vectorField& values
)
{
    Foam::primitiveFieldRef(displacement) = values;
    displacement.correctBoundaryConditions();

    stabilisation.updateVector(displacement, &dummyGrad);

    const volVectorField& cellStabilisation =
        stabilisation.cellVector(&diffusivity, true);
    const scalarField& volumes = displacement.mesh().V();
    vectorField residual(values.size(), vector::zero);

    forAll(residual, cellI)
    {
        residual[cellI] = volumes[cellI]*cellStabilisation[cellI];
    }

    return residual;
}


vectorField finiteDifferenceJacobianVectorProduct
(
    const stabilisationModel& stabilisation,
    volVectorField& displacement,
    const volTensorField& dummyGrad,
    const surfaceScalarField& diffusivity,
    const vectorField& direction,
    const scalar epsilon
)
{
    const vectorField original(displacement.internalField());
    vectorField plusValues(original);
    vectorField minusValues(original);

    plusValues += epsilon*direction;
    minusValues -= epsilon*direction;

    const vectorField plusResidual = explicitStabilisationResidual
    (
        stabilisation,
        displacement,
        dummyGrad,
        diffusivity,
        plusValues
    );
    const vectorField minusResidual = explicitStabilisationResidual
    (
        stabilisation,
        displacement,
        dummyGrad,
        diffusivity,
        minusValues
    );

    Foam::primitiveFieldRef(displacement) = original;
    displacement.correctBoundaryConditions();

    vectorField product(plusResidual);
    product -= minusResidual;
    product /= 2.0*epsilon;

    return product;
}


vectorField analyticalJacobianVectorProduct
(
    Mat jacobian,
    const foamPetscSnesHelper& petscHelper,
    const leastSquaresScheme& reconstruction,
    const volVectorField& displacement,
    const surfaceScalarField& diffusivity,
    const scalar scaleFactor,
    const vectorField& direction,
    const labelList& components
)
{
    hofvm::insertAlphaStabIntoPETScMatrix
    (
        jacobian,
        petscHelper,
        reconstruction,
        displacement,
        diffusivity,
        scaleFactor
    );

    AssertPETSc(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY));
    AssertPETSc(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY));

    VecHandle directionVec;
    VecHandle productVec;
    AssertPETSc(MatCreateVecs(jacobian, directionVec.put(), productVec.put()));
    AssertPETSc(VecSet(directionVec.get(), 0.0));
    petscHelper.InsertFieldComponents<vector>
    (
        direction,
        directionVec.get(),
        0,
        components
    );

    AssertPETSc(MatMult(jacobian, directionVec.get(), productVec.get()));

    vectorField product(direction.size(), vector::zero);
    petscHelper.ExtractFieldComponents<vector>
    (
        productVec.get(),
        product,
        0,
        components
    );

    return product;
}


jacobianError calculateError
(
    const vectorField& analytical,
    const vectorField& finiteDifference,
    const scalar absoluteTolerance,
    const scalar relativeTolerance
)
{
    jacobianError error;

    forAll(analytical, cellI)
    {
        error.maximumAbsolute = max
        (
            error.maximumAbsolute,
            mag(analytical[cellI] - finiteDifference[cellI])
        );
        error.maximumReference = max
        (
            error.maximumReference,
            mag(finiteDifference[cellI])
        );
    }

    reduce(error.maximumAbsolute, maxOp<scalar>());
    reduce(error.maximumReference, maxOp<scalar>());

    error.relative =
        error.maximumAbsolute/max(error.maximumReference, VSMALL);
    error.passed =
        error.maximumAbsolute
     <= absoluteTolerance + relativeTolerance*error.maximumReference;

    return error;
}


jacobianError runCheck
(
    const stabilisationModel& stabilisation,
    volVectorField& displacement,
    const volTensorField& dummyGrad,
    const leastSquaresScheme& reconstruction,
    const foamPetscSnesHelper& petscHelper,
    Mat jacobian,
    const surfaceScalarField& diffusivity,
    const vectorField& direction,
    const labelList& components,
    const scalar epsilon,
    const scalar absoluteTolerance,
    const scalar relativeTolerance
)
{
    const vectorField finiteDifference =
        finiteDifferenceJacobianVectorProduct
        (
            stabilisation,
            displacement,
            dummyGrad,
            diffusivity,
            direction,
            epsilon
        );
    const vectorField analytical = analyticalJacobianVectorProduct
    (
        jacobian,
        petscHelper,
        reconstruction,
        displacement,
        diffusivity,
        stabilisation.scaleFactor(),
        direction,
        components
    );

    return calculateError
    (
        analytical,
        finiteDifference,
        absoluteTolerance,
        relativeTolerance
    );
}


void printResult(const word& name, const jacobianError& error)
{
    Info<< name << nl
        << "    maximum absolute error : " << error.maximumAbsolute << nl
        << "    maximum reference norm : " << error.maximumReference << nl
        << "    maximum relative error : " << error.relative << nl
        << "    result                 : "
        << (error.passed ? "PASSED" : "FAILED") << nl << endl;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    if (Pstream::parRun())
    {
        Info<< "Test-alphaStabJacobian: SKIPPED in parallel until "
            << "processor-face matrix assembly is implemented"
            << nl << endl;

        return 0;
    }

    autoPtr<solidModel> solidPtr = solidModel::New
    (
        runTime,
        dynamicFvMesh::defaultRegion
    );
    solidModel& solid = solidPtr();
    dynamicFvMesh& mesh = solid.mesh();

    const label nDimensions = mesh.nGeometricD();
    if (nDimensions != 2 && nDimensions != 3)
    {
        FatalErrorInFunction
            << "This test requires a two- or three-dimensional mesh"
            << abort(FatalError);
    }

    if (!solid.highOrderResidual())
    {
        FatalErrorInFunction
            << "Test-alphaStabJacobian requires highOrderResidual true"
            << abort(FatalError);
    }

    const solidModel& constSolid = solid;
    const dictionary& solidProperties = constSolid.solidProperties();
    const word solidModelType(solidProperties.lookup("solidModel"));
    const dictionary& solidModelCoeffs =
        solidProperties.subDict(solidModelType + "Coeffs");
    const dictionary& momentumStabilisationCoeffs =
        solidModelCoeffs.subDict("stabilisation").subDict("momentum");
    autoPtr<stabilisationModel> stabilisationPtr = stabilisationModel::New
    (
        mesh,
        momentumStabilisationCoeffs,
        dimless
    );
    const stabilisationModel& stabilisation = stabilisationPtr();
    if (stabilisation.type() != "alpha")
    {
        FatalErrorInFunction
            << "Test-alphaStabJacobian requires momentum stabilisation type "
            << "alphaStab, but selected " << stabilisation.type()
            << abort(FatalError);
    }

    if (mag(stabilisation.scaleFactor()) < SMALL)
    {
        FatalErrorInFunction
            << "The alphaStab scaleFactor must be non-zero for this test"
            << abort(FatalError);
    }

    volVectorField& displacement = solid.solutionD();
    const leastSquaresScheme& reconstruction =
        solid.displacementLeastSquares();
    const foamPetscSnesHelper* petscHelperPtr =
        dynamic_cast<const foamPetscSnesHelper*>(&solid);

    if (petscHelperPtr == nullptr)
    {
        FatalErrorInFunction
            << "The selected solid model does not provide a PETSc SNES helper"
            << abort(FatalError);
    }

    const foamPetscSnesHelper& petscHelper = *petscHelperPtr;

    labelList components(nDimensions);
    forAll(components, componentI)
    {
        components[componentI] = componentI;
    }

    vectorField direction(mesh.nCells(), vector::zero);
    const vectorField& cellCentres = mesh.C();
    forAll(direction, cellI)
    {
        const scalar index = scalar(cellI + 1);
        direction[cellI] = vector
        (
            0.7 + Foam::sin(0.37*index) + 0.11*cellCentres[cellI].x(),
           -0.4 + Foam::cos(0.23*index) - 0.13*cellCentres[cellI].y(),
            nDimensions == 3
          ? 0.2 + Foam::sin(0.19*index) + 0.17*cellCentres[cellI].z()
          : 0.0
        );
    }

    volTensorField dummyGrad
    (
        IOobject
        (
            "alphaStabJacobianDummyGrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor
        (
            "zero",
            displacement.dimensions()/dimLength,
            tensor::zero
        )
    );

    surfaceScalarField internalDiffusivity
    (
        IOobject
        (
            "alphaStabInternalDiffusivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    );
    Foam::primitiveFieldRef(internalDiffusivity) = 1.0;
    Foam::boundaryFieldRef(internalDiffusivity) = 0.0;

    surfaceScalarField boundaryDiffusivity
    (
        IOobject
        (
            "alphaStabBoundaryDiffusivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    );

    label nFixedValueFaces = 0;
    auto& boundaryDiffusivityBf =
        Foam::boundaryFieldRef(boundaryDiffusivity);
    forAll(displacement.boundaryField(), patchI)
    {
        if (displacement.boundaryField()[patchI].fixesValue())
        {
            boundaryDiffusivityBf[patchI] = 1.0;
            nFixedValueFaces += mesh.boundary()[patchI].size();
        }
    }

    MatHandle internalJacobian;
    hofvm::initialiseJacobian
    (
        internalJacobian.m,
        petscHelper,
        reconstruction,
        displacement,
        nDimensions
    );

    const scalar epsilon = 1e-5;
    const scalar absoluteTolerance = 1e-10;
    const scalar relativeTolerance = 1e-8;
    bool passed = true;

    Info<< nl
        << "High-order alphaStab Jacobian-vector product test" << nl
        << "    reconstruction type : " << reconstruction.type() << nl
        << "    polynomial order    : "
        << reconstruction.polynomialOrder() << nl
        << "    alpha scale factor  : " << stabilisation.scaleFactor() << nl
        << "    finite-difference epsilon : " << epsilon << nl
        << "    absolute tolerance       : " << absoluteTolerance << nl
        << "    relative tolerance       : " << relativeTolerance
        << nl << endl;

    const jacobianError internalError = runCheck
    (
        stabilisation,
        displacement,
        dummyGrad,
        reconstruction,
        petscHelper,
        internalJacobian.get(),
        internalDiffusivity,
        direction,
        components,
        epsilon,
        absoluteTolerance,
        relativeTolerance
    );
    printResult("Internal-face contribution", internalError);
    passed = passed && internalError.passed;

    if (nFixedValueFaces)
    {
        MatHandle boundaryJacobian;
        hofvm::initialiseJacobian
        (
            boundaryJacobian.m,
            petscHelper,
            reconstruction,
            displacement,
            nDimensions
        );

        const jacobianError boundaryError = runCheck
        (
            stabilisation,
            displacement,
            dummyGrad,
            reconstruction,
            petscHelper,
            boundaryJacobian.get(),
            boundaryDiffusivity,
            direction,
            components,
            epsilon,
            absoluteTolerance,
            relativeTolerance
        );
        printResult("Fixed-value-boundary contribution", boundaryError);
        passed = passed && boundaryError.passed;
    }
    else
    {
        Info<< "Fixed-value-boundary contribution" << nl
            << "    result                 : SKIPPED (no fixed-value patch)"
            << nl << endl;
    }

    Info<< "Overall result: " << (passed ? "PASSED" : "FAILED")
        << nl << endl;

    return passed ? 0 : 1;
}


// ************************************************************************* //

#endif
