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
    Test-highOrderGrad

Description
    Tests cell-centred and face-quadrature gradients reconstructed by
    leastSquareScheme (kExactLeastSquares,  movingLeastSquares, ...)
    for linear, quadratic and cubic scalar and vector polynomials on a
    two- or three-dimensional mesh. The reconstruction settings are read from
    constant/solidProperties.

    For kExactLeastSquares, the unknowns are initialized with exact cell
    averages. For movingLeastSquares, they are initialized with polynomial
    values at the cell centres.

    In 2-D, the manufactured polynomials contain every monomial in x and y up to
    the tested order. In 3-D, the same polynomials are extended with every
    z-dependent monomial, so every monomial in x, y and z through cubic order is
    represented. This exercises all reconstructed derivative components rather
    than testing a 2-D polynomial on a 3-D mesh.

    For every polynomial, the utility prints the maximum scalar and vector
    `grad()` errors over all cells and the maximum scalar and vector `fGrad()`
    errors over all internal-face quadrature points. A polynomial is expected
    to be reproduced exactly when its degree is no greater than the configured
    `polynomialOrder`. Higher-degree errors are still printed but are marked
    `REPORTED ONLY` and do not determine the exit status.

    The utility exits with zero when all expected-exact tests satisfy the
    built-in tolerances and exits with one otherwise. A failure can indicate a
    reconstruction error or an insufficient or rank-deficient stencil on the
    chosen mesh. It can also indicate inconsistent cell quadrature.
    In particular, the k-exact construction assumes that `mesh.C()` is the
    volume centroid, so the first central moments are zero, and that the cell
    quadrature weights sum to `mesh.V()`.
    Run `Test-cellMoments` on the same mesh to check both conditions.

    In parallel, the utility also verifies that cell and internal-face
    reconstruction stencils use off-processor cell values. Scalar and vector
    gradients are also checked directly at processor-face quadrature points.

    Scalar and vector gradients are checked at physical traction/Neumann face
    quadrature points in serial and parallel. These patches are identified by
    the fixedGradient base type of the primary displacement boundary field.

    In serial, fixed-displacement boundary tests use patches whose exact
    runtime type is fixedDisplacement in the primary displacement field
    selected by the solid model (D or DD). Derived analytical boundary
    conditions are not overwritten. The utility replaces the selected cell
    and patch values in memory with its manufactured constant or normal
    polynomial; it does not modify the case files.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "compatibilityFunctions.H"
#include "fvMeshQuadrature.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedGradientFvPatchFields.H"
#include "leastSquaresReconstruction.H"
#include "processorPolyPatch.H"
#include "solidModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

scalar polynomialValue
(
    const point& p,
    const label order,
    const bool threeD
)
{
    const scalar x = p.x();
    const scalar y = p.y();
    const scalar z = p.z();

    scalar value = 1.0 + 2.0*x - 3.0*y;

    if (threeD)
    {
        value += 4.0*z;
    }

    if (order >= 2)
    {
        value += 0.5*x*x + 0.7*x*y - 0.2*y*y;

        if (threeD)
        {
            value += 0.23*x*z - 0.31*y*z + 0.4*z*z;
        }
    }

    if (order >= 3)
    {
        value +=
            0.11*x*x*x - 0.13*x*x*y
          + 0.17*x*y*y - 0.19*y*y*y;

        if (threeD)
        {
            value +=
                0.07*x*x*z - 0.09*x*y*z + 0.12*x*z*z
              + 0.14*y*y*z - 0.16*y*z*z + 0.18*z*z*z;
        }
    }

    return value;
}


vector polynomialGrad
(
    const point& p,
    const label order,
    const bool threeD
)
{
    const scalar x = p.x();
    const scalar y = p.y();
    const scalar z = p.z();

    vector grad(2.0, -3.0, threeD ? 4.0 : 0.0);

    if (order >= 2)
    {
        grad.x() += x + 0.7*y;
        grad.y() += 0.7*x - 0.4*y;

        if (threeD)
        {
            grad.x() += 0.23*z;
            grad.y() -= 0.31*z;
            grad.z() += 0.23*x - 0.31*y + 0.8*z;
        }
    }

    if (order >= 3)
    {
        grad.x() += 0.33*x*x - 0.26*x*y + 0.17*y*y;
        grad.y() += -0.13*x*x + 0.34*x*y - 0.57*y*y;

        if (threeD)
        {
            grad.x() += 0.14*x*z - 0.09*y*z + 0.12*z*z;
            grad.y() += -0.09*x*z + 0.28*y*z - 0.16*z*z;
            grad.z() +=
                0.07*x*x - 0.09*x*y + 0.24*x*z
              + 0.14*y*y - 0.32*y*z + 0.54*z*z;
        }
    }

    return grad;
}


struct polynomialErrors
{
    scalar grad;
    scalar vectorGrad;
    scalar fGrad;
    scalar vectorFGrad;
    scalar processorFGrad;
    scalar processorVectorFGrad;
    scalar neumannFGrad;
    scalar neumannVectorFGrad;

    polynomialErrors()
    :
        grad(0.0),
        vectorGrad(0.0),
        fGrad(0.0),
        vectorFGrad(0.0),
        processorFGrad(0.0),
        processorVectorFGrad(0.0),
        neumannFGrad(0.0),
        neumannVectorFGrad(0.0)
    {}
};


struct processorBoundaryCoverage
{
    label remoteValues;
    label cellsUsingRemoteValues;
    label internalFacesUsingRemoteValues;
    label processorFaces;

    processorBoundaryCoverage()
    :
        remoteValues(0),
        cellsUsingRemoteValues(0),
        internalFacesUsingRemoteValues(0),
        processorFaces(0)
    {}
};


processorBoundaryCoverage calculateProcessorBoundaryCoverage
(
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction
)
{
    processorBoundaryCoverage coverage;
    const leastSquaresStencil& stencil = reconstruction.stencil();
    const globalIndex& globalCells = stencil.globalCells();
    const List<labelList>& remoteCells = stencil.remoteCellsPerProc();

    forAll(remoteCells, procI)
    {
        coverage.remoteValues += remoteCells[procI].size();
    }

    auto& cellStencils =
        compactListListCRef(stencil.cellsStencil());
    forAll(cellStencils, cellI)
    {
        forAll(cellStencils[cellI], stencilI)
        {
            if (!globalCells.isLocal(cellStencils[cellI][stencilI]))
            {
                ++coverage.cellsUsingRemoteValues;
                break;
            }
        }
    }

    auto& faceStencils = compactListListCRef
    (
        reconstruction.faceGradStencil()
    );

    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        forAll(faceStencils[faceI], stencilI)
        {
            if (!globalCells.isLocal(faceStencils[faceI][stencilI]))
            {
                ++coverage.internalFacesUsingRemoteValues;
                break;
            }
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& processorPatch =
                refCast<const processorPolyPatch>(patch);

            if (Pstream::myProcNo() < processorPatch.neighbProcNo())
            {
                coverage.processorFaces += patch.size();
            }
        }
    }

    reduce(coverage.remoteValues, sumOp<label>());
    reduce(coverage.cellsUsingRemoteValues, sumOp<label>());
    reduce(coverage.internalFacesUsingRemoteValues, sumOp<label>());
    reduce(coverage.processorFaces, sumOp<label>());

    return coverage;
}


boolList physicalNeumannPatchMask
(
    const fvMesh& mesh,
    const volVectorField& displacement
)
{
    boolList patchMask(mesh.boundary().size(), false);

    forAll(patchMask, patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];
        const fvPatchVectorField& patchField =
            displacement.boundaryField()[patchI];

        patchMask[patchI] =
           !patch.coupled()
         && isA<fixedGradientFvPatchVectorField>(patchField);
    }

    return patchMask;
}


label countPatchQuadraturePoints
(
    const boolList& patchMask,
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction
)
{
    auto& faceQuadPoints = compactListListCRef
    (
        reconstruction.quadrature().faceQuadPoints()
    );
    label nPoints = 0;

    forAll(patchMask, patchI)
    {
        if (!patchMask[patchI])
        {
            continue;
        }

        const polyPatch& patch = mesh.boundaryMesh()[patchI];
        forAll(patch, patchFaceI)
        {
            nPoints += faceQuadPoints[patch.start() + patchFaceI].size();
        }
    }

    reduce(nPoints, sumOp<label>());

    return nPoints;
}


scalar normalPolynomialValue
(
    const point& p,
    const point& origin,
    const vector& normal,
    const label order
)
{
    const scalar s = normal & (p - origin);
    scalar value = 1.0 + 2.0*s;

    if (order >= 2)
    {
        value += 0.5*s*s;
    }
    if (order >= 3)
    {
        value += 0.11*s*s*s;
    }

    return value;
}


vector normalPolynomialGrad
(
    const point& p,
    const point& origin,
    const vector& normal,
    const label order
)
{
    const scalar s = normal & (p - origin);
    scalar derivative = 2.0;

    if (order >= 2)
    {
        derivative += s;
    }
    if (order >= 3)
    {
        derivative += 0.33*s*s;
    }

    return derivative*normal;
}


scalar calculateFixedDisplacementBoundaryError
(
    const label polynomialOrder,
    const bool cellAverageUnknown,
    const label patchI,
    const point& origin,
    const vector& normal,
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction,
    const fvMeshQuadrature& exactQuadrature,
    volVectorField& displacement
)
{
    const vector componentScales(1.0, -0.7, 1.3);
    auto& cellQuadPoints = compactListListCRef
    (
        exactQuadrature.cellQuadPoints()
    );
    auto& cellQuadWeights = compactListListCRef
    (
        exactQuadrature.cellQuadWeights()
    );
    const scalarField& volumes = mesh.V();

    forAll(displacement, cellI)
    {
        if (cellAverageUnknown)
        {
            vector integral = vector::zero;

            forAll(cellQuadPoints[cellI], qpI)
            {
                integral +=
                    cellQuadWeights[cellI][qpI]
                   *normalPolynomialValue
                    (
                        cellQuadPoints[cellI][qpI],
                        origin,
                        normal,
                        polynomialOrder
                    )
                   *componentScales;
            }

            displacement[cellI] = integral/volumes[cellI];
        }
        else
        {
            displacement[cellI] =
                normalPolynomialValue
                (
                    mesh.C()[cellI],
                    origin,
                    normal,
                    polynomialOrder
                )
               *componentScales;
        }
    }

    fvPatchVectorField& fixedPatch =
        boundaryFieldRef(displacement)[patchI];
    const polyPatch& pp = mesh.boundaryMesh()[patchI];
    auto& faceQuadPoints = compactListListCRef
    (
        reconstruction.quadrature().faceQuadPoints()
    );

    forAll(fixedPatch, patchFaceI)
    {
        const label faceI = pp.start() + patchFaceI;
        fixedPatch[patchFaceI] =
            normalPolynomialValue
            (
                faceQuadPoints[faceI][0],
                origin,
                normal,
                polynomialOrder
            )
           *componentScales;
    }

    CompactListList<tensor> faceGrad(faceQuadPoints.sizes());
    reconstruction.fGrad(displacement, faceGrad);

    scalar error = 0.0;

    forAll(pp, patchFaceI)
    {
        const label faceI = pp.start() + patchFaceI;

        forAll(faceQuadPoints[faceI], qpI)
        {
            error = max
            (
                error,
                mag
                (
                    faceGrad[faceI][qpI]
                  - normalPolynomialGrad
                    (
                        faceQuadPoints[faceI][qpI],
                        origin,
                        normal,
                        polynomialOrder
                    )
                   *componentScales
                )
            );
        }
    }

    return error;
}


scalar calculateConstantFixedDisplacementBoundaryError
(
    const boolList& fixedPatchMask,
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction,
    volVectorField& displacement
)
{
    const vector constantValue(1.0, -0.7, 1.3);

    forAll(displacement, cellI)
    {
        displacement[cellI] = constantValue;
    }

    forAll(fixedPatchMask, patchI)
    {
        if (!fixedPatchMask[patchI])
        {
            continue;
        }

        fvPatchVectorField& fixedPatch =
            boundaryFieldRef(displacement)[patchI];

        forAll(fixedPatch, patchFaceI)
        {
            fixedPatch[patchFaceI] = constantValue;
        }
    }

    auto& faceQuadPoints = compactListListCRef
    (
        reconstruction.quadrature().faceQuadPoints()
    );
    CompactListList<tensor> faceGrad(faceQuadPoints.sizes());
    reconstruction.fGrad(displacement, faceGrad);
    scalar error = 0.0;

    forAll(fixedPatchMask, patchI)
    {
        if (!fixedPatchMask[patchI])
        {
            continue;
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        forAll(pp, patchFaceI)
        {
            const label faceI = pp.start() + patchFaceI;

            forAll(faceGrad[faceI], qpI)
            {
                error = max(error, mag(faceGrad[faceI][qpI]));
            }
        }
    }

    return error;
}


void setPolynomialFields
(
    const label polynomialOrder,
    const bool cellAverageUnknown,
    const fvMesh& mesh,
    const fvMeshQuadrature& exactQuadrature,
    const bool threeD,
    volScalarField& phi,
    volVectorField& vectorPhi
)
{
    auto& cellQuadPoints = compactListListCRef
    (
        exactQuadrature.cellQuadPoints()
    );
    auto& cellQuadWeights = compactListListCRef
    (
        exactQuadrature.cellQuadWeights()
    );
    const scalarField& volumes = mesh.V();

    const vector componentScales(1.0, -0.7, 1.3);

    // Set the field unknowns according to the selected reconstruction meaning
    forAll(phi, cellI)
    {
        if (cellAverageUnknown)
        {
            scalar integral = 0.0;

            forAll(cellQuadPoints[cellI], qI)
            {
                integral +=
                    cellQuadWeights[cellI][qI]
                   *polynomialValue
                    (
                        cellQuadPoints[cellI][qI],
                        polynomialOrder,
                        threeD
                    );
            }

            phi[cellI] = integral/volumes[cellI];
        }
        else
        {
            phi[cellI] = polynomialValue
            (
                mesh.C()[cellI],
                polynomialOrder,
                threeD
            );
        }

        vectorPhi[cellI] = phi[cellI]*componentScales;
    }
}


void calculateCellGradErrors
(
    const label polynomialOrder,
    const bool threeD,
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction,
    const volScalarField& phi,
    const volVectorField& vectorPhi,
    polynomialErrors& errors
)
{
    const vector componentScales(1.0, -0.7, 1.3);
    const vectorField& cellCentres = mesh.C();

    const tmp<volVectorField> tGrad = reconstruction.grad(phi);
    const volVectorField& grad = tGrad();
    const tmp<volTensorField> tVectorGrad = reconstruction.grad(vectorPhi);
    const volTensorField& vectorGrad = tVectorGrad();

    forAll(cellCentres, cellI)
    {
        errors.grad =
            max
            (
                errors.grad,
                mag
                (
                    grad[cellI]
                  - polynomialGrad
                    (
                        cellCentres[cellI],
                        polynomialOrder,
                        threeD
                    )
                )
            );

        errors.vectorGrad =
            max
            (
                errors.vectorGrad,
                mag
                (
                    vectorGrad[cellI]
                  - polynomialGrad
                    (
                        cellCentres[cellI],
                        polynomialOrder,
                        threeD
                    )*componentScales
                )
            );
    }
}


void calculateFaceGradErrors
(
    const label polynomialOrder,
    const bool threeD,
    const boolList& physicalNeumannPatches,
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction,
    const volScalarField& phi,
    const volVectorField& vectorPhi,
    polynomialErrors& errors
)
{
    const vector componentScales(1.0, -0.7, 1.3);
    auto& faceQuadPoints = compactListListCRef
    (
        reconstruction.quadrature().faceQuadPoints()
    );
    CompactListList<vector> faceGrad(faceQuadPoints.sizes());
    reconstruction.fGrad(phi, faceGrad);
    CompactListList<tensor> vectorFaceGrad(faceQuadPoints.sizes());
    reconstruction.fGrad(vectorPhi, vectorFaceGrad);

    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        forAll(faceQuadPoints[faceI], qI)
        {
            errors.fGrad =
                max
                (
                    errors.fGrad,
                    mag
                    (
                        faceGrad[faceI][qI]
                      - polynomialGrad
                        (
                            faceQuadPoints[faceI][qI],
                            polynomialOrder,
                            threeD
                        )
                    )
                );

            errors.vectorFGrad =
                max
                (
                    errors.vectorFGrad,
                    mag
                    (
                        vectorFaceGrad[faceI][qI]
                      - polynomialGrad
                        (
                            faceQuadPoints[faceI][qI],
                            polynomialOrder,
                            threeD
                        )*componentScales
                    )
                );
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (!isA<processorPolyPatch>(patch))
        {
            continue;
        }

        forAll(patch, patchFaceI)
        {
            const label faceI = patch.start() + patchFaceI;

            forAll(faceQuadPoints[faceI], qI)
            {
                errors.processorFGrad =
                    max
                    (
                        errors.processorFGrad,
                        mag
                        (
                            faceGrad[faceI][qI]
                          - polynomialGrad
                            (
                                faceQuadPoints[faceI][qI],
                                polynomialOrder,
                                threeD
                            )
                        )
                    );

                errors.processorVectorFGrad =
                    max
                    (
                        errors.processorVectorFGrad,
                        mag
                        (
                            vectorFaceGrad[faceI][qI]
                          - polynomialGrad
                            (
                                faceQuadPoints[faceI][qI],
                                polynomialOrder,
                                threeD
                            )*componentScales
                        )
                    );
            }
        }
    }

    forAll(physicalNeumannPatches, patchI)
    {
        if (!physicalNeumannPatches[patchI])
        {
            continue;
        }

        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        forAll(patch, patchFaceI)
        {
            const label faceI = patch.start() + patchFaceI;

            forAll(faceQuadPoints[faceI], qI)
            {
                errors.neumannFGrad =
                    max
                    (
                        errors.neumannFGrad,
                        mag
                        (
                            faceGrad[faceI][qI]
                          - polynomialGrad
                            (
                                faceQuadPoints[faceI][qI],
                                polynomialOrder,
                                threeD
                            )
                        )
                    );

                errors.neumannVectorFGrad =
                    max
                    (
                        errors.neumannVectorFGrad,
                        mag
                        (
                            vectorFaceGrad[faceI][qI]
                          - polynomialGrad
                            (
                                faceQuadPoints[faceI][qI],
                                polynomialOrder,
                                threeD
                            )*componentScales
                        )
                    );
            }
        }
    }
}


polynomialErrors calculateErrors
(
    const label polynomialOrder,
    const bool cellAverageUnknown,
    const boolList& physicalNeumannPatches,
    const fvMesh& mesh,
    const leastSquaresScheme& reconstruction,
    const fvMeshQuadrature& exactQuadrature,
    const bool threeD
)
{
    volScalarField phi
    (
        IOobject
        (
            "testPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        "zeroGradient"
    );

    volVectorField vectorPhi
    (
        IOobject
        (
            "testVectorPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero),
        "zeroGradient"
    );

    setPolynomialFields
    (
        polynomialOrder,
        cellAverageUnknown,
        mesh,
        exactQuadrature,
        threeD,
        phi,
        vectorPhi
    );

    polynomialErrors errors;

    calculateCellGradErrors
    (
        polynomialOrder,
        threeD,
        mesh,
        reconstruction,
        phi,
        vectorPhi,
        errors
    );
    calculateFaceGradErrors
    (
        polynomialOrder,
        threeD,
        physicalNeumannPatches,
        mesh,
        reconstruction,
        phi,
        vectorPhi,
        errors
    );

    reduce(errors.grad, maxOp<scalar>());
    reduce(errors.vectorGrad, maxOp<scalar>());
    reduce(errors.fGrad, maxOp<scalar>());
    reduce(errors.vectorFGrad, maxOp<scalar>());
    reduce(errors.processorFGrad, maxOp<scalar>());
    reduce(errors.processorVectorFGrad, maxOp<scalar>());
    reduce(errors.neumannFGrad, maxOp<scalar>());
    reduce(errors.neumannVectorFGrad, maxOp<scalar>());

    return errors;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

#ifdef FOAMEXTEND
    autoPtr<dynamicFvMesh> meshPtr = dynamicFvMesh::New
    (
        IOobject
        (
            dynamicFvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
    dynamicFvMesh& mesh = autoPtrRef(meshPtr);
    volVectorField displacement
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
#else
    autoPtr<solidModel> solidPtr = solidModel::New
    (
        runTime,
        dynamicFvMesh::defaultRegion
    );
    solidModel& solid = solidPtr();
    dynamicFvMesh& mesh = solid.mesh();
    volVectorField& displacement = solid.solutionD();
#endif

    if (mesh.nGeometricD() != 2 && mesh.nGeometricD() != 3)
    {
        FatalErrorInFunction
            << "This test requires a two- or three-dimensional mesh"
            << abort(FatalError);
    }

    if (mesh.nGeometricD() == 2 && mesh.solutionD()[vector::Z] != -1)
    {
        FatalErrorInFunction
            << "For a two-dimensional mesh, this test requires z as the "
            << "empty direction" << abort(FatalError);
    }

    const bool threeD = mesh.nGeometricD() == 3;

    IOdictionary solidProperties
    (
        IOobject
        (
            "solidProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word solidModelType(solidProperties.lookup("solidModel"));
    const dictionary& solidModelCoeffs =
        solidProperties.subDict(solidModelType + "Coeffs");

    if (!solidModelCoeffs.found("highOrderCoeffs"))
    {
        FatalErrorInFunction
            << "The " << solidModelType << "Coeffs dictionary does not "
            << "contain highOrderCoeffs" << abort(FatalError);
    }

    const dictionary& highOrderCoeffs =
        solidModelCoeffs.subDict("highOrderCoeffs");

    if (!highOrderCoeffs.found("displacement"))
    {
        FatalErrorInFunction
            << "The highOrderCoeffs dictionary does not contain displacement"
            << abort(FatalError);
    }

    const dictionary& displacementCoeffs =
        highOrderCoeffs.subDict("displacement");
    const boolList includePatchInStencils(mesh.boundary().size(), false);
    const leastSquaresReconstruction& reconstructions =
        leastSquaresReconstruction::New(mesh);
    const leastSquaresScheme& reconstruction = reconstructions.scheme
    (
        "displacement",
        includePatchInStencils,
        displacementCoeffs
    );

    const bool cellAverageUnknown =
        reconstruction.type() == "kExactLeastSquares";

    const boolList physicalNeumannPatches =
        physicalNeumannPatchMask(mesh, displacement);
    const label nNeumannQuadraturePoints = countPatchQuadraturePoints
    (
        physicalNeumannPatches,
        mesh,
        reconstruction
    );

    if
    (
        !cellAverageUnknown
     && reconstruction.type() != "movingLeastSquares"
    )
    {
        FatalErrorInFunction
            << "Test-highOrderGrad supports kExactLeastSquares and "
            << "movingLeastSquares, but selected " << reconstruction.type()
            << abort(FatalError);
    }

    // Use third-order quadrature to calculate exact cell averages for all
    // three polynomial degrees, independently of the reconstruction order
    const fvMeshQuadrature exactQuadrature(mesh, 3, 2, true);

    const scalar gradTolerance = 1e-10;
    const scalar fGradTolerance = 1e-10;
    bool passed = true;

    const processorBoundaryCoverage parallelCoverage =
        calculateProcessorBoundaryCoverage(mesh, reconstruction);

    if (Pstream::parRun())
    {
        const bool parallelCoveragePassed =
            parallelCoverage.remoteValues > 0
         && parallelCoverage.cellsUsingRemoteValues > 0
         && parallelCoverage.internalFacesUsingRemoteValues > 0
         && parallelCoverage.processorFaces > 0;

        passed = passed && parallelCoveragePassed;

        Info<< nl
            << "Processor-boundary stencil coverage" << nl
            << "    exchanged remote cell values      : "
            << parallelCoverage.remoteValues << nl
            << "    cells using remote values          : "
            << parallelCoverage.cellsUsingRemoteValues << nl
            << "    internal faces using remote values : "
            << parallelCoverage.internalFacesUsingRemoteValues << nl
            << "    processor faces tested             : "
            << parallelCoverage.processorFaces << nl
            << "    result                             : "
            << (parallelCoveragePassed ? "PASSED" : "FAILED")
            << nl << endl;
    }

    Info<< nl
        << "Least-squares polynomial reconstruction test" << nl
        << "    reconstruction type    : " << reconstruction.type() << nl
        << "    unknown interpretation : "
        << (cellAverageUnknown ? "cell average" : "cell-centre point value")
        << nl
        << "    geometric dimensions  : " << mesh.nGeometricD() << nl
        << "    reconstruction order   : "
        << reconstruction.polynomialOrder() << nl
        << "    grad tolerance          : " << gradTolerance << nl
        << "    fGrad tolerance         : " << fGradTolerance
        << nl
        << "    traction/Neumann boundary quadrature points: "
        << nNeumannQuadraturePoints
        << nl << endl;

    for (label order = 1; order <= 3; ++order)
    {
        const polynomialErrors errors = calculateErrors
        (
            order,
            cellAverageUnknown,
            physicalNeumannPatches,
            mesh,
            reconstruction,
            exactQuadrature,
            threeD
        );
        const bool expectedExact = order <= reconstruction.polynomialOrder();
        const bool orderPassed =
            errors.grad <= gradTolerance
         && errors.vectorGrad <= gradTolerance
         && errors.fGrad <= fGradTolerance
         && errors.vectorFGrad <= fGradTolerance
         && errors.processorFGrad <= fGradTolerance
         && errors.processorVectorFGrad <= fGradTolerance
         &&
            (
                nNeumannQuadraturePoints == 0
             ||
                (
                    errors.neumannFGrad <= fGradTolerance
                 && errors.neumannVectorFGrad <= fGradTolerance
                )
            );

        if (expectedExact)
        {
            passed = passed && orderPassed;
        }

        Info<< "    polynomial order " << order << nl
            << "        maximum scalar grad error: "
            << errors.grad << nl
            << "        maximum vector grad error: "
            << errors.vectorGrad << nl
            << "        maximum scalar fGrad error: "
            << errors.fGrad << nl
            << "        maximum vector fGrad error: "
            << errors.vectorFGrad << nl
            << "        maximum processor-face scalar fGrad error: "
            << errors.processorFGrad << nl
            << "        maximum processor-face vector fGrad error: "
            << errors.processorVectorFGrad << nl
            << "        maximum traction/Neumann-face scalar fGrad error: "
            << errors.neumannFGrad << nl
            << "        maximum traction/Neumann-face vector fGrad error: "
            << errors.neumannVectorFGrad << nl
            << "        expected exact     : "
            << (expectedExact ? "yes" : "no") << nl;

        if (expectedExact)
        {
            Info<< "        result             : "
                << (orderPassed ? "PASSED" : "FAILED") << nl;
        }
        else
        {
            Info<< "        result             : REPORTED ONLY" << nl;
        }

        Info<< endl;
    }

    if (nNeumannQuadraturePoints == 0)
    {
        Info<< "Physical traction/Neumann boundary gradient test: SKIPPED"
            << nl
            << "    No fixedGradient boundary quadrature points were found"
            << nl << endl;
    }

    if (!Pstream::parRun())
    {
        label fixedPatchI = -1;
        boolList allFixedPatchMask(mesh.boundary().size(), false);

        forAll(displacement.boundaryField(), patchI)
        {
            if
            (
                displacement.boundaryField()[patchI].type()
             == fixedDisplacementFvPatchVectorField::typeName
             && mesh.boundaryMesh()[patchI].size()
            )
            {
                allFixedPatchMask[patchI] = true;

                if (fixedPatchI < 0)
                {
                    fixedPatchI = patchI;
                }
            }
        }

        if (fixedPatchI >= 0)
        {
            autoPtr<leastSquaresScheme> allFixedReconstructionPtr =
                leastSquaresScheme::New
                (
                    mesh,
                    allFixedPatchMask,
                    displacementCoeffs
                );
            const leastSquaresScheme& allFixedReconstruction =
                allFixedReconstructionPtr();
            const scalar constantError =
                calculateConstantFixedDisplacementBoundaryError
                (
                    allFixedPatchMask,
                    mesh,
                    allFixedReconstruction,
                    displacement
                );
            const bool constantPassed = constantError <= fGradTolerance;
            passed = passed && constantPassed;

            Info<< nl
                << "Fixed-displacement constant-field test" << nl
                << "    maximum boundary fGrad error: "
                << constantError << nl
                << "    result                      : "
                << (constantPassed ? "PASSED" : "FAILED")
                << nl << endl;

            boolList singleFixedPatchMask(mesh.boundary().size(), false);
            singleFixedPatchMask[fixedPatchI] = true;

            autoPtr<leastSquaresScheme> fixedReconstructionPtr =
                leastSquaresScheme::New
                (
                    mesh,
                    singleFixedPatchMask,
                    displacementCoeffs
                );
            const leastSquaresScheme& fixedReconstruction =
                fixedReconstructionPtr();
            const polyPatch& fixedPatch =
                mesh.boundaryMesh()[fixedPatchI];
            auto& fixedFaceQuadPoints = compactListListCRef
            (
                fixedReconstruction.quadrature().faceQuadPoints()
            );
            const point origin =
                fixedFaceQuadPoints[fixedPatch.start()][0];
            vector normal = sum(fixedPatch.faceAreas());
            normal /= mag(normal) + VSMALL;

            scalar maxOutOfPlane = 0.0;
            scalar maxDistance = 0.0;

            forAll(fixedPatch, patchFaceI)
            {
                const label faceI = fixedPatch.start() + patchFaceI;

                forAll(fixedFaceQuadPoints[faceI], qpI)
                {
                    const vector delta =
                        fixedFaceQuadPoints[faceI][qpI] - origin;
                    maxOutOfPlane =
                        max(maxOutOfPlane, mag(normal & delta));
                    maxDistance = max(maxDistance, mag(delta));
                }
            }

            const bool planar =
                maxOutOfPlane <= 1e-10*max(scalar(1), maxDistance);

            Info<< nl
                << "Fixed-displacement boundary gradient test" << nl
                << "    patch                   : "
                << fixedPatch.name() << nl;

            if (planar)
            {
                for
                (
                    label order = 1;
                    order <= fixedReconstruction.polynomialOrder();
                    ++order
                )
                {
                    const scalar error =
                        calculateFixedDisplacementBoundaryError
                        (
                            order,
                            cellAverageUnknown,
                            fixedPatchI,
                            origin,
                            normal,
                            mesh,
                            fixedReconstruction,
                            exactQuadrature,
                            displacement
                        );
                    const bool orderPassed = error <= fGradTolerance;
                    passed = passed && orderPassed;

                    Info<< "    polynomial order " << order << nl
                        << "        maximum boundary fGrad error: "
                        << error << nl
                        << "        result             : "
                        << (orderPassed ? "PASSED" : "FAILED")
                        << nl;
                }
            }
            else
            {
                Info<< "    result                  : SKIPPED" << nl
                    << "    reason                  : patch is not planar"
                    << nl;
            }

            Info<< endl;
        }
        else
        {
            Info<< nl
                << "Fixed-displacement boundary gradient test: SKIPPED" << nl
                << "    No patch with runtime type fixedDisplacement was found"
                << nl << endl;
        }
    }
    else
    {
        Info<< nl
            << "Fixed-displacement boundary gradient test: SKIPPED" << nl
            << "    This focused test currently runs in serial" << nl << endl;
    }
    Info<< "Overall result: " << (passed ? "PASSED" : "FAILED")
        << nl << endl;

    return passed ? 0 : 1;
}


// ************************************************************************* //
