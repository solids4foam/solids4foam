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
    Test-kExactGrad

Description
    Tests cell-centred and internal-face gradients reconstructed by
    kExactLeastSquares from exact cell averages of linear, quadratic and cubic
    polynomials on a two- or three-dimensional mesh. The reconstruction
    settings are read from constant/solidProperties.

    The utility supports a two-dimensional mesh whose empty direction is z, or a
    three-dimensional mesh. It uses third-order cell quadrature to generate the
    reference averages, independently of the polynomial order selected for the
    reconstruction.

    In 2-D, the manufactured polynomials contain every monomial in x and y up to
    the tested order. In 3-D, the same polynomials are extended with every
    z-dependent monomial, so every monomial in x, y and z through cubic order is
    represented. This exercises all reconstructed derivative components rather
    than testing a 2-D polynomial on a 3-D mesh.

    For every polynomial, the utility prints the maximum `grad()` error over all
    cells and the maximum `fGrad()` error over all internal-face quadrature
    points. A polynomial is expected to be reproduced exactly when its degree is
    no greater than the configured `polynomialOrder`. Higher-degree errors are
    still printed but are marked `REPORTED ONLY` and do not determine the exit
    status.

    The utility exits with zero when all expected-exact tests satisfy the
    built-in tolerances and exits with one otherwise. A failure can indicate a
    reconstruction error or an insufficient or rank-deficient stencil on the
    chosen mesh. It can also indicate inconsistent cell quadrature.
    In particular, the k-exact construction assumes that `mesh.C()` is the
    volume centroid, so the first central moments are zero, and that the cell
    quadrature weights sum to `mesh.V()`.
    Run `Test-cellMoments` on the same mesh to check both conditions.

Author
    Ivan Batistic, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMeshQuadrature.H"
#include "leastSquaresReconstruction.H"

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
    scalar fGrad;

    polynomialErrors()
    :
        grad(0.0),
        fGrad(0.0)
    {}
};


polynomialErrors calculateErrors
(
    const label polynomialOrder,
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

    const CompactListList<point>& cellQuadPoints =
        exactQuadrature.cellQuadPoints();
    const CompactListList<scalar>& cellQuadWeights =
        exactQuadrature.cellQuadWeights();
    const scalarField& volumes = mesh.V();

    // Set the field unknowns to exact cell averages of the polynomial
    forAll(phi, cellI)
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

    const tmp<volVectorField> tGrad = reconstruction.grad(phi);
    const volVectorField& grad = tGrad();
    const vectorField& cellCentres = mesh.C();

    polynomialErrors errors;

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
    }

    const CompactListList<point>& faceQuadPoints =
        reconstruction.quadrature().faceQuadPoints();
    CompactListList<vector> faceGrad(faceQuadPoints.sizes());
    reconstruction.fGrad(phi, faceGrad);

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
        }
    }

    reduce(errors.grad, maxOp<scalar>());
    reduce(errors.fGrad, maxOp<scalar>());

    return errors;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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

    if (reconstruction.type() != "kExactLeastSquares")
    {
        FatalErrorInFunction
            << "Test-kExactGrad requires 'type kExactLeastSquares' in "
            << "highOrderCoeffs.displacement, but selected "
            << reconstruction.type() << abort(FatalError);
    }

    // Use third-order quadrature to calculate exact cell averages for all
    // three test fields, independently of the selected reconstruction order
    const fvMeshQuadrature exactQuadrature(mesh, 3, 2, true);

    const scalar gradTolerance = 1e-8;
    const scalar fGradTolerance = 1e-8;
    bool passed = true;

    Info<< nl
        << "K-exact polynomial cell-average gradient test" << nl
        << "    geometric dimensions  : " << mesh.nGeometricD() << nl
        << "    reconstruction order   : "
        << reconstruction.polynomialOrder() << nl
        << "    grad tolerance          : " << gradTolerance << nl
        << "    fGrad tolerance         : " << fGradTolerance << nl << endl;

    for (label order = 1; order <= 3; ++order)
    {
        const polynomialErrors errors = calculateErrors
        (
            order,
            mesh,
            reconstruction,
            exactQuadrature,
            threeD
        );
        const bool expectedExact = order <= reconstruction.polynomialOrder();
        const bool orderPassed =
            errors.grad <= gradTolerance
         && errors.fGrad <= fGradTolerance;

        if (expectedExact)
        {
            passed = passed && orderPassed;
        }

        Info<< "    polynomial order " << order << nl
            << "        maximum grad error : " << errors.grad << nl
            << "        maximum fGrad error: " << errors.fGrad << nl
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

    Info<< "Overall result: " << (passed ? "PASSED" : "FAILED")
        << nl << endl;

    return passed ? 0 : 1;
}


// ************************************************************************* //
