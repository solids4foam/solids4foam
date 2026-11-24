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
    checkTaylorGreenVortexReconstruction

Description
    Checks the accuracy of the high-order least-squares reconstruction of the
    decaying Taylor-Green vortex velocity and pressure fields.

    The cell values of U and p are set to the exact cell averages at the
    current time, computed with a high-order cell quadrature. U and p are
    then reconstructed with the scheme in system/reconstructionCheckDict,
    with polynomialOrder set to 1, 2 and 3 in turn, and the reconstructed
    values and gradients are compared with the exact solution at

    - the cell centres (cellValue and cellGrad, weighted by the cell volumes),
    - the internal-face quadrature points (internalFaceValue, taken from both
      sides of each face, and internalFaceGrad, weighted by the quadrature
      weights),
    - the quadrature points of the non-empty boundary faces
      (boundaryFaceValue and boundaryFaceGrad).

    For each measure, the weighted L2 norm and the maximum (LInf) of the error
    magnitude are printed. A k-exact reconstruction gives errors of order k+1
    for the values and k for the gradients; verification/Allverify runs the
    mesh refinement study and checks the observed orders.

    The fixed-value patches of U (the walls) are in the U reconstruction
    stencils, with the exact velocity at the face quadrature points given by
    the decayingTaylorGreenVortexVelocity condition. The pressure walls have
    a fixed gradient, so they are not in the pressure stencils and the
    pressure reconstruction is one-sided at the walls.

Author
    Amirhossein Taran, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "emptyFvPatch.H"
#include "fvMeshQuadrature.H"
#include "leastSquaresScheme.H"
#include "compatibilityFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

const scalar pi = constant::mathematical::pi;

// Reynolds number, as in the decayingTaylorGreenVortexVelocity condition
const scalar Re = 10;

// Pressure offset, as in initialiseTaylorGreenVortex
const scalar pOffset = 0.5;

vector exactU(const point& c, const scalar t)
{
    const scalar A = Foam::exp(-2.0*sqr(pi)*t/Re);

    return A*vector
    (
        Foam::sin(pi*c.x())*Foam::cos(pi*c.y()),
      - Foam::cos(pi*c.x())*Foam::sin(pi*c.y()),
        0
    );
}

// Velocity gradient, with (grad U)_ij = d U_j/d x_i as in OpenFOAM
tensor exactGradU(const point& c, const scalar t)
{
    const scalar A = Foam::exp(-2.0*sqr(pi)*t/Re);
    const scalar sx = Foam::sin(pi*c.x());
    const scalar cx = Foam::cos(pi*c.x());
    const scalar sy = Foam::sin(pi*c.y());
    const scalar cy = Foam::cos(pi*c.y());

    return A*pi*tensor
    (
        cx*cy, sx*sy, 0,
      - sx*sy, -cx*cy, 0,
        0, 0, 0
    );
}

scalar exactP(const point& c, const scalar t)
{
    const scalar A = Foam::exp(-4.0*sqr(pi)*t/Re);

    return
        A*0.25*(Foam::cos(2*pi*c.x()) + Foam::cos(2*pi*c.y()))
      + pOffset;
}

vector exactGradP(const point& c, const scalar t)
{
    const scalar A = Foam::exp(-4.0*sqr(pi)*t/Re);

    return -A*0.5*pi*vector
    (
        Foam::sin(2*pi*c.x()),
        Foam::sin(2*pi*c.y()),
        0
    );
}


//- Weighted L2 norm and maximum of an error magnitude
class errorNorm
{
    // Private Data

        //- Weighted sum of the squared errors
        scalar sumSqr_;

        //- Sum of the weights
        scalar sumWeight_;

        //- Maximum error
        scalar max_;

public:

    // Constructors

        //- Construct null
        errorNorm()
        :
            sumSqr_(0),
            sumWeight_(0),
            max_(0)
        {}


    // Member Functions

        //- Add an error with its weight
        void add(const scalar error, const scalar weight)
        {
            sumSqr_ += weight*sqr(error);
            sumWeight_ += weight;
            max_ = Foam::max(max_, error);
        }

        //- Return the weighted L2 norm
        scalar L2() const
        {
            return sumWeight_ > 0 ? Foam::sqrt(sumSqr_/sumWeight_) : 0;
        }

        //- Return the maximum
        scalar LInf() const
        {
            return max_;
        }
};


void printNorm
(
    const word& fieldName,
    const word& measure,
    const errorNorm& norm
)
{
    Info<< "    " << fieldName << " " << measure
        << " L2 = " << norm.L2()
        << " LInf = " << norm.LInf() << nl;
}


template<class Type, class ExactValue, class ExactGrad>
void checkField
(
    const word& fieldName,
    const leastSquaresScheme& scheme,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ExactValue& exactValue,
    const ExactGrad& exactGrad
)
{
    const fvMesh& mesh = vf.mesh();
    const vectorField& C = mesh.C();
    const scalarField& V = mesh.V();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    // Quadrature points and weights are indexed by global face labels
    auto& faceQuadPoints =
        compactListListCRef(scheme.quadrature().faceQuadPoints());
    auto& faceQuadWeights =
        compactListListCRef(scheme.quadrature().faceQuadWeights());

    errorNorm cellValue;
    errorNorm cellGrad;
    errorNorm internalFaceValue;
    errorNorm internalFaceGrad;
    errorNorm boundaryFaceValue;
    errorNorm boundaryFaceGrad;

    // Cell centres
    const auto tGrad(scheme.grad(vf));
    const auto& grad = tGrad();

    forAll(C, cellI)
    {
        cellValue.add
        (
            mag(scheme.valueAtPoint(vf, cellI, C[cellI]) - exactValue(C[cellI])),
            V[cellI]
        );

        cellGrad.add(mag(grad[cellI] - exactGrad(C[cellI])), V[cellI]);
    }

    // Face quadrature points
    CompactListList<typename outerProduct<vector, Type>::type> faceGrad
    (
        faceQuadPoints.sizes()
    );
    scheme.fGrad(vf, faceGrad);

    for (label faceI = 0; faceI < mesh.nInternalFaces(); ++faceI)
    {
        forAll(faceQuadPoints[faceI], qI)
        {
            const point& q = faceQuadPoints[faceI][qI];
            const scalar w = faceQuadWeights[faceI][qI];
            const Type exact(exactValue(q));

            internalFaceValue.add
            (
                mag(scheme.valueAtPoint(vf, own[faceI], q) - exact), w
            );
            internalFaceValue.add
            (
                mag(scheme.valueAtPoint(vf, nei[faceI], q) - exact), w
            );
            internalFaceGrad.add(mag(faceGrad[faceI][qI] - exactGrad(q)), w);
        }
    }

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];

        if (isA<emptyFvPatch>(patch) || patch.coupled())
        {
            continue;
        }

        const labelUList& faceCells = patch.faceCells();

        forAll(faceCells, i)
        {
            const label faceI = patch.start() + i;

            forAll(faceQuadPoints[faceI], qI)
            {
                const point& q = faceQuadPoints[faceI][qI];
                const scalar w = faceQuadWeights[faceI][qI];

                boundaryFaceValue.add
                (
                    mag
                    (
                        scheme.valueAtPoint(vf, faceCells[i], q)
                      - exactValue(q)
                    ),
                    w
                );
                boundaryFaceGrad.add
                (
                    mag(faceGrad[faceI][qI] - exactGrad(q)), w
                );
            }
        }
    }

    printNorm(fieldName, "cellValue", cellValue);
    printNorm(fieldName, "cellGrad", cellGrad);
    printNorm(fieldName, "internalFaceValue", internalFaceValue);
    printNorm(fieldName, "internalFaceGrad", internalFaceGrad);
    printNorm(fieldName, "boundaryFaceValue", boundaryFaceValue);
    printNorm(fieldName, "boundaryFaceGrad", boundaryFaceGrad);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    leastSquaresScheme::checkMeshSupported(mesh);

    const scalar t = runTime.value();

    IOdictionary checkDict
    (
        IOobject
        (
            "reconstructionCheckDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Set the cell values to the exact cell averages. The cell quadrature is
    // exact for polynomials of degree 6, so its error is far below that of
    // the reconstructions being tested
    {
        const fvMeshQuadrature quad(mesh, 6, 6);
        auto& cellQuadPoints = compactListListCRef(quad.cellQuadPoints());
        auto& cellQuadWeights = compactListListCRef(quad.cellQuadWeights());
        const scalarField& V = mesh.V();

        forAll(U, cellI)
        {
            vector UIntegral = Zero;
            scalar pIntegral = 0;

            forAll(cellQuadPoints[cellI], qI)
            {
                const point& q = cellQuadPoints[cellI][qI];
                const scalar w = cellQuadWeights[cellI][qI];

                UIntegral += w*exactU(q, t);
                pIntegral += w*exactP(q, t);
            }

            U.primitiveFieldRef()[cellI] = UIntegral/V[cellI];
            p.primitiveFieldRef()[cellI] = pIntegral/V[cellI];
        }
    }

    U.correctBoundaryConditions();

    // Fixed-value patches are included in the reconstruction stencils
    boolList UPatchesInStencils(mesh.boundary().size(), false);
    boolList pPatchesInStencils(mesh.boundary().size(), false);

    forAll(mesh.boundary(), patchI)
    {
        UPatchesInStencils[patchI] = U.boundaryField()[patchI].fixesValue();
        pPatchesInStencils[patchI] = p.boundaryField()[patchI].fixesValue();
    }

    Info<< nl << "Number of cells " << mesh.nCells() << nl
        << "Geometric dimensions " << mesh.nGeometricD() << nl
        << "Time " << t << nl;

    for (label k = 1; k <= 3; ++k)
    {
        dictionary schemeDict(checkDict);
        schemeDict.set("polynomialOrder", k);

        const autoPtr<leastSquaresScheme> UScheme
        (
            leastSquaresScheme::New(mesh, UPatchesInStencils, schemeDict)
        );

        const autoPtr<leastSquaresScheme> pScheme
        (
            leastSquaresScheme::New(mesh, pPatchesInStencils, schemeDict)
        );

        Info<< nl << "polynomialOrder " << k << nl;

        checkField
        (
            "U",
            UScheme(),
            U,
            [t](const point& c){ return exactU(c, t); },
            [t](const point& c){ return exactGradU(c, t); }
        );

        checkField
        (
            "p",
            pScheme(),
            p,
            [t](const point& c){ return exactP(c, t); },
            [t](const point& c){ return exactGradP(c, t); }
        );
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
