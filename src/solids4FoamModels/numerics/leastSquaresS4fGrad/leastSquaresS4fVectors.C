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

\*---------------------------------------------------------------------------*/

#include "leastSquaresS4fVectors.H"
#include "volFields.H"
#include "symmetryPolyPatch.H"
#include "compatibilityFunctions.H"
#include "cellZoneInterface.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "symmetryPlanePolyPatch.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresS4fVectors, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{
    // The smallest and largest eigenvalues of a symmetric tensor.
    //
    // OpenFOAM.org and foam-extend find eigenvalues as the roots of the
    // characteristic cubic. A repeated eigenvalue is a double root, which
    // rounding can push off the real axis, and both then report it as complex
    // and return zero for it. The empty-direction fill in
    // calcWideStencilVectors gives a repeated eigenvalue in every uniform 2-D
    // cell, so there a healthy stencil would read as rank-deficient. The
    // trigonometric form used here is exact for a symmetric tensor, repeated
    // eigenvalues included, and gives the same answer on every fork
    void minMaxEigenValues
    (
        const symmTensor& T,
        scalar& lambdaMin,
        scalar& lambdaMax
    )
    {
        const scalar q = tr(T)/3.0;

        const scalar p =
            sqrt
            (
                (
                    sqr(T.xx() - q) + sqr(T.yy() - q) + sqr(T.zz() - q)
                  + 2.0*(sqr(T.xy()) + sqr(T.xz()) + sqr(T.yz()))
                )/6.0
            );

        if (p <= SMALL*mag(q) || p < VSMALL)
        {
            // Isotropic, or zero
            lambdaMin = q;
            lambdaMax = q;
            return;
        }

        symmTensor B(T);
        B.xx() -= q;
        B.yy() -= q;
        B.zz() -= q;
        B /= p;

        const scalar phi = acos(max(min(0.5*det(B), 1.0), -1.0))/3.0;

        // The three eigenvalues are q + 2p cos(phi + 2k pi/3), k = 0, 1, 2, and
        // with phi in [0, pi/3] the largest is k = 0 and the smallest k = 1
        lambdaMax = q + 2.0*p*cos(phi);
        lambdaMin =
            q - p*(cos(phi) + sqrt(3.0)*sin(phi));
    }


    // Flag cells whose stencil does not span three directions, or spans one
    // of them poorly.
    //
    // dd is symmetric positive semi-definite, so its eigenvalues are real and
    // non-negative, and the smallest one vanishes exactly when the stencil is
    // confined to a plane or a line. This ratio sets how poorly conditioned a
    // face stencil may be before the cell is given the wider stencil. It is
    // not used to decide how a tensor is inverted, which pseudoInverse does
    // from the eigenvalues alone
    const scalar minEigenRatio = 1e-4;

    // A 2-D or axisymmetric mesh carries no information in its empty
    // direction, and that is not a degeneracy. The unit tensor of the empty
    // directions is added, scaled by the mean eigenvalue, before testing, as
    // inv(symmTensorField) does through safeInv()
    symmTensor emptyDirections(const fvMesh& mesh)
    {
        const Vector<label> geometricD = mesh.geometricD();

        symmTensor emptyDirs(symmTensor::zero);
        emptyDirs.xx() = (geometricD.x() == -1 ? 1.0 : 0.0);
        emptyDirs.yy() = (geometricD.y() == -1 ? 1.0 : 0.0);
        emptyDirs.zz() = (geometricD.z() == -1 ? 1.0 : 0.0);

        return emptyDirs;
    }

    bool illConditioned
    (
        const symmTensor& dd,
        const symmTensor& emptyDirs,
        scalar& lambdaMin,
        scalar& lambdaMax
    )
    {
        const symmTensor ddc(dd + (tr(dd)/3.0)*emptyDirs);
        minMaxEigenValues(ddc, lambdaMin, lambdaMax);

        return lambdaMax < VSMALL || lambdaMin <= minEigenRatio*lambdaMax;
    }


    // The eigen-decomposition of a symmetric tensor by cyclic Jacobi
    // rotations: T = sum_i lambda_i e_i e_i^T. Each rotation zeroes one
    // off-diagonal entry exactly, and for a positive semi-definite tensor
    // the small eigenvalues come out with high relative accuracy, which the
    // closed-form roots of the characteristic cubic cannot give
    void jacobiEigenDecomposition
    (
        const symmTensor& T,
        vector& lambda,
        List<vector>& e
    )
    {
        scalar a[3][3] =
        {
            {T.xx(), T.xy(), T.xz()},
            {T.xy(), T.yy(), T.yz()},
            {T.xz(), T.yz(), T.zz()}
        };

        scalar v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

        for (label sweep = 0; sweep < 50; sweep++)
        {
            const scalar offDiag =
                sqr(a[0][1]) + sqr(a[0][2]) + sqr(a[1][2]);
            const scalar diag = sqr(a[0][0]) + sqr(a[1][1]) + sqr(a[2][2]);

            if (offDiag == 0 || offDiag <= sqr(VSMALL)*diag)
            {
                break;
            }

            for (label p = 0; p < 2; p++)
            {
                for (label q = p + 1; q < 3; q++)
                {
                    // An entry already negligible next to its diagonal is
                    // zeroed rather than rotated away: its rotation angle
                    // would be huge, and its square can overflow
                    if
                    (
                        mag(a[p][q])
                     <= VSMALL*(mag(a[p][p]) + mag(a[q][q]))
                    )
                    {
                        a[p][q] = 0;
                        a[q][p] = 0;
                        continue;
                    }

                    // The tangent of the rotation angle that zeroes a[p][q],
                    // in its small-angle form when the angle is tiny
                    const scalar theta = (a[q][q] - a[p][p])/(2.0*a[p][q]);
                    scalar t = 0;
                    if (mag(theta) > 1.0/ROOTVSMALL)
                    {
                        t = 0.5/theta;
                    }
                    else
                    {
                        t =
                            (theta >= 0 ? 1.0 : -1.0)
                           /(mag(theta) + sqrt(sqr(theta) + 1.0));
                    }
                    const scalar cs = 1.0/sqrt(sqr(t) + 1.0);
                    const scalar sn = t*cs;

                    for (label k = 0; k < 3; k++)
                    {
                        const scalar akp = a[k][p];
                        const scalar akq = a[k][q];
                        a[k][p] = cs*akp - sn*akq;
                        a[k][q] = sn*akp + cs*akq;
                    }

                    for (label k = 0; k < 3; k++)
                    {
                        const scalar apk = a[p][k];
                        const scalar aqk = a[q][k];
                        a[p][k] = cs*apk - sn*aqk;
                        a[q][k] = sn*apk + cs*aqk;
                    }

                    for (label k = 0; k < 3; k++)
                    {
                        const scalar vkp = v[k][p];
                        const scalar vkq = v[k][q];
                        v[k][p] = cs*vkp - sn*vkq;
                        v[k][q] = sn*vkp + cs*vkq;
                    }
                }
            }
        }

        lambda = vector(a[0][0], a[1][1], a[2][2]);

        e.setSize(3);
        for (label i = 0; i < 3; i++)
        {
            e[i] = vector(v[0][i], v[1][i], v[2][i]);
        }
    }
}
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::symmTensor Foam::leastSquaresS4fVectors::pseudoInverse
(
    const symmTensor& dd,
    label& rank
)
{
    // A direction whose eigenvalue is below this fraction of the sum of
    // them is treated as one the stencil does not span: the inverse along it
    // would amplify round-off by its reciprocal
    const scalar rankTol = 1e-12;

    // Sum of the eigenvalues, all non-negative
    const scalar scale = tr(dd);

    if (scale < VSMALL)
    {
        rank = 0;
        return symmTensor::zero;
    }

    vector lambda(vector::zero);
    List<vector> e;
    jacobiEigenDecomposition(dd, lambda, e);

    const scalar tol = rankTol*scale;

    rank = 0;
    symmTensor pinv(symmTensor::zero);

    for (label i = 0; i < 3; i++)
    {
        if (lambda[i] > tol)
        {
            rank++;
            pinv += (1.0/lambda[i])*sqr(e[i]);
        }
    }

    return pinv;
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresS4fVectors::leastSquaresS4fVectors
(
    const fvMesh& mesh,
    const boolList& useBoundaryFaceValues_
)
:
#ifdef OPENFOAM_NOT_EXTEND
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresS4fVectors>(mesh),
#else
    MeshObject<fvMesh, leastSquaresS4fVectors>(mesh),
#endif
    useBoundaryFaceValues_(useBoundaryFaceValues_),
    pVectors_
    (
        IOobject
        (
            "LeastSquaresP",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector("0", dimless/dimLength, vector::zero)
    ),
    nVectors_
    (
        IOobject
        (
            "LeastSquaresN",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector("0", dimless/dimLength, vector::zero)
    )
{
    calcLeastSquaresVectors();
}


Foam::leastSquaresS4fVectors::leastSquaresS4fVectors
(
    const word& objName,
    const fvMesh& mesh,
    const boolList& useBoundaryFaceValues_
)
:
#ifdef OPENFOAM_COM
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresS4fVectors>
    (
        objName, mesh
    ),
#elif defined(OPENFOAM_ORG)
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresS4fVectors>(mesh),
#else
    MeshObject<fvMesh, leastSquaresS4fVectors>(mesh),
#endif
    useBoundaryFaceValues_(useBoundaryFaceValues_),
    pVectors_
    (
        IOobject
        (
            "LeastSquaresP",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    ),
    nVectors_
    (
        IOobject
        (
            "LeastSquaresN",
            mesh.pointsInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    )
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresS4fVectors::~leastSquaresS4fVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresS4fVectors::calcLeastSquaresVectors() const
{
    DebugInFunction
        << "Calculating least square gradient vectors" << nl;

    const fvMesh& mesh = this->mesh();

    // Set local references to mesh data
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceScalarField& magSf = mesh.magSf();

    const Field<bool> interface(cellZoneInterface(mesh));

    // A material interface that runs along a processor boundary is an internal
    // face on neither side, so it has to be skipped here as well or the
    // stencil would reach across it in parallel and not in serial
    const List<boolList> interfaceCoupled(cellZoneInterfaceCoupled(mesh));

    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh.nCells(), symmTensor::zero);
    forAll(owner, facei)
    {
        if (interface[facei])
        {
            // Skip contributions across interfaces
            continue;
        }

        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        symmTensor wdd = (magSf[facei]/magSqr(d))*sqr(d);

        dd[own] += (1 - w[facei])*wdd;
        dd[nei] += w[facei]*wdd;
    }


    auto& pVectorsBf = boundaryFieldRef(pVectors_);

    forAll(pVectorsBf, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.patch().faceCells();
        const vectorField& Cf = p.Cf();

        // Build the d-vectors
        // In OF.com/OF.org, p.delta are the orthogonal components of the real d
        // vectors, so we need to build them ourselves
        //vectorField pd(p.delta());

        if (pw.coupled())
        {
            // Coupled d vectors
            const vectorField pd(p.delta());

            forAll(pd, patchFacei)
            {
                if (interfaceCoupled[patchi][patchFacei])
                {
                    // Skip contributions across interfaces
                    continue;
                }

                const vector& d = pd[patchFacei];

                dd[faceCells[patchFacei]] +=
                    ((1 - pw[patchFacei])*pMagSf[patchFacei]/magSqr(d))*sqr(d);
            }
        }
        else if
        (
            isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchi])
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchi])
#endif
        )
        {
            // Treat symmetry planes consistently with internal faces
            // Use the mirrored face-cell values rather than the patch face
            // values
            // See https://doi.org/10.1080/10407790.2022.2105073
            const vectorField nHat(p.nf());
            forAll(nHat, patchFacei)
            {
                // Vector from the face-cell centre to the mirrored face-cell
                // centred
                //const vector& d = pd[patchFacei];
                const vector d =
                    transform
                    (
                        I - 2.0*sqr(nHat[patchFacei]),
                        C[faceCells[patchFacei]]
                    )
                  - C[faceCells[patchFacei]];

                // Multiply by 0.5 to be consistent with the internal field
                // where we multiply by the interpolation weight
                dd[faceCells[patchFacei]] +=
                    0.5*(pMagSf[patchFacei]/magSqr(d))*sqr(d);
            }
        }
        else if (useBoundaryFaceValues_[patchi])
        {
            vectorField pd(p.size(), vector::zero);
            forAll(pd, faceI)
            {
                pd[faceI] = Cf[faceI] - C[faceCells[faceI]];
            }

            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                dd[faceCells[patchFacei]] +=
                    (pMagSf[patchFacei]/magSqr(d))*sqr(d);
            }
        }
    }

    // Invert the dd tensor. A cell whose face stencil does not span three
    // directions has a singular dd, whose determinant can round to exactly
    // zero, and the division then raises a floating point exception on a
    // build that traps them (issue #417). The inverse is not needed there:
    // calcWideStencilVectors replaces the face stencil in every such cell and
    // switches these vectors off. The unit tensor of the mesh's directions is
    // inverted instead. It keeps the empty directions of a 2-D mesh empty,
    // which matters on OpenFOAM.org and foam-extend, where
    // inv(symmTensorField) reads the empty directions off the first tensor
    // and applies them to the whole field
    symmTensorField ddInvertible(dd);
    {
        const symmTensor emptyDirs(emptyDirections(mesh));
        const symmTensor unitSpanned(symmTensor(I) - emptyDirs);

        forAll(dd, cellI)
        {
            scalar lambdaMin = 0;
            scalar lambdaMax = 0;

            if (illConditioned(dd[cellI], emptyDirs, lambdaMin, lambdaMax))
            {
                ddInvertible[cellI] = unitSpanned;
            }
        }
    }
    const symmTensorField invDd(inv(ddInvertible));

    // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
    forAll(owner, facei)
    {
        if (interface[facei])
        {
            // Set face contribution to zero across interfaces
            pVectors_[facei] = vector::zero;
            nVectors_[facei] = vector::zero;

            continue;
        }

        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        scalar magSfByMagSqrd = magSf[facei]/magSqr(d);

        pVectors_[facei] = (1 - w[facei])*magSfByMagSqrd*(invDd[own] & d);
        nVectors_[facei] = -w[facei]*magSfByMagSqrd*(invDd[nei] & d);
    }

    forAll(pVectorsBf, patchi)
    {
        fvsPatchVectorField& patchLsP = pVectorsBf[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();
        const vectorField& Cf = p.Cf();

        // Build the d-vectors
        // In OF.com/OF.org, p.delta are the orthogonal components of the real d
        // vectors, so we need to build them ourselves
        //vectorField pd(p.delta());

        if (pw.coupled())
        {
            // Coupled d vectors
            const vectorField pd(p.delta());

            forAll(pd, patchFacei)
            {
                if (interfaceCoupled[patchi][patchFacei])
                {
                    // Set face contribution to zero across interfaces
                    patchLsP[patchFacei] = vector::zero;
                    continue;
                }

                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    ((1 - pw[patchFacei])*pMagSf[patchFacei]/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
            }
        }
        else if
        (
            isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchi])
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchi])
#endif
        )
        {
            // Treat symmetry planes consistently with internal faces
            // Use the mirrored face-cell values rather than the patch face
            // values
            // See https://doi.org/10.1080/10407790.2022.2105073
            const vectorField nHat(p.nf());
            forAll(nHat, patchFacei)
            {
                // Vector from the face-cell centre to the mirrored face-cell
                // centred
                const vector d =
                    transform
                    (
                        I - 2.0*sqr(nHat[patchFacei]),
                        C[faceCells[patchFacei]]
                    )
                  - C[faceCells[patchFacei]];

                // Multiply by 0.5 to be consistent with the internal field
                // where we multiply by the interpolation weight
                patchLsP[patchFacei] =
                    0.5*pMagSf[patchFacei]*(1.0/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
            }
        }
        else if (useBoundaryFaceValues_[patchi])
        {
            vectorField pd(p.size(), vector::zero);
            forAll(pd, faceI)
            {
                pd[faceI] = Cf[faceI] - C[faceCells[faceI]];
            }

            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    pMagSf[patchFacei]*(1.0/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
            }
        }
    }

    // Replace the face stencil where it is rank-deficient
    calcWideStencilVectors(dd);

    DebugInfo
        << "Finished calculating least square gradient vectors" << nl;
}


void Foam::leastSquaresS4fVectors::calcWideStencilVectors
(
    const symmTensorField& dd
) const
{
    const fvMesh& mesh = this->mesh();

    // Set local references to mesh data
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();

    // Flag cells whose face stencil does not span three directions, with the
    // empty directions of a 2-D or axisymmetric mesh filled first
    const symmTensor emptyDirs(emptyDirections(mesh));

    // A simplex cell carries the smallest face stencil a mesh can offer: only
    // nDim + 1 face neighbours, i.e. one more equation than unknowns. Such a
    // stencil is formally full rank, so the eigenvalue test above sees nothing
    // wrong with it, but it is far too small to reconstruct a cell-centred
    // gradient accurately. Widen it in every simplex cell, independently of the
    // eigenvalue test. Faces on empty patches carry no neighbour and so are
    // discounted: a tetrahedron has 4 faces, and a triangular prism in a
    // one-cell-thick 2-D mesh has 3 + 2 = 5
    const label nEmptyDirs = 3 - mesh.nGeometricD();
    const label maxSimplexFaces = mesh.nGeometricD() + 1 + 2*nEmptyDirs;

    const cellList& cells = mesh.cells();

    labelList cellToWide(mesh.nCells(), -1);
    label nWide = 0;
    label nSimplex = 0;
    scalar maxCond = 0.0;

    forAll(dd, cellI)
    {
        scalar lambdaMin = 0;
        scalar lambdaMax = 0;

        const bool deficient =
            illConditioned(dd[cellI], emptyDirs, lambdaMin, lambdaMax);

        const bool simplex = (cells[cellI].size() <= maxSimplexFaces);

        if (deficient || simplex)
        {
            cellToWide[cellI] = nWide++;

            if (simplex && !deficient)
            {
                nSimplex++;
            }
        }
        else
        {
            maxCond = max(maxCond, lambdaMax/lambdaMin);
        }
    }

    if (debug)
    {
        Info<< "    leastSquaresS4fVectors: "
            << returnReduce(nWide, sumOp<label>()) << " of "
            << returnReduce(mesh.nCells(), sumOp<label>())
            << " cells use the wide stencil, of which "
            << returnReduce(nSimplex, sumOp<label>())
            << " are simplex cells with a full-rank face stencil,"
            << " max cond(dd) = "
            << returnReduce(maxCond, maxOp<scalar>()) << endl;
    }

    wideCells_.setSize(nWide);
    wideStencil_.setSize(nWide);
    wideVectors_.setSize(nWide);

    // Decided globally, not per processor: what follows exchanges data across
    // coupled patches and reduces, and a processor that returned here would
    // leave the others waiting on it
    if (returnReduce(nWide, sumOp<label>()) == 0)
    {
        return;
    }

    forAll(cellToWide, cellI)
    {
        if (cellToWide[cellI] != -1)
        {
            wideCells_[cellToWide[cellI]] = cellI;
        }
    }

    // Collect the cells sharing at least one point with each flagged cell.
    // Unlike a face neighbour, a point neighbour is not removed by the
    // traction rule, as it is not reached across a boundary face
    const labelListList& cellPoints = mesh.cellPoints();
    const labelListList& pointCells = mesh.pointCells();

    // A point neighbour of a different material is no more admissible than a
    // face neighbour of one: either way the gradient would be built from two
    // materials at once. Every cell is -1 when the case has a single material,
    // so this costs nothing there
    const labelList materialID(cellMaterialID(mesh));

    // As for the face stencil, an interface running along a processor
    // boundary has to be skipped explicitly
    const List<boolList> interfaceCoupled(cellZoneInterfaceCoupled(mesh));

    forAll(wideCells_, wcI)
    {
        const label cellI = wideCells_[wcI];
        const labelList& curCellPoints = cellPoints[cellI];
        const label cellMaterial = materialID[cellI];

        labelHashSet stencil;
        forAll(curCellPoints, cpI)
        {
            const labelList& curPointCells = pointCells[curCellPoints[cpI]];

            forAll(curPointCells, pcI)
            {
                if (materialID[curPointCells[pcI]] == cellMaterial)
                {
                    stencil.insert(curPointCells[pcI]);
                }
            }
        }
        stencil.erase(cellI);

        wideStencil_[wcI] = stencil.toc();
    }

    // Assemble the wide moment tensor. All contributions use an inverse
    // distance squared weight: the face-area weight used by the face stencil
    // has no meaning for a point neighbour, which shares no face
    symmTensorField wideDd(nWide, symmTensor::zero);

    forAll(wideCells_, wcI)
    {
        const label cellI = wideCells_[wcI];
        const labelList& stencil = wideStencil_[wcI];

        forAll(stencil, i)
        {
            const vector d = C[stencil[i]] - C[cellI];

            wideDd[wcI] += (1.0/magSqr(d))*sqr(d);
        }
    }

    // Boundary faces of the flagged cells contribute as they do for the face
    // stencil: coupled and known-value patches give the face, symmetry planes
    // give the mirrored cell, and traction patches give nothing
    forAll(mesh.boundary(), patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();
        const vectorField& Cf = p.Cf();

        if (pw.coupled())
        {
            const vectorField pd(p.delta());

            forAll(pd, patchFacei)
            {
                const label wcI = cellToWide[faceCells[patchFacei]];

                if (wcI != -1 && !interfaceCoupled[patchi][patchFacei])
                {
                    const vector& d = pd[patchFacei];

                    wideDd[wcI] += (1.0/magSqr(d))*sqr(d);
                }
            }
        }
        else if
        (
            isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchi])
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchi])
#endif
        )
        {
            const vectorField nHat(p.nf());

            forAll(nHat, patchFacei)
            {
                const label cellI = faceCells[patchFacei];
                const label wcI = cellToWide[cellI];

                if (wcI != -1)
                {
                    const vector d =
                        transform(I - 2.0*sqr(nHat[patchFacei]), C[cellI])
                      - C[cellI];

                    wideDd[wcI] += (1.0/magSqr(d))*sqr(d);
                }
            }
        }
        else if (useBoundaryFaceValues_[patchi])
        {
            forAll(faceCells, patchFacei)
            {
                const label cellI = faceCells[patchFacei];
                const label wcI = cellToWide[cellI];

                if (wcI != -1)
                {
                    const vector d = Cf[patchFacei] - C[cellI];

                    wideDd[wcI] += (1.0/magSqr(d))*sqr(d);
                }
            }
        }
    }

    // Check whether any cell is still rank-deficient after widening
    label nStillSingular = 0;
    label worstCell = -1;
    scalar worstRatio = GREAT;
    boolList stillSingular(nWide, false);

    forAll(wideDd, wcI)
    {
        scalar lambdaMin = 0;
        scalar lambdaMax = 0;

        if (illConditioned(wideDd[wcI], emptyDirs, lambdaMin, lambdaMax))
        {
            nStillSingular++;
            stillSingular[wcI] = true;

            const scalar ratio =
                lambdaMax < VSMALL ? 0.0 : lambdaMin/lambdaMax;

            if (ratio < worstRatio)
            {
                worstRatio = ratio;
                worstCell = wideCells_[wcI];
            }
        }
    }

    reduce(nStillSingular, sumOp<label>());

    // Every cell is -1 with a single material, on every processor
    bool multiMaterial = false;
    forAll(materialID, cellI)
    {
        if (materialID[cellI] != -1)
        {
            multiMaterial = true;
            break;
        }
    }
    reduce(multiMaterial, orOp<bool>());

    if (nStillSingular > 0 && !multiMaterial)
    {
        // With a single material nothing has been filtered, so the stencil is
        // the one this scheme used before the filtering arrived. A body one
        // cell thick in some direction produces this in serial - the striker
        // in pipeCrush is a single row of cells - and a point-cell stencil
        // truncated at a processor boundary can produce it in parallel.
        // Neither is new, so neither is fatal: the gradient is reconstructed
        // in the directions the stencil spans and is zero along the others
        WarningInFunction
            << nStillSingular
            << " cells remain rank-deficient after widening the gradient"
            << " stencil to point neighbours." << nl
            << "    Worst on this processor: cell " << worstCell
            << ", smallest eigenvalue ratio " << worstRatio << nl
            << "    Either the body is one cell thick in some direction or, in"
            << " parallel, the point-cell stencil is truncated at a processor"
            << " boundary." << endl;
    }
    else if (nStillSingular > 0)
    {
        // A cell whose stencil cannot span the mesh's directions has no
        // gradient that can be reconstructed from it, and with several
        // materials there is no sound fallback: widening again would not help,
        // and falling back on an unfiltered stencil would build the gradient
        // from two materials at once, which is what the filtering exists to
        // prevent.
        //
        // Two geometries produce this. The cell's own material does not
        // surround it - a material one cell thick, or a cell at a material
        // corner or tip - or the point-cell stencil has been truncated at a
        // processor boundary. Both mean the same thing about the answer, so
        // both are fatal
        FatalErrorInFunction
            << nStillSingular
            << " cells remain rank-deficient after widening the gradient"
            << " stencil to point neighbours of the same material." << nl
            << "    Worst on this processor: cell " << worstCell
            << ", material " << (worstCell >= 0 ? materialID[worstCell] : -1)
            << ", smallest eigenvalue ratio " << worstRatio << nl
            << "    Either the cell's own material does not surround it in"
            << " every direction, in which case refine so that at least two"
            << " cells span the material everywhere; or the stencil is"
            << " truncated at a processor boundary, in which case decompose"
            << " so that it is not."
            << exit(FatalError);
    }

    // Invert the wide tensors. A flagged cell is kept out of the field
    // inversion, whose division by a determinant that rounds to exactly zero
    // raises a floating point exception on a build that traps them (issue
    // #417), and whose treatment of the empty directions of a 2-D mesh is
    // read off the first tensor on OpenFOAM.org and foam-extend: it is given
    // the unit tensor of the mesh's directions there, and its own inverse
    // below. That is the exact inverse when the tensor is invertible, so a
    // linear field is still reproduced in every direction the stencil spans,
    // and the pseudo-inverse when it is not. A cell with no admissible
    // stencil at all has nothing to reconstruct a gradient from
    symmTensorField wideDdInvertible(wideDd);
    {
        const symmTensor unitSpanned(symmTensor(I) - emptyDirs);

        forAll(wideDd, wcI)
        {
            if (stillSingular[wcI])
            {
                wideDdInvertible[wcI] = unitSpanned;
            }
        }
    }
    symmTensorField wideInvDd(inv(wideDdInvertible));

    forAll(wideDd, wcI)
    {
        if (stillSingular[wcI])
        {
            label rank = 0;
            wideInvDd[wcI] = pseudoInverse(wideDd[wcI], rank);

            if (rank == 0)
            {
                FatalErrorInFunction
                    << "Cell " << wideCells_[wcI] << " has no admissible "
                    << "gradient stencil: no face or point neighbour of its "
                    << "own material, and no known-value boundary face."
                    << exit(FatalError);
            }
        }
    }

    forAll(wideCells_, wcI)
    {
        const label cellI = wideCells_[wcI];
        const labelList& stencil = wideStencil_[wcI];
        vectorList& lsVecs = wideVectors_[wcI];
        lsVecs.setSize(stencil.size(), vector::zero);

        forAll(stencil, i)
        {
            const vector d = C[stencil[i]] - C[cellI];

            lsVecs[i] = (1.0/magSqr(d))*(wideInvDd[wcI] & d);
        }
    }

    // Switch off the face-stencil contribution in the flagged cells. The
    // vectors are stored per side, so silencing one cell leaves the cell across
    // the face untouched
    forAll(owner, facei)
    {
        if (cellToWide[owner[facei]] != -1)
        {
            pVectors_[facei] = vector::zero;
        }

        if (cellToWide[neighbour[facei]] != -1)
        {
            nVectors_[facei] = vector::zero;
        }
    }

    // Rebuild the boundary vectors of the flagged cells with the wide tensor
    auto& pVectorsBf = boundaryFieldRef(pVectors_);

    forAll(pVectorsBf, patchi)
    {
        fvsPatchVectorField& patchLsP = pVectorsBf[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();
        const vectorField& Cf = p.Cf();

        if (pw.coupled())
        {
            const vectorField pd(p.delta());

            forAll(pd, patchFacei)
            {
                const label wcI = cellToWide[faceCells[patchFacei]];

                if (wcI != -1 && !interfaceCoupled[patchi][patchFacei])
                {
                    const vector& d = pd[patchFacei];

                    patchLsP[patchFacei] =
                        (1.0/magSqr(d))*(wideInvDd[wcI] & d);
                }
            }
        }
        else if
        (
            isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchi])
#ifdef OPENFOAM_NOT_EXTEND
         || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchi])
#endif
        )
        {
            const vectorField nHat(p.nf());

            forAll(nHat, patchFacei)
            {
                const label cellI = faceCells[patchFacei];
                const label wcI = cellToWide[cellI];

                if (wcI != -1)
                {
                    const vector d =
                        transform(I - 2.0*sqr(nHat[patchFacei]), C[cellI])
                      - C[cellI];

                    patchLsP[patchFacei] =
                        (1.0/magSqr(d))*(wideInvDd[wcI] & d);
                }
            }
        }
        else if (useBoundaryFaceValues_[patchi])
        {
            forAll(faceCells, patchFacei)
            {
                const label cellI = faceCells[patchFacei];
                const label wcI = cellToWide[cellI];

                if (wcI != -1)
                {
                    const vector d = Cf[patchFacei] - C[cellI];

                    patchLsP[patchFacei] =
                        (1.0/magSqr(d))*(wideInvDd[wcI] & d);
                }
            }
        }
    }

    DebugInfo
        << "    wide point-cell stencil applied" << nl;
}


#ifdef OPENFOAM_NOT_EXTEND

bool Foam::leastSquaresS4fVectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}

#else

bool Foam::leastSquaresS4fVectors::movePoints() const
{
    calcLeastSquaresVectors();
    return true;
}

bool Foam::leastSquaresS4fVectors::updateMesh(const mapPolyMesh&) const
{
    calcLeastSquaresVectors();
    return true;
}

#endif // ifdef OPENFOAM_NOT_EXTEND

// ************************************************************************* //
