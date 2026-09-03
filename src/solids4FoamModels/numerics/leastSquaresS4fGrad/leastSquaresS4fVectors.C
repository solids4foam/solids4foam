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

    // Invert the dd tensor
    const symmTensorField invDd(inv(dd));

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

    // Flag cells whose face stencil does not span three directions.
    // dd is symmetric positive semi-definite, so its eigenvalues are real and
    // non-negative, and the smallest one vanishes exactly when the stencil is
    // confined to a plane or a line. Healthy and degenerate cells are separated
    // by many orders of magnitude, so this ratio is not sensitive: any value
    // between 1e-4 and 1e-12 selects the same cells
    const scalar minEigenRatio = 1e-4;

    // A 2-D or axisymmetric mesh carries no information in its empty
    // direction, and that is not a degeneracy. Fill the empty directions
    // before testing, as inv(symmTensorField) does through safeInv()
    const Vector<label> geometricD = mesh.geometricD();
    symmTensor emptyDirs(symmTensor::zero);
    emptyDirs.xx() = (geometricD.x() == -1 ? 1.0 : 0.0);
    emptyDirs.yy() = (geometricD.y() == -1 ? 1.0 : 0.0);
    emptyDirs.zz() = (geometricD.z() == -1 ? 1.0 : 0.0);

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
        const symmTensor ddc(dd[cellI] + (tr(dd[cellI])/3.0)*emptyDirs);
        const vector eVals = eigenValues(ddc);
        const scalar lambdaMin = eVals[vector::X];
        const scalar lambdaMax = eVals[vector::Z];

        const bool rankDeficient =
        (
            lambdaMax < VSMALL || lambdaMin <= minEigenRatio*lambdaMax
        );

        const bool simplex = (cells[cellI].size() <= maxSimplexFaces);

        if (rankDeficient || simplex)
        {
            cellToWide[cellI] = nWide++;

            if (simplex && !rankDeficient)
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

    if (nWide == 0)
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

    forAll(wideCells_, wcI)
    {
        const label cellI = wideCells_[wcI];
        const labelList& curCellPoints = cellPoints[cellI];

        labelHashSet stencil;
        forAll(curCellPoints, cpI)
        {
            const labelList& curPointCells = pointCells[curCellPoints[cpI]];

            forAll(curPointCells, pcI)
            {
                stencil.insert(curPointCells[pcI]);
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

                if (wcI != -1)
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

    // Warn if any cell is still rank-deficient after widening: in parallel the
    // point-cell stencil is truncated at processor boundaries
    label nStillSingular = 0;
    forAll(wideDd, wcI)
    {
        const symmTensor ddc(wideDd[wcI] + (tr(wideDd[wcI])/3.0)*emptyDirs);
        const vector eVals = eigenValues(ddc);

        if
        (
            eVals[vector::Z] < VSMALL
         || eVals[vector::X] <= minEigenRatio*eVals[vector::Z]
        )
        {
            nStillSingular++;
        }
    }

    if (returnReduce(nStillSingular, sumOp<label>()) > 0)
    {
        WarningInFunction
            << returnReduce(nStillSingular, sumOp<label>())
            << " cells remain rank-deficient after widening the gradient"
            << " stencil to point neighbours." << nl
            << "    In parallel this indicates the point-cell stencil is"
            << " truncated at a processor boundary." << endl;
    }

    const symmTensorField wideInvDd(inv(wideDd));

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

                if (wcI != -1)
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
