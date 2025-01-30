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
#include "symmetryPlanePolyPatch.H"

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
    MeshObject<fvMesh, Foam::MoveableMeshObject, leastSquaresS4fVectors>(mesh),
    useBoundaryFaceValues_(useBoundaryFaceValues_),
    pVectors_
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, Zero)
    ),
    nVectors_
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, Zero)
    )
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresS4fVectors::~leastSquaresS4fVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresS4fVectors::calcLeastSquaresVectors()
{
    DebugInFunction
        << "Calculating least square gradient vectors" << nl;

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceScalarField& magSf = mesh.magSf();


    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), Zero);
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        symmTensor wdd = (magSf[facei]/magSqr(d))*sqr(d);

        dd[own] += (1 - w[facei])*wdd;
        dd[nei] += w[facei]*wdd;
    }


    surfaceVectorField::Boundary& pVectorsBf =
        pVectors_.boundaryFieldRef();

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
         || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchi])
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

                dd[faceCells[patchFacei]] +=
                    (pMagSf[patchFacei]/magSqr(d))*sqr(d);
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
         || isA<symmetryPlanePolyPatch>(mesh.boundaryMesh()[patchi])
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

                patchLsP[patchFacei] =
                    pMagSf[patchFacei]*(1.0/magSqr(d))
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

    DebugInfo
        << "Finished calculating least square gradient vectors" << nl;
}


bool Foam::leastSquaresS4fVectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}


// ************************************************************************* //
