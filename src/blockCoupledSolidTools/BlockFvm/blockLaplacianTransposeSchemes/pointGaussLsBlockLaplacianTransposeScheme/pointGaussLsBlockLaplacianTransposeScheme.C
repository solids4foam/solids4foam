/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pointGaussLsBlockLaplacianTransposeScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

#include "BlockFvmDivSigma.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "solidPolyMesh.H"
#include "blockLduSolvers.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "emptyFvPatchFields.H"
#include "blockFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(fv::pointGaussLsBlockLaplacianTransposeScheme, 0);
    fv::blockLaplacianTranspose::addIstreamConstructorToTable
    <
        fv::pointGaussLsBlockLaplacianTransposeScheme
    > addPointGaussLsBlockLaplacianTransposeSchemeIstreamConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void pointGaussLsBlockLaplacianTransposeScheme::insertCoeffsNorm
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB,
    BlockLduMatrix<vector>& blockM
)
{
    // Const reference to fvMesh
    const fvMesh& mesh = U.mesh();

    // Grab block diagonal
    Field<tensor>& d = blockM.diag().asSquare();

    // Grab linear off-diagonal
    Field<tensor>& l = blockM.lower().asSquare();
    Field<tensor>& u = blockM.upper().asSquare();

    // Normal derivative terms

    // Note that the non-orthogonal correction for these normal derivative terms
    // is added in the insertCoeffTang

    if (fv::debug)
    {
        Info<< nl << "Adding normal derivatve terms to matrix" << endl;
    }

    const unallocLabelList& fvOwner = mesh.owner();
    const unallocLabelList& fvNeighbour = mesh.neighbour();
    const surfaceVectorField n = mesh.Sf()/mesh.magSf();
    const vectorField nI = n.internalField();
    const scalarField& mufI = muf.internalField();
    const scalarField& magSfI = mesh.magSf().internalField();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs().internalField();
    const labelList& fvMap = solidMesh.fvMeshAddressingMap();

    // Internal faces

    forAll(fvOwner, faceI)
    {
        const label own = fvOwner[faceI];
        const label nei = fvNeighbour[faceI];

        const vector& faceN = nI[faceI];
        const scalar faceMu = mufI[faceI];
        const scalar faceMagSf = magSfI[faceI];
        const scalar faceDeltaCoeff = deltaCoeffsI[faceI];

        // Normal derivative terms
        const tensor coeff = sqr(faceN)*faceMu*faceMagSf*faceDeltaCoeff;

        d[own] -= coeff;
        d[nei] -= coeff;

        const label varI = fvMap[faceI];
        u[varI] += coeff;
        l[varI] += coeff;

        if (fv::debug > 1)
        {
            Info<< nl << "face " << faceI << nl
                << "    d[" << own << "] " << -coeff << nl
                << "    d[" << nei << "] " << -coeff << nl
                << "    u[" << fvMap[faceI] << "] " << coeff << nl
                << "    l[" << fvMap[faceI] << "] " << coeff << endl;
        }
    }

    // Boundary faces

    forAll(mesh.boundaryMesh(), patchI)
    {
        const word& patchType = mesh.boundaryMesh()[patchI].type();

        if (patchType == emptyPolyPatch::typeName)
        {
            // Do nothing
        }
        else if (patchType == processorPolyPatch::typeName)
        {
            const vectorField& pN = n.boundaryField()[patchI];
            const scalarField& pMu = muf.boundaryField()[patchI];
            const scalarField& pDeltaCoeffs =
                mesh.deltaCoeffs().boundaryField()[patchI];
            const scalarField& pMagSf =
                mesh.magSf().boundaryField()[patchI];
            const unallocLabelList& faceCells =
                mesh.boundary()[patchI].faceCells();

            Field<tensor>& coupleUpper =
                blockM.coupleUpper()[patchI].asSquare();

            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
                const label own = faceCells[faceI];
                const vector& faceN = pN[faceI];
                const scalar faceMu = pMu[faceI];
                const scalar faceMagSf = pMagSf[faceI];
                const scalar faceDeltaCoeff = pDeltaCoeffs[faceI];

                // Normal derivative terms
                const tensor coeff = sqr(faceN)*faceMu*faceMagSf*faceDeltaCoeff;

                // Add diagonal contribution to owner cell
                d[own] -= coeff;

                // Add upper contribution for owner cell from shadow cell across
                // the proc boundary
                coupleUpper[faceI] += coeff;
            }
        }
        else
        {
            const vectorField& pN = n.boundaryField()[patchI];
            const scalarField& pMu = muf.boundaryField()[patchI];
            const scalarField& pDeltaCoeffs =
                mesh.deltaCoeffs().boundaryField()[patchI];
            const scalarField& pMagSf =
                mesh.magSf().boundaryField()[patchI];
            const unallocLabelList& faceCells =
                mesh.boundary()[patchI].faceCells();

            const label start = mesh.boundaryMesh()[patchI].start();

            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
                const label own = faceCells[faceI];
                const vector& faceN = pN[faceI];
                const scalar faceMu = pMu[faceI];
                const scalar faceMagSf = pMagSf[faceI];
                const scalar faceDeltaCoeff = pDeltaCoeffs[faceI];

                const label curFaceID = start + faceI;

                // Normal derivative terms
                const tensor coeff = sqr(faceN)*faceMu*faceMagSf*faceDeltaCoeff;


                // Add diagonal contribution to owner cell
                d[own] -= coeff;

                const label curImEdgeID =
                    solidMesh.findCellFaceImplicitBond(own, curFaceID);

                if (curImEdgeID == -1)
                {
                    FatalErrorIn
                    (
                        "pointGaussLsBlockLaplacianTransposeScheme::"
                        "insertCoeffNorm()"
                    )   << "curImEdge not found"
                        << abort(FatalError);
                }

                // Add upper contribution for owner cell
                u[curImEdgeID] += coeff;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //


pointGaussLsBlockLaplacianTransposeScheme::
pointGaussLsBlockLaplacianTransposeScheme
(
    const fvMesh& mesh,
    Istream& is
)
:
    blockLaplacianTranspose(mesh, is)
{}


// * * * * * * * * * Public member functions  * * * * * * * * * * * * * * * * //

tmp<BlockLduMatrix<vector> >
pointGaussLsBlockLaplacianTransposeScheme::fvmBlockLaplacianTranspose
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB
)
{
    tmp<BlockLduMatrix<vector> > tBlockM
    (
        new BlockLduMatrix<vector>(solidMesh)
    );
    BlockLduMatrix<vector>& blockM = tBlockM();

    // Grab block diagonal and set it to zero
    Field<tensor>& d = blockM.diag().asSquare();
    d = tensor::zero;

    // Grab linear off-diagonal and set it to zero
    Field<tensor>& l = blockM.lower().asSquare();
    Field<tensor>& u = blockM.upper().asSquare();
    u = tensor::zero;
    l = tensor::zero;

    // Check size of source
    if (blockB.size() != solidMesh.nVariables())
    {
        FatalErrorIn("pointGaussLsBlockLaplacianTransposeScheme::fvmDivSigma")
            << "The source vector has the wrong length" << abort(FatalError);
    }

    // Initialise LDU interfaces
    blockM.interfaces() = U.boundaryField().blockInterfaces();

    // Reset all global coefficients to zero, they will be non-zero from
    // previous time-steps: this is done is the solver at the moment but we
    // should have a better way
    //solidMesh.clearOutGlobalCoeffs();

    // Insert coeffs due to normal derivative terms
    insertCoeffsNorm(solidMesh, muf, U, blockB, blockM);

    // Insert coeffs due to tangential derivative terms
    insertCoeffsTang(solidMesh, muf, U, blockB, blockM);

    // Insert coeffs due to the boundary conditions
    //insertCoeffsBc(solidMesh, muf, U, blockB, blockM);

    return tBlockM;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
