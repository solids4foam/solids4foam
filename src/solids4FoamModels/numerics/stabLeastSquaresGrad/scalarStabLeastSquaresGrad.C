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

Description
    Specialisation of leastSquaresGrad for scalars. Needed for implicit fvmGrad
    operator for block coupled systems.

\*---------------------------------------------------------------------------*/

#include "stabLeastSquaresGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, vector> > stabLeastSquaresGrad<scalar>::fvmGrad
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();

    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    volScalarField cellV
    (
        IOobject
        (
            "cellV",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    );
    cellV.internalField() = mesh.V();
    cellV.correctBoundaryConditions();
    const scalarField& cellVIn = cellV.internalField();

    const surfaceScalarField& w = mesh.weights();

    tmp<BlockLduSystem<vector, vector> > tbs
    (
        new BlockLduSystem<vector, vector>(mesh)
    );
    BlockLduSystem<vector, vector>& bs = tbs();
    vectorField& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
    CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
    CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

    // Get references to least square vectors
    const stabLeastSquaresVectors& lsv = stabLeastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    // Internal field
    const vectorField& ownLsIn = ownLs.internalField();
    const vectorField& neiLsIn = neiLs.internalField();

    label owner, neighbour;

    forAll (nei, faceI)
    {
        owner = own[faceI];
        neighbour = nei[faceI];

        u[faceI] = cellVIn[owner]*ownLsIn[faceI];
        l[faceI] = cellVIn[neighbour]*neiLsIn[faceI];

        // Caution - this is NOT negSumDiag(). VV, 17/July/2014.
        d[owner] -= u[faceI];
        d[neighbour] -= l[faceI];
    }

    // Boundary contributions
    forAll (vf.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& pownLs = ownLs.boundaryField()[patchI];
        const fvsPatchScalarField& pw = w.boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        // Part of diagonal contribution irrespective of the patch type
        forAll (pf, faceI)
        {
            const label cellI = fc[faceI];
            d[cellI] -= cellVIn[cellI]*pownLs[faceI];
        }

        if (patch.coupled())
        {
            const vectorField& pneiLs = neiLs.boundaryField()[patchI];
            const scalarField cellVInNei =
                cellV.boundaryField()[patchI].patchNeighbourField();

            CoeffField<vector>::linearTypeField& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();
            CoeffField<vector>::linearTypeField& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            // Coupling  and diagonal contributions
            forAll (pf, faceI)
            {
                pcoupleUpper[faceI] -= cellVIn[fc[faceI]]*pownLs[faceI];
                pcoupleLower[faceI] -= cellVInNei[faceI]*pneiLs[faceI];
            }
        }
        else
        {
            const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            // Diagonal and source contributions depending on the patch type
            forAll (pf, faceI)
            {
                const label cellI = fc[faceI];

                d[cellI] += cellVIn[cellI]*pownLs[faceI]*internalCoeffs[faceI];

                source[cellI] -= cellVIn[cellI]*pownLs[faceI]*
                    boundaryCoeffs[faceI];
            }
        }
    }

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
