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

#include "stabLeastSquaresGrad.H"
#include "stabLeastSquaresVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
stabLeastSquaresGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vsf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vsf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad();

    // Get reference to least square vectors
    const stabLeastSquaresVectors& lsv = stabLeastSquaresVectors::New(mesh);

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    // Get access to internal field

    const Field<Type>& vsfIn = vsf.internalField();

    Field<GradType>& lsGradIn = lsGrad.internalField();

    const vectorField& ownLsIn = ownLs.internalField();
    const vectorField& neiLsIn = neiLs.internalField();

    label ownFaceI, neiFaceI;

    forAll (own, facei)
    {
        ownFaceI = own[facei];
        neiFaceI = nei[facei];

        Type deltaVsf = vsfIn[neiFaceI] - vsfIn[ownFaceI];

        lsGradIn[ownFaceI] += ownLsIn[facei]*deltaVsf;
        lsGradIn[neiFaceI] -= neiLsIn[facei]*deltaVsf;
    }

    // Boundary faces
    forAll (vsf.boundaryField(), patchi)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];

        const unallocLabelList& faceCells =
            lsGrad.boundaryField()[patchi].patch().faceCells();

        if (vsf.boundaryField()[patchi].coupled())
        {
            Field<Type> neiVsf =
                vsf.boundaryField()[patchi].patchNeighbourField();

            forAll (neiVsf, patchFaceI)
            {
                lsGrad[faceCells[patchFaceI]] +=
                    patchOwnLs[patchFaceI]
                   *(neiVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
            }
        }
        else
        {
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];

            forAll (patchVsf, patchFaceI)
            {
                lsGrad[faceCells[patchFaceI]] +=
                     patchOwnLs[patchFaceI]
                    *(patchVsf[patchFaceI] - vsf[faceCells[patchFaceI]]);
            }
        }
    }

    lsGrad.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}

template<class Type>
tmp
<
    BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
> stabLeastSquaresGrad<Type>::fvmGrad
(
   const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
   FatalErrorIn
   (
       "tmp<BlockLduSystem> fvmGrad\n"
       "(\n"
       "    GeometricField<Type, fvPatchField, volMesh>&"
       ")\n"
   )   << "Implicit gradient operator defined only for scalar."
       << abort(FatalError);

   typedef typename outerProduct<vector, Type>::type GradType;

   tmp<BlockLduSystem<vector, GradType> > tbs
   (
       new BlockLduSystem<vector, GradType>(vf.mesh())
   );

   return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
