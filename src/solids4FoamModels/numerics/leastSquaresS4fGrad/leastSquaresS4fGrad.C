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

#include "leastSquaresS4fGrad.H"
#include "leastSquaresS4fVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "solidTractionFvPatchVectorField.H"
#include "boolIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::leastSquaresS4fGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vsf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tlsGrad
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
            dimensioned<GradType>(vsf.dimensions()/dimLength, Zero),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad.ref();

    // Lookup the useBoundaryFaceValues list
    // This list needs to be created before calling the gradient operation
    // The list defines which patch values should be used in the least squares
    // calculation.
    // There must be a nice way to define this on the fly, but for now this
    // solution is OK, as it avoids misuse
    const word useBndFcValName("useBoundaryFaceValues_" + vsf.name());
    if (!mesh.foundObject<boolIOList>(useBndFcValName))
    {
        FatalErrorInFunction
            << useBndFcValName << " boolIOList not found! " << nl
            << "To use the leastSquaresS4fGrad scheme, you must first define an"
            << " boolIOList called " << useBndFcValName << ", which indicates"
            << " which patches should be used in the least squares problems for"
            << " the " << vsf.name() << " field."
            << abort(FatalError);
    }
    const boolIOList& useBoundaryFaceValues =
        mesh.lookupObject<boolIOList>("useBoundaryFaceValues_" + vsf.name());

    // Get reference to least square vectors
    const leastSquaresS4fVectors& lsv =
        leastSquaresS4fVectors::New
        (
            "leastSquaresVectors" + vsf.name(), mesh, useBoundaryFaceValues
        );

    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    forAll(own, facei)
    {
        label ownFacei = own[facei];
        label neiFacei = nei[facei];

        Type deltaVsf = vsf[neiFacei] - vsf[ownFacei];

        lsGrad[ownFacei] += ownLs[facei]*deltaVsf;
        lsGrad[neiFacei] -= neiLs[facei]*deltaVsf;
    }

    // Boundary faces
    forAll(vsf.boundaryField(), patchi)
    {
        const fvsPatchVectorField& patchOwnLs = ownLs.boundaryField()[patchi];

        const labelUList& faceCells =
            vsf.boundaryField()[patchi].patch().faceCells();

        if (vsf.boundaryField()[patchi].coupled())
        {
            const Field<Type> neiVsf
            (
                vsf.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(neiVsf, patchFacei)
            {
                lsGrad[faceCells[patchFacei]] +=
                    patchOwnLs[patchFacei]
                   *(neiVsf[patchFacei] - vsf[faceCells[patchFacei]]);
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
            const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];
            const vectorField nHat(patchVsf.patch().nf());
            forAll(patchVsf, patchFacei)
            {
                lsGrad[faceCells[patchFacei]] +=
                     patchOwnLs[patchFacei]
                    *(
                        transform
                        (
                            I - 2.0*sqr(nHat[patchFacei]),
                            vsf[faceCells[patchFacei]]
                        )
                      - vsf[faceCells[patchFacei]]
                    );
            }
        }
        else
        {
            if (useBoundaryFaceValues[patchi])
            {
                const fvPatchField<Type>& patchVsf =
                    vsf.boundaryField()[patchi];
                forAll(patchVsf, patchFacei)
                {
                    lsGrad[faceCells[patchFacei]] +=
                         patchOwnLs[patchFacei]
                        *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);
                }
            }
        }
    }


    lsGrad.correctBoundaryConditions();

    // This causes convergence problems on non-orthogonal grids for the SNES
    // solver. This may be related to Rhie-Chow stabilisation
    // if (useNeumannBoundaryFaceValues_)
    // {
    //     // Replace the normal gradient on boundary faces
    //     //gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);
    // }

    return tlsGrad;
}


// ************************************************************************* //
