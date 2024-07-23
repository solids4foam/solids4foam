/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "leastSquaresS4fGrad.H"
#include "leastSquaresS4fVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "extrapolatedCalculatedFvPatchField.H"

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

    // Get reference to least square vectors
    const leastSquaresS4fVectors& lsv = leastSquaresS4fVectors::New(mesh);

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
        // PC: do not include boundary faces
        // else
        // {
        //     const fvPatchField<Type>& patchVsf = vsf.boundaryField()[patchi];

        //     forAll(patchVsf, patchFacei)
        //     {
        //         lsGrad[faceCells[patchFacei]] +=
        //              patchOwnLs[patchFacei]
        //             *(patchVsf[patchFacei] - vsf[faceCells[patchFacei]]);
        //     }
        // }
    }


    lsGrad.correctBoundaryConditions();
    // PC: do not replace the normal gradient on boundary faces
    // gaussGrad<Type>::correctBoundaryConditions(vsf, lsGrad);

    return tlsGrad;
}


// ************************************************************************* //
