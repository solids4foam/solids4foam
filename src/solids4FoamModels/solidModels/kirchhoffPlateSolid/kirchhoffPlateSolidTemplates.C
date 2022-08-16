/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "kirchhoffPlateSolid.H"
#include "GeometricFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::solidModels::kirchhoffPlateSolid::
mapAreaFieldToSingleLayerVolumeField
(
    const GeometricField<Type, faPatchField, areaMesh>& af,
    GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Take a references
    const polyMesh& pMesh = mesh();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();
    const labelList& faceLabels = aMesh_.faceLabels();
    const label areaPatchID = areaPatch().index();
    const label pMeshNFaces = pMesh.nFaces();
    const Field<Type> afI = af.internalField();
    FieldField<fvPatchField, Type>& bf = vf.boundaryField();


    // 1. Map af to vf areaPatch

    forAll(faceLabels, aFaceI)
    {
        const label faceID = faceLabels[aFaceI];

        // Escape if face is beyond active faces, eg belongs to a face zone
        if (faceID < pMeshNFaces)
        {
            const label localFaceID = bm[areaPatchID].whichFace(faceID);
            bf[areaPatchID][localFaceID] = afI[aFaceI];
        }
    }


    // 2. Next, we map to the vf internal field

    {
        const unallocLabelList& faceCells = areaPatch().faceCells();
        const Field<Type>& patchField = vf.boundaryField()[areaPatchID];

        forAll(faceCells, faceI)
        {
            vf[faceCells[faceI]] = patchField[faceI];
        }
    }


    // 3. Finally, we map to the areaShadowPatch

    {
        const unallocLabelList& faceCells = areaShadowPatch().faceCells();
        Field<Type>& patchField = vf.boundaryField()[areaShadowPatch().index()];

        forAll(faceCells, faceI)
        {
            patchField[faceI] = vf[faceCells[faceI]];
        }
    }
}

// ************************************************************************* //
