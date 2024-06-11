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

#ifdef OPENFOAM_COM

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
#ifdef OPENFOAM_NOT_EXTEND
    FieldField<fvPatchField, Type>& bf = vf.boundaryFieldRef();
#else
    FieldField<fvPatchField, Type>& bf = vf.boundaryField();
#endif

    // Check that bf area patch is not of type empty
    if (bf[areaPatchID].type() == "empty")
    {
        FatalErrorIn("mapAreaFieldToSingleLayerVolumeField(...)")
            << "Patch " << bf[areaPatchID].patch().name()
            << " must not be of type 'empty'"
            << abort(FatalError);
    }

    // 1. Map af to vf areaPatch
    Field<Type>& bfAreaPatch = bf[areaPatchID];
    forAll(faceLabels, aFaceI)
    {
        const label faceID = faceLabels[aFaceI];

        // Escape if face is beyond active faces, eg belongs to a face zone
        if (faceID < pMeshNFaces)
        {
            const label localFaceID = bm[areaPatchID].whichFace(faceID);
            bfAreaPatch[localFaceID] = afI[aFaceI];
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
#ifdef OPENFOAM_NOT_EXTEND
        Field<Type>& patchField =
            vf.boundaryFieldRef()[areaShadowPatch().index()];
#else
        Field<Type>& patchField =
            vf.boundaryField()[areaShadowPatch().index()];
#endif

        forAll(faceCells, faceI)
        {
            patchField[faceI] = vf[faceCells[faceI]];
        }
    }
}

#endif // OPENFOAM_COM

// ************************************************************************* //
