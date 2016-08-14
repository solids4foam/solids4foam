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

#include "mechanicalModel.H"
#include "processorFvsPatchField.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Type>
void Foam::mechanicalModel::mapSubMeshVolField
(
    const label meshI,
    const GeometricField<Type, fvPatchField, volMesh>& subMeshField,
    GeometricField<Type, fvPatchField, volMesh>& baseMeshField
) const
{
    // Map iternal field from sub-mesh to global mesh

    const labelList& cellMap = subMeshes()[meshI].cellMap();

    forAll(subMeshField, cellI)
    {
        baseMeshField[cellMap[cellI]] = subMeshField[cellI];
    }

    // Map boundary field

    const labelList& patchMap = subMeshes()[meshI].patchMap();
    const labelList& faceMap = subMeshes()[meshI].faceMap();

    forAll(subMeshField.boundaryField(), patchI)
    {
        const fvPatchField<Type>& subMeshFieldP =
            subMeshField.boundaryField()[patchI];

        const label start = subMeshFieldP.patch().patch().start();

        if (patchMap[patchI] != -1)
        {
            fvPatchField<Type>& baseMeshFieldP =
                baseMeshField.boundaryField()[patchMap[patchI]];

            if (!baseMeshFieldP.coupled())
            {
                forAll(subMeshFieldP, faceI)
                {
                    const label globalGlobalMeshFace = faceMap[start + faceI];

                    const label curGlobalMeshPatchFace =
                        globalGlobalMeshFace
                      - mesh().boundaryMesh()[patchMap[patchI]].start();

                    baseMeshFieldP[curGlobalMeshPatchFace] =
                        subMeshFieldP[faceI];
                }
            }
        }
        else // interface faces shared by two materials
        {
            forAll(subMeshFieldP, faceI)
            {
                const label globalGlobalMeshFace = faceMap[start + faceI];

                if (globalGlobalMeshFace < mesh().nInternalFaces())
                {
                    // Face value will be the average from both sides
                    baseMeshField[globalGlobalMeshFace] +=
                        0.5*subMeshFieldP[faceI];
                }
                else
                {
                    FatalErrorIn
                    (
                        "template<class Type>\n"
                        "void Foam::mechanicalModel::mapSubMeshField<Type> "
                        "subMeshes\n"
                        "(\n"
                        "    const label meshI,\n"
                        "    const GeometricField<Type, fvPatchField, volMesh>&"
                        " subMeshField,\n"
                        "    GeometricField<Type, fvPatchField, volMesh>& "
                        "baseMeshField\n"
                        ") const\n"
                    )   << "What do we do with this face?" << abort(FatalError);
                }
            }
        }
    }

    baseMeshField.correctBoundaryConditions();
}


template<class Type>
void Foam::mechanicalModel::mapSubMeshSurfaceField
(
    const label meshI,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& subMeshField,
    GeometricField<Type, fvsPatchField, surfaceMesh>& baseMeshField
) const
{
    // Map iternal field from sub-mesh to global mesh

    const labelList& faceMap = subMeshes()[meshI].faceMap();

    forAll(subMeshField, faceI)
    {
        baseMeshField[faceMap[faceI]] = subMeshField[faceI];
    }

    // Map boundary field

    const labelList& patchMap = subMeshes()[meshI].patchMap();

    forAll(subMeshField.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& subMeshFieldP =
            subMeshField.boundaryField()[patchI];

        const label start = subMeshFieldP.patch().patch().start();

        if (patchMap[patchI] != -1)
        {
            fvsPatchField<Type>& baseMeshFieldP =
                baseMeshField.boundaryField()[patchMap[patchI]];

            // Note: unlike volFields, we do map on the coupled patches for
            // surface fields
            forAll(subMeshFieldP, faceI)
            {
                const label globalGlobalMeshFace = faceMap[start + faceI];

                const label curGlobalMeshPatchFace =
                    globalGlobalMeshFace
                  - mesh().boundaryMesh()[patchMap[patchI]].start();

                baseMeshFieldP[curGlobalMeshPatchFace] = subMeshFieldP[faceI];
            }
        }
        else // interface faces shared by two materials
        {
            forAll(subMeshFieldP, faceI)
            {
                const label globalGlobalMeshFace = faceMap[start + faceI];

                if (globalGlobalMeshFace < mesh().nInternalFaces())
                {
                    // Face value will be the average from both sides
                    baseMeshField[globalGlobalMeshFace] +=
                        0.5*subMeshFieldP[faceI];
                }
                else
                {
                    const label curPatch =
                        mesh().boundaryMesh().whichPatch
                        (
                            globalGlobalMeshFace
                        );

                    const label curPatchFace =
                        globalGlobalMeshFace
                      - mesh().boundaryMesh()[curPatch].start();

                    baseMeshField.boundaryField()[curPatch][curPatchFace] =
                        subMeshFieldP[faceI];
                }
            }
        }
    }

    baseMeshField.correctBoundaryConditions();

    // Make sure the field is consistent across processor patches
    forAll(baseMeshField.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& baseMeshFieldP =
            baseMeshField.boundaryField()[patchI];

        if (baseMeshFieldP.type() == processorFvsPatchField<Type>::typeName)
        {
            const Field<Type>& patchField = baseMeshFieldP;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            OPstream::write
            (
                Pstream::blocking,
                procPatch.neighbProcNo(),
                reinterpret_cast<const char*>(patchField.begin()),
                patchField.byteSize()
            );
        }
    }

    forAll(baseMeshField.boundaryField(), patchI)
    {
        fvsPatchField<Type>& baseMeshFieldP =
            baseMeshField.boundaryField()[patchI];

        if (baseMeshFieldP.type() == processorFvsPatchField<Type>::typeName)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh().boundaryMesh()[patchI]
                );

            Field<Type> ngbPatchField(procPatch.size(), pTraits<Type>::zero);

            IPstream::read
            (
                Pstream::blocking,
                procPatch.neighbProcNo(),
                reinterpret_cast<char*>(ngbPatchField.begin()),
                ngbPatchField.byteSize()
            );

            Field<Type>& patchField = baseMeshFieldP;

            patchField = 0.5*(patchField + ngbPatchField);
        }
    }
}


// ************************************************************************* //
