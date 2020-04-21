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

#include "meshSmoother.H"
//#include "conservativeMeshToMesh.H"
#include "IOobjectList.H"
#include "patchToPatchInterpolation.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// template<class Type>
// void meshSmoother::MapConservativeVolFields
// (
//     const conservativeMeshToMesh& meshToMeshInterp,
//     const label method
// )
// {
//     Info<< "    MapConservativeVolFields" << endl;

//     const fvMesh& meshSource = meshToMeshInterp.srcMesh();
//     const fvMesh& meshTarget = meshToMeshInterp.tgtMesh();

//     word fieldClassName
//     (
//         GeometricField<Type, fvPatchField, volMesh>::typeName
//     );

//     HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields =
//     (
//         mesh_.thisDb().objectRegistry::template lookupClass
//         <GeometricField<Type, fvPatchField, volMesh> >()
//     );

//     for
//     (
//      typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
//             iterator fieldIter = fields.begin();
//         fieldIter != fields.end();
//         ++fieldIter
//     )
//     {
//         //if (debug)
//         {
//             Info<< "    Interpolating " << fieldIter()->name();
//         }

//         // Reference to the field
//         GeometricField<Type, fvPatchField, volMesh>& field =
//             const_cast<GeometricField<Type, fvPatchField, volMesh>&>
//             (*fieldIter());

//         // Take a copy of the patch fields as we perform a separate boundary
//         // interpolation below after interpolating the internal field
//         List<Field<Type> > prevPatchField(field.boundaryField().size());
//         forAll(prevPatchField, patchI)
//         {
//             prevPatchField[patchI] = field.boundaryField()[patchI];
//         }

//         // Compute integral of source field
//         Type intSource = gSum(meshSource.V()*field.internalField());
//         //Info<< "Integral source: " << intSource << endl;

//         // Update the field and make sure the name remains unchanged
//         const word fieldName = field.name();
//         field = meshToMeshInterp.interpolate(field, method);
//         field.rename(fieldName);

//         // Compute integral of mapped field
//         Type intTarget = gSum(meshTarget.V()*field.internalField());
//         //Info<< "Integral target: " << intTarget << endl;

//         //if (debug)
//         {
//             Info<< ", error = "
//                 << 100.0*mag(intSource - intTarget)/(mag(intSource) + SMALL)
//                 << "%" << endl;
//         }

//         // Update the boundary fields because conservativeMeshToMesh does not
//         // interpolate on the boundaries; it justs performs a simple map.
//         // We will use the first order inverse distance method in
//         // the patchToPatchInterpolation class

//         if
//         (
//             meshSource.boundaryMesh().size()
//          != meshTarget.boundaryMesh().size()
//         )
//         {
//             FatalErrorIn
//             (
//                 "template<class Type>\n"
//                 "void meshSmoother::MapConservativeVolFields\n"
//                 "(\n"
//                 "    const conservativeMeshToMesh& meshToMeshInterp,\n"
//                 "    const label method\n"
//                 ")"
//             )   << "The source and target meshes have a different number "
//                 << "of patches!" << abort(FatalError);
//         }

//         forAll(meshSource.boundaryMesh(), patchI)
//         {
//             if (meshSource.boundaryMesh()[patchI].type() != "empty")
//             {
//                 // Create the patch interpolator
//                 // Note: it would be better here to use a globalPolyPatch to
//                 // ensure consistent treatment in parallel: actually the
//               // internal field mapping is currently not parallelised so this
//                 // is fine for serial
//                 patchToPatchInterpolation patchInterp
//                 (
//                     meshSource.boundaryMesh()[patchI], // fromPatch
//                     meshTarget.boundaryMesh()[patchI]  // toPatch
//                 );

//                 field.boundaryField()[patchI] =
//                     patchInterp.faceInterpolate(prevPatchField[patchI]);
//             }
//         }
//     }
// }


template<class Type>
void meshSmoother::AdvectVolFields
(
    const surfaceScalarField& sweptVol,
    const volScalarField& volOld,
    const volScalarField& volNew,
    const wordList& ignoreFields
)
{
    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    Info<< "    AdvectVolFields: " << fieldClassName << endl;

    // Get a list of the Type volFields
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields =
    (
        mesh_.thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >()
    );

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Check if we should skip this field
        // It might be more efficient to use a hashSet but assuming there are
        // only a small number of fields to ignore, then this is fine
        bool skipField = false;
        forAll(ignoreFields, fI)
        {
            if (ignoreFields[fI] == fieldIter()->name())
            {
                skipField = true;
                break;
            }
        }

        if (skipField)
        {
            continue;
        }

        //if (debug)
        {
            Info<< "        " << fieldIter()->name() << endl;
        }

        // Reference to the field: we use const_cast to update it
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        // Take a copy of the boundary patch snGrad
        List< Field<Type> > bsngPrev(field.boundaryField().size());
        forAll(bsngPrev, patchI)
        {
            bsngPrev[patchI] = field.boundaryField()[patchI].snGrad();
        }

        // Explicitly advect the field
        // Note: we are assuming that the field is an intrinsic property i.e.
        // per unit volume
        // Note2: fvc::div divides by the current volume so we multiply by it to
        // cancel it out
        field =
            (1.0/volNew)
           *(
               field*volOld
             - volNew*fvc::div(-sweptVol, field, "advectFields")
           );

        // Extrapolate to the boundaries using the previous snGrad
        if (Switch(dict().lookup("extrapolateBoundaries")))
        {
            forAll(bsngPrev, patchI)
            {
                if (!field.boundaryField()[patchI].coupled())
                {
                    // Extrapolate to the boundary
                    field.boundaryField()[patchI] =
                        field.boundaryField()[patchI].patchInternalField()
                      + bsngPrev[patchI]
                       /mesh().boundary()[patchI].deltaCoeffs();
                }
            }
        }
        else if (Switch(dict().lookup("zeroGradientBoundaries")))
        {
            forAll(bsngPrev, patchI)
            {
                if (!field.boundaryField()[patchI].coupled())
                {
                    // Zero-gradient extrapolation to the boundary
                    field.boundaryField()[patchI] =
                        field.boundaryField()[patchI].patchInternalField();
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
