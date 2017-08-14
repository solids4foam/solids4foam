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

#include "tetPointFieldDecomposer.H"
#include "processorTetPolyPatchFields.H"
#include "globalTetPolyPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
tetPointFieldDecomposer::decomposeField
(
    const GeometricField<Type, tetPolyPatchField, tetPointMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.internalField(), directAddressing());

    // Create and map the patch field values
    PtrList<tetPolyPatchField<Type> > patchFields
    (
        boundaryAddressing_.size() + 1
    );

    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFields.set
            (
                patchI,
                tetPolyPatchField<Type>::New
                (
                    field.boundaryField()
                        [boundaryAddressing_[patchI]],
                    processorMesh_.boundary()[patchI],
                    DimensionedField<Type, tetPointMesh>::null(),
                    *patchFieldDecompPtrs_[patchI]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                new ProcessorPointPatchField
                <
                    tetPolyPatchField,
                    tetPointMesh,
                    tetPolyPatch,
                    processorTetPolyPatch,
                    tetFemMatrix,
                    Type
                >
                (
                    processorMesh_.boundary()[patchI],
                    DimensionedField<Type, tetPointMesh>::null()
                )
            );
        }
    }

    // Add the global patch by hand.  This needs to be present on
    // all processors
    patchFields.set
    (
        patchFields.size() - 1,
        new GlobalPointPatchField
        <
            tetPolyPatchField,
            tetPointMesh,
            tetPolyPatch,
            globalTetPolyPatch,
            tetFemMatrix,
            Type
        >
        (
            processorMesh_.boundary().globalPatch(),
            DimensionedField<Type, tetPointMesh>::null()
        )
    );

    // Create the field for the processor
    return tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
    (
        new GeometricField<Type, tetPolyPatchField, tetPointMesh>
        (
            IOobject
            (
                field.name(),
                processorMesh_().time().timeName(),
                processorMesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            processorMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
tmp<GeometricField<Type, elementPatchField, elementMesh> >
tetPointFieldDecomposer::decomposeField
(
    const GeometricField<Type, elementPatchField, elementMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.internalField(), cellAddressing_);

    // Create and map the patch field values
    PtrList<elementPatchField<Type> > patchFields
    (
        boundaryAddressing_.size() + 1
    );

    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFields.set
            (
                patchI,
                elementPatchField<Type>::New
                (
                    field.boundaryField()
                        [boundaryAddressing_[patchI]],
                    processorMesh_.boundary()[patchI],
                    DimensionedField<Type, elementMesh>::null(),
                    *patchFieldDecompPtrs_[patchI]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                new ProcessorPointPatchField
                <
                    elementPatchField,
                    elementMesh,
                    tetPolyPatch,
                    processorTetPolyPatch,
                    DummyMatrix,
                    Type
                >
                (
                    processorMesh_.boundary()[patchI],
                    DimensionedField<Type, elementMesh>::null()
                )
            );
        }
    }

    // Add the global patch by hand.  This needs to be present on
    // all processors
    patchFields.set
    (
        patchFields.size() - 1,
        new GlobalPointPatchField
        <
            elementPatchField,
            elementMesh,
            tetPolyPatch,
            globalTetPolyPatch,
            DummyMatrix,
            Type
        >
        (
            processorMesh_.boundary().globalPatch(),
            DimensionedField<Type, elementMesh>::null()
        )
    );

    // Create the field for the processor
    return tmp<GeometricField<Type, elementPatchField, elementMesh> >
    (
        new GeometricField<Type, elementPatchField, elementMesh>
        (
            IOobject
            (
                field.name(),
                processorMesh_().time().timeName(),
                processorMesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            processorMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class GeoField>
void tetPointFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    forAll (fields, fieldI)
    {
        decomposeField(fields[fieldI])().write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
