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

#include "pointFieldDecomposer.H"
#include "processorPointPatchFields.H"
#include "globalPointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, pointPatchField, pointMesh> >
pointFieldDecomposer::decomposeField
(
    const GeometricField<Type, pointPatchField, pointMesh>& field
) const
{
    // Create and map the internal field values
//     Field<Type> internalField(field.internalField(), pointAddressing_);
    Field<Type> internalField(procMesh_.nPoints(), pTraits<Type>::zero);
    const Field<Type>& fieldI = field.internalField();
    for(label i=0; i<procMesh_.nPoints(); i++)
    {
        internalField[i] = fieldI[pointAddressing_[i]];
    }


    // Create a list of pointers for the patchFields including one extra
    // for the global patch
    PtrList<pointPatchField<Type> > patchFields
    (
        boundaryAddressing_.size() + 1
    );

    // Create and map the patch field values
    forAll (boundaryAddressing_, patchi)
    {
        if (patchFieldDecomposerPtrs_[patchi])
        {
            patchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, pointMesh>::null(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new ProcessorPointPatchField
                <
                    pointPatchField,
                    pointMesh,
                    pointPatch,
                    processorPointPatch,
                    DummyMatrix,
                    Type
                >
                (
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
    }

    // Add the global patch
    patchFields.set
    (
        boundaryAddressing_.size(),
        new GlobalPointPatchField
        <
            pointPatchField,
            pointMesh,
            pointPatch,
            globalPointPatch,
            DummyMatrix,
            Type
        >
        (
            procMesh_.boundary().globalPatch(),
            DimensionedField<Type, pointMesh>::null()
        )
    );

    // Create the field for the processor
    return tmp<GeometricField<Type, pointPatchField, pointMesh> >
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_().time().timeName(),
                procMesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class GeoField>
void pointFieldDecomposer::decomposeFields
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
