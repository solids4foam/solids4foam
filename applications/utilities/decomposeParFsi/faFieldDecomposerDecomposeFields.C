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

#include "faFieldDecomposer.H"
#include "processorFaPatchField.H"
#include "processorFaePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
faFieldDecomposer::decomposeField
(
    const GeometricField<Type, faPatchField, areaMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.internalField(), faceAddressing_);

    // Create and map the patch field values
    PtrList<faPatchField<Type> > patchFields(boundaryAddressing_.size());

    forAll (boundaryAddressing_, patchi)
    {
        if (boundaryAddressing_[patchi] >= 0)
        {
            patchFields.set
            (
                patchi,
                faPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, areaMesh>::null(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new processorFaPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, areaMesh>::null(),
                    Field<Type>
                    (
                        field.internalField(),
                        *processorAreaPatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
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


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
faFieldDecomposer::decomposeField
(
    const GeometricField<Type, faePatchField, edgeMesh>& field
) const
{
    labelList mapAddr
    (
        labelList::subList
        (
            edgeAddressing_,
            procMesh_.nInternalEdges()
        )
    );
    forAll (mapAddr, i)
    {
        mapAddr[i] -= 1;
    }

    // Create and map the internal field values
    Field<Type> internalField
    (
        field.internalField(),
        mapAddr
    );

    // Problem with addressing when a processor patch picks up both internal
    // edges and edges from cyclic boundaries. This is a bit of a hack, but
    // I cannot find a better solution without making the internal storage
    // mechanism for edgeFields correspond to the one of edges in polyMesh
    // (i.e. using slices)
    Field<Type> allEdgeField(field.mesh().nEdges());

    forAll (field.internalField(), i)
    {
        allEdgeField[i] = field.internalField()[i];
    }

    forAll (field.boundaryField(), patchi)
    {
        const Field<Type> & p = field.boundaryField()[patchi];

        const label patchStart = field.mesh().boundary()[patchi].start();

        forAll (p, i)
        {
            allEdgeField[patchStart + i] = p[i];
        }
    }

    // Create and map the patch field values
    PtrList<faePatchField<Type> > patchFields(boundaryAddressing_.size());

    forAll (boundaryAddressing_, patchi)
    {
        if (boundaryAddressing_[patchi] >= 0)
        {
            patchFields.set
            (
                patchi,
                faePatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, edgeMesh>::null(),
                    *patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new processorFaePatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, edgeMesh>::null(),
                    Field<Type>
                    (
                        allEdgeField,
                        *processorEdgePatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, faePatchField, edgeMesh> >
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
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
void faFieldDecomposer::decomposeFields
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
