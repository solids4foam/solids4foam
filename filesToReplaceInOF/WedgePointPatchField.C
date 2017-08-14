/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "WedgePointPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
WedgePointPatchField
<PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>::
WedgePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
WedgePointPatchField
<PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>::
WedgePointPatchField
(
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const dictionary& dict
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<WedgePointPatch>(p))
    {
        FatalIOErrorIn
        (
            "WedgePointPatchField"
            "<PatchField, Mesh, PointPatch, WedgePointPatch, "
            "MatrixType, Type>::WedgePointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not wedge type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
WedgePointPatchField
<PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>::
WedgePointPatchField
(
    const WedgePointPatchField
        <PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>&,
    const PointPatch& p,
    const DimensionedField<Type, Mesh>& iF,
    const PointPatchFieldMapper&
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<WedgePointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "WedgePointPatchField"
            "<PatchField, Mesh, PointPatch, WedgePointPatch, "
            "MatrixType, Type>::WedgePointPatchField\n"
            "(\n"
            "    const WedgePointPatchField"
            "    <PatchField, Mesh, PointPatch, WedgePointPatch, "
            "MatrixType, Type>&,\n"
            "    const PointPatch& p,\n"
            "    const DimensionedField<Type, Mesh>& iF,\n"
            "    const PointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
WedgePointPatchField
<PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>::
WedgePointPatchField
(
    const WedgePointPatchField
    <PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>& ptf
)
:
    PatchField<Type>(ptf)
{}


template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
WedgePointPatchField
<PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>::
WedgePointPatchField
(
    const WedgePointPatchField
    <PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>& ptf,
    const DimensionedField<Type, Mesh>& iF
)
:
    PatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate patch field
template
<
    template<class> class PatchField,
    class Mesh,
    class PointPatch,
    class WedgePointPatch,
    template<class> class MatrixType,
    class Type
>
void
WedgePointPatchField
<PatchField, Mesh, PointPatch, WedgePointPatch, MatrixType, Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // ZT, 26/02/2017: Size of the patch could be zero in parallel runs
    if (this->patch().meshPoints().size())
    {
        // In order to ensure that the wedge patch is always flat, take the
        // normal vector from the first point
        const vector& nHat = this->patch().pointNormals()[0];

        tmp<Field<Type> > tvalues =
            transform(I - nHat*nHat, this->patchInternalField());
        const Field<Type>& values = tvalues();

        // Get internal field to insert values into
        Field<Type>& iF = const_cast<Field<Type>&>(this->internalField());

        // Get addressing
        const labelList& meshPoints = this->patch().meshPoints();

        forAll (meshPoints, pointI)
        {
            iF[meshPoints[pointI]] = values[pointI];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
