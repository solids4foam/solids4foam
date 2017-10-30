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

#include "wedgeFvPatch.H"
#include "newWedgeFvPatchField.H"
#include "transformField.H"
#include "symmTransform.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
newWedgeFvPatchField<Type>::newWedgeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF)
{}


template<class Type>
newWedgeFvPatchField<Type>::newWedgeFvPatchField
(
    const newWedgeFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "newWedgeFvPatchField<Type>::newWedgeFvPatchField\n"
            "(\n"
            "    const newWedgeFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
newWedgeFvPatchField<Type>::newWedgeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF, dict)
{
    if (!isType<wedgeFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "newWedgeFvPatchField<Type>::newWedgeFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    evaluate();
}


template<class Type>
newWedgeFvPatchField<Type>::newWedgeFvPatchField
(
    const newWedgeFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf)
{}


template<class Type>
newWedgeFvPatchField<Type>::newWedgeFvPatchField
(
    const newWedgeFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > newWedgeFvPatchField<Type>::snGrad() const
{
    Field<Type> pif = this->patchInternalField();
    return
    (
        transform(refCast<const wedgeFvPatch>(this->patch()).cellT(), pif) - pif
    )*(0.5*this->patch().deltaCoeffs());
}


template<class Type>
void newWedgeFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::operator==
    (
        transform
        (
            refCast<const wedgeFvPatch>(this->patch()).faceT(),
            this->patchInternalField()
        )
    );
}


template<class Type>
tmp<Field<Type> > newWedgeFvPatchField<Type>::snGradTransformDiag() const
{
    const diagTensor diagT =
        0.5*diag(I - refCast<const wedgeFvPatch>(this->patch()).cellT());

    const vector diagV(diagT.xx(), diagT.yy(), diagT.zz());

    return tmp<Field<Type> >
    (
        new Field<Type>
        (
            this->size(),
            transformMask<Type>
            (
                pow
                (
                    diagV,
                    pTraits<typename powProduct<vector, pTraits<Type>::rank>
                    ::type>::zero
                )
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
