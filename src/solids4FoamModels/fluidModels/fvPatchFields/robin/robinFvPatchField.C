/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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

#include "robinFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
robinFvPatchField<Type>::robinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    coeff0_(p.size(), 0),
    coeff1_(p.size(), 1),
    rhs_(p.size(), pTraits<Type>::zero)
{}


template<class Type>
robinFvPatchField<Type>::robinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    coeff0_("coeff0", dict, p.size()),
    coeff1_("coeff1", dict, p.size()),
    rhs_("rhs", dict, p.size())
{
    evaluate();
}


template<class Type>
robinFvPatchField<Type>::robinFvPatchField
(
    const robinFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
#ifdef OPENFOAM_ORG
    coeff0_(mapper(ptf.coeff0_)),
    coeff1_(mapper(ptf.coeff1_)),
    rhs_(mapper(ptf.rhs_))
#else
    coeff0_(ptf.coeff0_, mapper),
    coeff1_(ptf.coeff1_, mapper),
    rhs_(ptf.rhs_, mapper)
#endif
{}

#ifndef OPENFOAM_ORG
template<class Type>
robinFvPatchField<Type>::robinFvPatchField
(
    const robinFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    coeff0_(ptf.coeff0_),
    coeff1_(ptf.coeff1_),
    rhs_(ptf.rhs_)
{}
#endif

template<class Type>
robinFvPatchField<Type>::robinFvPatchField
(
    const robinFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    coeff0_(ptf.coeff0_),
    coeff1_(ptf.coeff1_),
    rhs_(ptf.rhs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void robinFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
#ifdef OPENFOAM_ORG
    m(coeff0_, coeff0_);
    m(coeff1_, coeff1_);
    m(rhs_, rhs_);
#else
    coeff0_.autoMap(m);
    coeff1_.autoMap(m);
    rhs_.autoMap(m);
#endif
}


template<class Type>
void robinFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const robinFvPatchField<Type>& mptf =
        refCast<const robinFvPatchField<Type> >(ptf);

    coeff0_.rmap(mptf.coeff0_, addr);
    coeff1_.rmap(mptf.coeff1_, addr);
    rhs_.rmap(mptf.rhs_, addr);
}


template<class Type>
void robinFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const scalarField dn(1.0/this->patch().deltaCoeffs());

    Field<Type>::operator=
    (
        dn*rhs_/(dn*coeff0_+coeff1_)
      + coeff1_*this->patchInternalField()
       /(coeff0_*dn+coeff1_)
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > robinFvPatchField<Type>::snGrad() const
{
    return (*this - this->patchInternalField())*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > robinFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        pTraits<Type>::one*coeff1_
       /(coeff1_ + coeff0_/this->patch().deltaCoeffs());
}


template<class Type>
tmp<Field<Type> > robinFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        (rhs_/this->patch().deltaCoeffs())
       /(coeff1_ + coeff0_/this->patch().deltaCoeffs());
}


template<class Type>
tmp<Field<Type> > robinFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one
        *(coeff0_/(coeff1_ + coeff0_/this->patch().deltaCoeffs()));
}


template<class Type>
tmp<Field<Type> > robinFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return rhs_/(coeff1_ + coeff0_/this->patch().deltaCoeffs());
}


template<class Type>
void robinFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
#ifdef OPENFOAM_ORG
    writeEntry(os, "coeff0", coeff0_);
    writeEntry(os, "coeff1", coeff1_);
    writeEntry(os, "rhs", rhs_);
    writeEntry(os, "value", *this);
#else
    coeff0_.writeEntry("coeff0", os);
    coeff1_.writeEntry("coeff1", os);
    rhs_.writeEntry("rhs", os);
    this->writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
