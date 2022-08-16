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

#include "cantileverTractionPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(p, iF),
    P_(0.0),
    E_(0.0),
    nu_(0.0),
    L_(0.0),
    D_(0.0),
    I_(0.0),
    curTimeIndex_(-1)
{}


cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    solidTractionPointPatchVectorField(p, iF),
    P_(readScalar(dict.lookup("P"))),
    E_(readScalar(dict.lookup("E"))),
    nu_(readScalar(dict.lookup("nu"))),
    L_(readScalar(dict.lookup("L"))),
    D_(readScalar(dict.lookup("D"))),
    I_(Foam::pow(D_, 3.0)/12.0),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        solidTractionPointPatchVectorField::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        solidTractionPointPatchVectorField::operator==(vector::zero);
    }
}


cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const cantileverTractionPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    solidTractionPointPatchVectorField(p, iF),
    P_(ptf.P_),
    E_(ptf.E_),
    nu_(ptf.nu_),
    L_(ptf.L_),
    D_(ptf.D_),
    I_(ptf.I_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


#ifndef OPENFOAMFOUNDATION
cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const cantileverTractionPointPatchVectorField& ptf
)
:
    solidTractionPointPatchVectorField(ptf),
    P_(ptf.P_),
    E_(ptf.E_),
    nu_(ptf.nu_),
    L_(ptf.L_),
    D_(ptf.D_),
    I_(ptf.I_),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif


cantileverTractionPointPatchVectorField::cantileverTractionPointPatchVectorField
(
    const cantileverTractionPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    solidTractionPointPatchVectorField(ptf, iF),
    P_(ptf.P_),
    E_(ptf.E_),
    nu_(ptf.nu_),
    L_(ptf.L_),
    D_(ptf.D_),
    I_(ptf.I_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void cantileverTractionPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    solidTractionPointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void cantileverTractionPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidTractionPointPatchVectorField::rmap(ptf, addr);
}


void cantileverTractionPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        // Patch point coordinates
        const vectorField& p = patch().localPoints();
        const scalarField x(p.component(vector::X));
        const scalarField y(p.component(vector::Y));

        // Calculate point traction field
        vectorField trac(p.size(), vector::zero);
        trac.replace(vector::X, P_*(L_ - x)*y/I_);
        trac.replace(vector::Y, -(P_/(2*I_))*(D_*D_/4 - y*y));
        
        // Set patch point traction field
        traction() = trac;
        pressure() = 0.0;
    }

    solidTractionPointPatchVectorField::initEvaluate(commsType);
}


void cantileverTractionPointPatchVectorField::write(Ostream& os) const
{
    solidTractionPointPatchVectorField::write(os);

    os.writeKeyword("P")
        << P_ << token::END_STATEMENT << nl;
    os.writeKeyword("E")
        << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("nu")
        << nu_ << token::END_STATEMENT << nl;
    os.writeKeyword("L")
        << L_ << token::END_STATEMENT << nl;
    os.writeKeyword("D")
        << D_ << token::END_STATEMENT << nl;
    os.writeKeyword("I")
        << I_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    cantileverTractionPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
