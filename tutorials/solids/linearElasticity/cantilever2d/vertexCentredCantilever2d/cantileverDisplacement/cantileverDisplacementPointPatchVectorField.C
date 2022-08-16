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

#include "cantileverDisplacementPointPatchVectorField.H"
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

cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    P_(0.0),
    E_(0.0),
    nu_(0.0),
    L_(0.0),
    D_(0.0),
    I_(0.0),
    curTimeIndex_(-1)
{}


cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
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
        fixedValuePointPatchVectorField::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        fixedValuePointPatchVectorField::operator==(vector::zero);
    }
}


cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const cantileverDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(p, iF),
    P_(ptf.P_),
    E_(ptf.E_),
    nu_(ptf.nu_),
    L_(ptf.L_),
    D_(ptf.D_),
    I_(ptf.I_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


#ifndef OPENFOAMFOUNDATION
cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const cantileverDisplacementPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    P_(ptf.P_),
    E_(ptf.E_),
    nu_(ptf.nu_),
    L_(ptf.L_),
    D_(ptf.D_),
    I_(ptf.I_),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif

cantileverDisplacementPointPatchVectorField::cantileverDisplacementPointPatchVectorField
(
    const cantileverDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
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
void cantileverDisplacementPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void cantileverDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}


void cantileverDisplacementPointPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        // Patch point coordinates
        const vectorField& p = patch().localPoints();

        // Calculate point displacement field
        vectorField disp(p.size(), vector::zero);

        const scalarField x(p.component(vector::X));
        const scalarField y(p.component(vector::Y));

        // Augarde and Deeks
        disp.replace(vector::X, (P_*y/(6*E_*I_))*(2 + nu_)*(y*y - D_*D_/4.0));
        disp.replace(vector::Y, -(P_/(6*E_*I_))*3*nu_*y*y*L_);

        // // Timoshenko and Goodier form
        // // Note: I assume origin is at centre of fixed end and positive y points
        // // up
        // // Take x and y coordinates
        // const scalarField x(L_ - p.component(vector::X));
        // const scalarField y(-p.component(vector::Y));

        // // Half the depth
        // const scalar c = D_/2.0;

        // // Shear modulus
        // const scalar G = E_/(2*(1 + nu_));

        // disp.replace
        // (
        //     vector::X,
        //     -P_*
        //     (
        //       - (x*x*y/2/E_)
        //       - (nu_*y*y*y/6/E_)
        //       + (y*y*y/6/G)
        //       + (L_*L_/2/E_ - c*c/2/G)*y
        //     )/I_
        // );
        // disp.replace
        // (
        //     vector::Y,
        //     -P_*
        //     (
        //         (nu_*x*y*y/2)
        //       + (x*x*x/6)
        //       - (L_*L_*x/2)
        //       + (L_*L_*L_/3)
        //     )/E_/I_
        // );
        
        // Set patch displacement field
        fixedValuePointPatchVectorField::operator==(disp);
    }

    fixedValuePointPatchVectorField::initEvaluate(commsType);
}


void cantileverDisplacementPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);

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
    cantileverDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
