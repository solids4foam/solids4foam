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

#include "linearSpatialDisplacementPointPatchVectorField.H"
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

linearSpatialDisplacementPointPatchVectorField::
linearSpatialDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    a_(vector::zero),
    b_(vector::zero)
{}


linearSpatialDisplacementPointPatchVectorField::
linearSpatialDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
    a_(dict.lookup("a")),
    b_(dict.lookup("b"))
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}


linearSpatialDisplacementPointPatchVectorField::
linearSpatialDisplacementPointPatchVectorField
(
    const linearSpatialDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper&
)
:
    fixedValuePointPatchVectorField(p, iF)
{}


#ifndef OPENFOAMFOUNDATION
linearSpatialDisplacementPointPatchVectorField::
linearSpatialDisplacementPointPatchVectorField
(
    const linearSpatialDisplacementPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    a_(ptf.a_),
    b_(ptf.b_)
{}
#endif


linearSpatialDisplacementPointPatchVectorField::
linearSpatialDisplacementPointPatchVectorField
(
    const linearSpatialDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    a_(ptf.a_),
    b_(ptf.b_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void linearSpatialDisplacementPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void linearSpatialDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}

void linearSpatialDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(a_ + cmptMultiply(b_, patch().localPoints()));
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void linearSpatialDisplacementPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);

    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("b")
        << b_ << token::END_STATEMENT << nl;
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    linearSpatialDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
