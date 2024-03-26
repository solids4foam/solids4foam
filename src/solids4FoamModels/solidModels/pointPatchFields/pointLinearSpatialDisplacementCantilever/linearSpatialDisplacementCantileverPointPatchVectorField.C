/*---------------------------------------------------------------------------*\
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

#include "linearSpatialDisplacementCantileverPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "pointPatchFields.H"
#include "pointBoundaryMesh.H"
#include "pointMesh.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "Time.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearSpatialDisplacementCantileverPointPatchVectorField::
linearSpatialDisplacementCantileverPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    v0_(0),
    L_(0),
    t_(0)
{}


linearSpatialDisplacementCantileverPointPatchVectorField::
linearSpatialDisplacementCantileverPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF),
    v0_(2),
    L_(20),
    t_(1)
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}


linearSpatialDisplacementCantileverPointPatchVectorField::
linearSpatialDisplacementCantileverPointPatchVectorField
(
    const linearSpatialDisplacementCantileverPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper&
)
:
    fixedValuePointPatchVectorField(p, iF),
    v0_(0),
    L_(0),
    t_(0)
{}


#ifndef OPENFOAM_NOT_EXTEND
linearSpatialDisplacementCantileverPointPatchVectorField::
linearSpatialDisplacementCantileverPointPatchVectorField
(
    const linearSpatialDisplacementCantileverPointPatchVectorField& ptf
)
:
    fixedValuePointPatchVectorField(ptf),
    v0_(ptf.v0_),
    L_(ptf.L_),
    t_(ptf.t_)
{}
#endif


linearSpatialDisplacementCantileverPointPatchVectorField::
linearSpatialDisplacementCantileverPointPatchVectorField
(
    const linearSpatialDisplacementCantileverPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ptf, iF),
    v0_(ptf.v0_),
    L_(ptf.L_),
    t_(ptf.t_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
void linearSpatialDisplacementCantileverPointPatchVectorField::autoMap
(
    const PointPatchFieldMapper& m
)
{
    fixedValuePointPatchVectorField::autoMap(m);
}


// Grab the values using rmap
void linearSpatialDisplacementCantileverPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchVectorField::rmap(ptf, addr);
}

void linearSpatialDisplacementCantileverPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(-patch().pointNormals() * ((v0_*patch().localPoints().component(vector::X)/L_)*t_));

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void linearSpatialDisplacementCantileverPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);

    os.writeKeyword("v0")
        << v0_ << token::END_STATEMENT << nl;
    os.writeKeyword("L")
        << L_ << token::END_STATEMENT << nl;
    os.writeKeyword("t")
        << t_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    linearSpatialDisplacementCantileverPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
