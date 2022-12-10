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

#include "linearSpatialDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearSpatialDisplacementFvPatchVectorField::linearSpatialDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    a_(vector::zero),
    b_(vector::zero)
{}


linearSpatialDisplacementFvPatchVectorField::linearSpatialDisplacementFvPatchVectorField
(
    const linearSpatialDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    a_(ptf.a_),
    b_(ptf.b_)
{}


linearSpatialDisplacementFvPatchVectorField::linearSpatialDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    a_(dict.lookup("a")),
    b_(dict.lookup("b"))
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}

#ifndef OPENFOAMFOUNDATION
linearSpatialDisplacementFvPatchVectorField::linearSpatialDisplacementFvPatchVectorField
(
    const linearSpatialDisplacementFvPatchVectorField& ptf
)
:
    fixedDisplacementFvPatchVectorField(ptf),
    a_(ptf.a_),
    b_(ptf.b_)
{}
#endif

linearSpatialDisplacementFvPatchVectorField::linearSpatialDisplacementFvPatchVectorField
(
    const linearSpatialDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(ptf, iF),
    a_(ptf.a_),
    b_(ptf.b_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void linearSpatialDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    totalDisp() = a_ + cmptMultiply(b_, patch().Cf());

    fixedDisplacementFvPatchVectorField::updateCoeffs();
}


void linearSpatialDisplacementFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("b")
        << b_ << token::END_STATEMENT << nl;

    fixedDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    linearSpatialDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
