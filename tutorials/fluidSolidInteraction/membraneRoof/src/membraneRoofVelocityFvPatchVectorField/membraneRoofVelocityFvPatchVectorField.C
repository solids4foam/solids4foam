/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "membraneRoofVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

membraneRoofVelocityFvPatchVectorField::
membraneRoofVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


membraneRoofVelocityFvPatchVectorField::
membraneRoofVelocityFvPatchVectorField
(
    const membraneRoofVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper)
{}


membraneRoofVelocityFvPatchVectorField::
membraneRoofVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}


#ifndef OPENFOAM_ORG
membraneRoofVelocityFvPatchVectorField::
membraneRoofVelocityFvPatchVectorField
(
    const membraneRoofVelocityFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf)
{}
#endif


membraneRoofVelocityFvPatchVectorField::
membraneRoofVelocityFvPatchVectorField
(
    const membraneRoofVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF)
{}


membraneRoofVelocityFvPatchVectorField::
~membraneRoofVelocityFvPatchVectorField()
{}


void membraneRoofVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void membraneRoofVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);
}


void membraneRoofVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Equation (29) in the reference paper is written in terms of z, but the
    // corresponding vertical coordinate in this tutorial mesh is y.
    const scalarField y(patch().Cf().component(vector::Y));
    const scalarField uy(Foam::pow(y/350.0, 0.22));

    // Current time
    const scalar t = db().time().value();
    const scalar pi = constant::mathematical::pi;

    // Time-based velocity scaling
    scalar ut = 1.0;
    if (t < 5.0)
    {
        ut = 0.5*(Foam::sin(pi*((t/5.0) - 0.5)) + 1.0);
    }
    operator==(100.0*ut*uy*vector(1, 0, 0));

    fixedValueFvPatchVectorField::updateCoeffs();
}


makePatchTypeField
(
    fvPatchVectorField,
    membraneRoofVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
