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

InClass
    cantileverDisplacementFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "cantileverDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cantileverDisplacementFvPatchVectorField::cantileverDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    analyticalSol_(),
    curTimeIndex_(-1)
{}


Foam::cantileverDisplacementFvPatchVectorField::cantileverDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   fixedDisplacementFvPatchVectorField(p, iF),
    analyticalSol_(dict),
    curTimeIndex_(-1)
{
    Info<< "Creating " << cantileverDisplacementFvPatchVectorField::typeName
        << " patch" << endl;

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=(patchInternalField());
    }
}


Foam::cantileverDisplacementFvPatchVectorField::cantileverDisplacementFvPatchVectorField
(
    const cantileverDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}

#ifndef OPENFOAM_ORG
Foam::cantileverDisplacementFvPatchVectorField::cantileverDisplacementFvPatchVectorField
(
    const cantileverDisplacementFvPatchVectorField& ptf
)
:
    fixedDisplacementFvPatchVectorField(ptf),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}
#endif


Foam::cantileverDisplacementFvPatchVectorField::cantileverDisplacementFvPatchVectorField
(
    const cantileverDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(ptf, iF),
    analyticalSol_(ptf.analyticalSol_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::cantileverDisplacementFvPatchVectorField::
~cantileverDisplacementFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cantileverDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedDisplacementFvPatchVectorField::autoMap(m);
}


void Foam::cantileverDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void Foam::cantileverDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        // Update old quantities at the start of a new time-step
        curTimeIndex_ = db().time().timeIndex();

        // Patch face coordinates
        const vectorField& f = patch().Cf();

        // Set the face displacement field
        fvPatchField<vector>::operator==(analyticalSol_.displacement(f));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::cantileverDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementFvPatchVectorField::write(os);

    analyticalSol_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cantileverDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
