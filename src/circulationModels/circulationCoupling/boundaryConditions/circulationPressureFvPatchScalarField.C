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

Description
    Sets the patch pressure of the circulationPressure boundary condition
    from the coupler, and writes the boundary condition.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationPressureFvPatchScalarField.H"
#include "circulationCoupler.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    propertiesName_("circulationProperties")
{}


Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    propertiesName_
    (
        dict.getOrDefault<word>("network", "circulationProperties")
    )
{}


Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const circulationPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    propertiesName_(ptf.propertiesName_)
{}


Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const circulationPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    propertiesName_(ptf.propertiesName_)
{}


Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const circulationPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    propertiesName_(ptf.propertiesName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::circulationPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    circulationCoupler& coupler =
        circulationCoupler::lookupOrCreate(mesh, propertiesName_);

    // Sample the fluxes and, on a new exchange, advance the network. Repeated
    // evaluations within one exchange reuse the result.
    coupler.sampleAndAdvance(mesh);

    // The network works in Pa; the incompressible solvers expect p/rho
    operator==
    (
        coupler.patchPressure(coupler.portOfPatch(patch().index()))
       /coupler.rho()
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::circulationPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "network", "circulationProperties", propertiesName_
    );
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        circulationPressureFvPatchScalarField
    );
}

// ************************************************************************* //
