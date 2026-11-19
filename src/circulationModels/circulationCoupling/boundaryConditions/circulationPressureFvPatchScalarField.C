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
    propertiesName_("circulationProperties"),
    phiName_("phi")
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
    ),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
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
    propertiesName_(ptf.propertiesName_),
    phiName_(ptf.phiName_)
{}


Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const circulationPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    propertiesName_(ptf.propertiesName_),
    phiName_(ptf.phiName_)
{}


Foam::circulationPressureFvPatchScalarField::
circulationPressureFvPatchScalarField
(
    const circulationPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    propertiesName_(ptf.propertiesName_),
    phiName_(ptf.phiName_)
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
        circulationCoupler::lookupOrCreate(mesh, propertiesName_, phiName_);

    // Sample the fluxes and, on a new exchange, advance the network. Repeated
    // evaluations within one exchange reuse the result.
    coupler.sampleAndAdvance(mesh);

    // The network works in Pa: divide by rho only for a kinematic p
    const dimensionSet& pDims = internalField().dimensions();
    scalar pScale = 1;

    if (pDims == dimPressure/dimDensity)
    {
        pScale = 1/coupler.rho();
    }
    else if (pDims != dimPressure)
    {
        FatalErrorInFunction
            << "The field " << internalField().name() << " on patch "
            << patch().name() << " has dimensions " << pDims
            << ", but a pressure " << dimPressure
            << " or a kinematic pressure " << dimPressure/dimDensity
            << " is required" << exit(FatalError);
    }

    operator==
    (
        pScale*coupler.patchPressure(coupler.portOfPatch(patch().index()))
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::circulationPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>
    (
        "network", "circulationProperties", propertiesName_
    );
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
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
