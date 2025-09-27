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

#include "decayingTaylorGreenVortexFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "patchCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

decayingTaylorGreenVortexFvPatchVectorField::decayingTaylorGreenVortexFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    nonOrthogonalCorrections_(true)
{}


decayingTaylorGreenVortexFvPatchVectorField::decayingTaylorGreenVortexFvPatchVectorField
(
    const decayingTaylorGreenVortexFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_)
{}


decayingTaylorGreenVortexFvPatchVectorField::decayingTaylorGreenVortexFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    nonOrthogonalCorrections_
    (
        dict.lookupOrDefault<Switch>("nonOrthogonalCorrections", true)
    )
{
    Info<< "Creating " << type() << " boundary condition" << endl;
}

decayingTaylorGreenVortexFvPatchVectorField::decayingTaylorGreenVortexFvPatchVectorField
(
    const decayingTaylorGreenVortexFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_)
{}

decayingTaylorGreenVortexFvPatchVectorField::decayingTaylorGreenVortexFvPatchVectorField
(
    const decayingTaylorGreenVortexFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    nonOrthogonalCorrections_(pvf.nonOrthogonalCorrections_)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //


decayingTaylorGreenVortexFvPatchVectorField::
~decayingTaylorGreenVortexFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void decayingTaylorGreenVortexFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void decayingTaylorGreenVortexFvPatchVectorField::rmap
(
    const fvPatchField<vector>& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);
}


void decayingTaylorGreenVortexFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar Re = 10;
    const scalar t = db().time().value();
    const scalar pi = constant::mathematical::pi;
    const scalar A = Foam::exp(-2.0*sqr(pi)*t/Re);

    const scalarField x(patch().Cf().component(vector::X));
    const scalarField y(patch().Cf().component(vector::Y));

    const vectorField velocity
    (
        A*Foam::sin(pi*x)*Foam::cos(pi*y)*vector(1, 0, 0)
      - A*Foam::cos(pi*x)*Foam::sin(pi*y)*vector(0, 1, 0)
    );

    operator==(velocity);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> >
decayingTaylorGreenVortexFvPatchVectorField::snGrad() const
{
    if
    (
        nonOrthogonalCorrections_
     && db().foundObject<volTensorField>("grad(" + internalField().name() + ")")
    )
    {
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + internalField().name() + ")"
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            *this
          - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        WarningInFunction
            << "Non-orthogonal corrections not applied on " << patch().name()
            << endl;

        return (*this - patchInternalField())*patch().deltaCoeffs();
    }
}

tmp<Field<vector> >
decayingTaylorGreenVortexFvPatchVectorField::gradientBoundaryCoeffs() const
{
    if
    (
        nonOrthogonalCorrections_
     && db().foundObject<volTensorField>("grad(" + internalField().name() + ")")
    )
    {
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + internalField().name() + ")"
            );

        // Non-orthogonal correction vectors
        const vectorField k(patchCorrectionVectors(patch()));

        return
        (
            patch().deltaCoeffs()
           *(*this - (k & gradField.patchInternalField()))
        );
    }
    else
    {
        WarningInFunction
            << "Non-orthogonal corrections not applied on " << patch().name()
            << endl;

        return (patch().deltaCoeffs()*(*this));
    }
}

void decayingTaylorGreenVortexFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("nonOrthogonalCorrections")
        << nonOrthogonalCorrections_ << token::END_STATEMENT << nl;

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    decayingTaylorGreenVortexFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
