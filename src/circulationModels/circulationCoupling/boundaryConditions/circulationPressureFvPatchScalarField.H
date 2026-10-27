/*---------------------------------------------------------------------------*\
  Author: David Kelly
  License: GNU General Public License v3 (as solids4foam)

Class
    Foam::circulationPressureFvPatchScalarField

Description
    Kinematic pressure p/rho on a patch coupled to a lumped circulation
    network. The one boundary condition for every coupling algorithm: whether
    the exchange is explicit or implicit, and of what order, is chosen in the
    "coupling" block of constant/circulationProperties, not by the BC type.

        outlet
        {
            type    circulationPressure;
            value   uniform 0;
            // network circulationProperties;   // optional: which description
        }

    The patch is bound to its port by the port's "patch" entry in the network
    description. Each time the boundary is evaluated the coupler is asked to
    sample the flux and, on a new exchange, advance the network; the pressure
    returned is the one the network wants imposed there, divided by rho.

\*---------------------------------------------------------------------------*/

#ifndef circulationPressureFvPatchScalarField_H
#define circulationPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class circulationPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    word propertiesName_;

public:

    TypeName("circulationPressure");

    circulationPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    circulationPressureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    circulationPressureFvPatchScalarField
    (
        const circulationPressureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    circulationPressureFvPatchScalarField
    (
        const circulationPressureFvPatchScalarField&
    );

    circulationPressureFvPatchScalarField
    (
        const circulationPressureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new circulationPressureFvPatchScalarField(*this)
        );
    }

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new circulationPressureFvPatchScalarField(*this, iF)
        );
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
