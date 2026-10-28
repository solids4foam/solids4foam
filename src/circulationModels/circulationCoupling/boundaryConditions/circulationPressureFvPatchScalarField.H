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

Author
    David Kelly, UCD.

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

// ************************************************************************* //
