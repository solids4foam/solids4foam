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
    Foam::windkesselPressureFvPatchScalarField

Description
    Pressure outlet boundary condition based on a three-element Windkessel
    model. The model represents the downstream circulation using a proximal
    characteristic resistance, a distal resistance and a compliance.

Usage
    \verbatim
    outlet
    {
        type            windkesselPressure;
        R               <distalResistance>;
        Rch             <characteristicResistance>;
        C               <compliance>;
        rho             <density>;
        value           uniform 0;
    }
    \endverbatim

    \c R and \c Rch have dimensions of pressure per volumetric flow rate;
    \c C has dimensions of volume per pressure. The optional state entries
    \c po, \c poo, \c Qo and \c Qoo are written automatically and may be
    supplied when restarting a case.

Author
    Zeljko Tukovic, FSB Zagreb.

    Edits by Philip Cardiff, UCD.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    windkesselPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef windkesselPressureFvPatchScalarField_H
#define windkesselPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class windkesselPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class windkesselPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:

    // Private data

        //- Outlet pressure at the previous time level
        scalar po_;

        //- Outlet pressure at the time level before \c po_
        scalar poo_;

        //- Total volumetric flow rate at the previous time level
        scalar Qo_;

        //- Total volumetric flow rate at the time level before \c Qo_
        scalar Qoo_;

        //- Distal resistance of the Windkessel model
        scalar R_;

        //- Characteristic (proximal) resistance of the Windkessel model
        scalar Rch_;

        //- Compliance of the Windkessel model
        scalar C_;

        //- Fluid density used to convert between kinematic and physical
        //  pressure
        scalar rho_;

        //- Time index at which the pressure and flow-rate histories were
        //  updated
        label timeIndex_;

public:

    //- Runtime type information
    TypeName("windkesselPressure");


    // Constructors

        //- Construct from patch and internal field
        windkesselPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        windkesselPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given windkesselPressureFvPatchScalarField
        //  onto a new patch
        windkesselPressureFvPatchScalarField
        (
            const windkesselPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        windkesselPressureFvPatchScalarField
        (
            const windkesselPressureFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new windkesselPressureFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        windkesselPressureFvPatchScalarField
        (
            const windkesselPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new windkesselPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
