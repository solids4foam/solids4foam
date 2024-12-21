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

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::tractionLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a traction boundary patch according
    to the Riemann solver. External boundary traction or pressure can also
    be specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        loadingType  | nature of loading: none,| no          | none
                     | traction or pressure_   |             |
        traction     | traction vector         | no          | vector(0,0,0)
        pressure     | pressure                | no          | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    free
    {
        type           tractionLinearMomentum;
        loadingType    traction;
        traction       (0 10 0);
    }
    free
    {
        type           tractionLinearMomentum;
        loadingType    pressure;
        pressure       5.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    tractionLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionLinearMomentumFvPatchVectorField_H
#define tractionLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class tractionLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class tractionLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

private:

    // Type of loading
    word loadingType_;

    // Boundary traction
    vector t_P_;

    // Boundary pressure
    scalar p_P_;

    // //- Pressure time series
    // interpolationTable<scalar> pressureSeries_;


public:

    // Runtime type information
    TypeName("tractionLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        tractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given tractionLinearMomentumFvPatchVectorField
        //  onto a new patch
        tractionLinearMomentumFvPatchVectorField
        (
            const tractionLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionLinearMomentumFvPatchVectorField
        (
            const tractionLinearMomentumFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionLinearMomentumFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionLinearMomentumFvPatchVectorField
        (
            const tractionLinearMomentumFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionLinearMomentumFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual word& loading()
            {
                return loadingType_;
            }

            virtual vector& traction()
            {
                return t_P_;
            }

            virtual scalar& pressure()
            {
                return p_P_;
            }

        // Mapping functions

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

        // Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
