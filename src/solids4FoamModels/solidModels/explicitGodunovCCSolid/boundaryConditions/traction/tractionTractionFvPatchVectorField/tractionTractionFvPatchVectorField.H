/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::tractionTractionFvPatchVectorField

Description
    Specifies contact traction for a traction boundary patch according to
    the Riemann solver. External boundary traction or pressure can also be
    specified according to type of loading.

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
	    type           tractionTraction;
	    loadingType    traction;
        traction       (0 10 0);
	}
    free
    {
        type           tractionTraction;
        loadingType    pressure;
        pressure       5.0;
    }
	\endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    tractionTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionTractionFvPatchVectorField_H
#define tractionTractionFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class tractionTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class tractionTractionFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

private:

    // Type of loading
    word    loadingType_;

    // Boundary traction
    vector t_P_;

    // Boundary pressure
    scalar p_P_;

    // //- Pressure time series
    // interpolationTable<scalar> pressureSeries_;

public:

    //- Runtime type information
    TypeName("tractionTraction");


    // Constructors

        //- Construct from patch and internal field
        tractionTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given tractionTractionFvPatchVectorField
        //  onto a new patch
        tractionTractionFvPatchVectorField
        (
            const tractionTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionTractionFvPatchVectorField
        (
            const tractionTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionTractionFvPatchVectorField
        (
            const tractionTractionFvPatchVectorField&,
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
                new tractionTractionFvPatchVectorField(*this, iF)
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

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

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
