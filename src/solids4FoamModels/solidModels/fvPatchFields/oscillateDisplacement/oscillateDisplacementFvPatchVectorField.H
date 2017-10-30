/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::oscillateDisplacementFvPatchVectorField

Description
    Sinusoidal oscillating displacement condition, specified by amplitude
    and frequency.

SourceFiles
    oscillateDisplacementFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef oscillateDisplacementFvPatchVectorField_H
#define oscillateDisplacementFvPatchVectorField_H

#include "fixedDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class oscillateDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class oscillateDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
    // Private Data

        //- Amplitude of the displacement sine curve
        const vector amplitude_;

        //- Frequency (in Hz) of the displacement sine curve
        const scalar frequency_;


public:

    //- Runtime type information
    TypeName("oscillateDisplacement");


    // Constructors

        //- Construct from patch and internal field
        oscillateDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oscillateDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oscillateDisplacementFvPatchVectorField
        //  onto a new patch
        oscillateDisplacementFvPatchVectorField
        (
            const oscillateDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oscillateDisplacementFvPatchVectorField
        (
            const oscillateDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new oscillateDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oscillateDisplacementFvPatchVectorField
        (
            const oscillateDisplacementFvPatchVectorField&,
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
                new oscillateDisplacementFvPatchVectorField(*this, iF)
            );
        }



    // Member functions


        // Evaluation functions

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
