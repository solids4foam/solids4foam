/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    fixedNormalDisplacementIncrementFvPatchVectorField

Description

SourceFiles
    fixedNormalDisplacementIncrementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedNormalDisplacementIncrementFvPatchVectorField_H
#define fixedNormalDisplacementIncrementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedDisplacementFvPatchVectorField.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedNormalDisplacementIncrementFvPatchVectorField
:
    public directionMixedDisplacementFvPatchVectorField
{
    // Private data

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("fixedNormalDisplacementIncrement");


    // Constructors

        //- Construct from patch and internal field
        fixedNormalDisplacementIncrementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedNormalDisplacementIncrementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedNormalDisplacementIncrementFvPatchVectorField
        //  onto a new patch
        fixedNormalDisplacementIncrementFvPatchVectorField
        (
            const fixedNormalDisplacementIncrementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedNormalDisplacementIncrementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedNormalDisplacementIncrementFvPatchVectorField
        (
            const fixedNormalDisplacementIncrementFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedNormalDisplacementIncrementFvPatchVectorField
                (*this, iF)
            );
        }

    // Destructor

        virtual ~fixedNormalDisplacementIncrementFvPatchVectorField()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
