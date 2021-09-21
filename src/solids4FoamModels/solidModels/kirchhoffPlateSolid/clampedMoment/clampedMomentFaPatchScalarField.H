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
    Foam::simplySupportedMomentFaPatchField

Description
    Moment sum boundary for a clamped edge.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef clampedMomentFaPatchField_H
#define clampedMomentFaPatchField_H

#include "fixedValueFaPatchField.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class clampedMomentFaPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class clampedMomentFaPatchScalarField
:
    public fixedValueFaPatchField<scalar>
{
    // Private data

        //- Relaxation factor
        const scalar relaxFac_;


public:

    //- Runtime type information
    TypeName("clampedMoment");


    // Constructors

        //- Construct from patch and internal field
        clampedMomentFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        clampedMomentFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        clampedMomentFaPatchScalarField
        (
            const clampedMomentFaPatchScalarField&,
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        clampedMomentFaPatchScalarField
        (
            const clampedMomentFaPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<scalar> > clone() const
        {
            return tmp<faPatchField<scalar> >
            (
                new clampedMomentFaPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        clampedMomentFaPatchScalarField
        (
            const clampedMomentFaPatchScalarField&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<scalar> > clone
        (
            const DimensionedField<scalar, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<scalar> >
            (
                new clampedMomentFaPatchScalarField(*this, iF)
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
