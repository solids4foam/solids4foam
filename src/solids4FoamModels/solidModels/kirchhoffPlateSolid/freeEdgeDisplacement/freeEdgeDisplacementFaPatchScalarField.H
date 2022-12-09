/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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
    Foam::simplySupportedMomentFaPatchField

Description
    Free plate edge boundary condition for transverse displacement

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef freeEdgeDisplacementFaPatchField_H
#define freeEdgeDisplacementFaPatchField_H

#include "fixedGradientFaPatchField.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class freeEdgeDisplacementFaPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class freeEdgeDisplacementFaPatchScalarField
:
    public fixedGradientFaPatchField<scalar>
{
    // Private data

        //- Shear force
        //const scalarField shearForce_;

        //- Relaxation factor
        const scalar relaxFac_;


public:

    //- Runtime type information
    TypeName("freeEdgeDisplacement");


    // Constructors

        //- Construct from patch and internal field
        freeEdgeDisplacementFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeEdgeDisplacementFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        freeEdgeDisplacementFaPatchScalarField
        (
            const freeEdgeDisplacementFaPatchScalarField&,
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        freeEdgeDisplacementFaPatchScalarField
        (
            const freeEdgeDisplacementFaPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<scalar> > clone() const
        {
            return tmp<faPatchField<scalar> >
            (
                new freeEdgeDisplacementFaPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        freeEdgeDisplacementFaPatchScalarField
        (
            const freeEdgeDisplacementFaPatchScalarField&,
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
                new freeEdgeDisplacementFaPatchScalarField(*this, iF)
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
