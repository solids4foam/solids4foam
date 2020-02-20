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
    Foam::freeSurfaceVelocityFvPatchVectorField

Description
    Velocity outlet boundary condition for patches where the pressure is
    specified. The outflow velocity is obtained by "fixedGradient" and then
    corrected from the flux.

SourceFiles
    freeSurfaceVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurfaceVelocityFvPatchVectorField_H
#define freeSurfaceVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class freeSurfaceVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class freeSurfaceVelocityFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private data

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Current time index
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("freeSurfaceVelocity");


    // Constructors

        //- Construct from patch and internal field
        freeSurfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeSurfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given freeSurfaceVelocityFvPatchVectorField
        //  onto a new patch
        freeSurfaceVelocityFvPatchVectorField
        (
            const freeSurfaceVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new freeSurfaceVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        freeSurfaceVelocityFvPatchVectorField
        (
            const freeSurfaceVelocityFvPatchVectorField&,
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
                new freeSurfaceVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
#ifdef OPENFOAMFOUNDATION
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
#else
                const Pstream::commsTypes commsType=Pstream::blocking
#endif
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
