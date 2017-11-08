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
    Foam::blockSolidVelocityFvPatchVectorField

Description
    BlockFixedDisplacement boundary condition for coupled solid solver.
    This is equivalent to a fixedValue boundary condition.

    Implements
        U = I & Ub

    where Ub is the user specified blockSolidVelocity

SourceFiles
    blockSolidVelocityFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef blockSolidVelocityFvPatchVectorField_H
#define blockSolidVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "blockFvPatchVectorField.H"
#include "interpolationTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blockSolidVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class blockSolidVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public blockFvPatchVectorField
{
    // Private Data

        //- Velocity (constant in time)
        vectorField velocity_;

        //- Velocity time series
        interpolationTable<vector> velocitySeries_;

        //- Patch interpolator to interpolate face values to the points
        mutable primitivePatchInterpolation* interpPtr_;

    // Private Functions

        //- Make the patch interpolator
        void makeInterp() const;

        //- Return the patch interpolator
        primitivePatchInterpolation& interp();

        //- Set the point displacement field boundary condition
        void setPointDisplacement(const vectorField& faceDisp);

public:

    //- Runtime type information
    TypeName("blockSolidVelocity");


    // Constructors

        //- Construct from patch and internal field
        blockSolidVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blockSolidVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blockSolidVelocityFvPatchVectorField
        //  onto a new patch
        blockSolidVelocityFvPatchVectorField
        (
            const blockSolidVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blockSolidVelocityFvPatchVectorField
        (
            const blockSolidVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blockSolidVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blockSolidVelocityFvPatchVectorField
        (
            const blockSolidVelocityFvPatchVectorField&,
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
                new blockSolidVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        // Block functions

            //- Insert coefficients into the block matrix
            virtual void insertBlockCoeffs
            (
                const solidPolyMesh& solidMesh,
                const surfaceScalarField& muf,
                const surfaceScalarField& lambdaf,
                const GeometricField<vector, fvPatchField, volMesh>& U,
                Field<vector>& blockB,
                BlockLduMatrix<vector>& blockM
            ) const;

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

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
