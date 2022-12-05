/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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
    Foam::blockFixedDisplacementZeroShearFvPatchVectorField

Description
    Fixed normal displacement and zero shear traction boundary condition for
    block coupled solid solver.
    This is equivalent to a symmetry condition but is also allows a non-zero
    normal displacement

    Implements
        Un = sqr(n) & U
        and
        0 = (I - sqr(n)) & (n & gradU)

    where Un is the user-specified normal displacement (THIS MAY BE CURRENTLY
    DISABLED AND SET TO ZERO).

SourceFiles
    blockFixedDisplacementZeroShearFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockFixedDisplacementZeroShearFvPatchVectorField_H
#define blockFixedDisplacementZeroShearFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

#include "blockFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blockFixedDisplacementZeroShearFvPatch Declaration
\*---------------------------------------------------------------------------*/

class blockFixedDisplacementZeroShearFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public blockFvPatchVectorField
{
    // Private Data

        //- Displacement time series
        interpolationTable<vector> dispSeries_;

public:

    //- Runtime type information
    TypeName("blockFixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        blockFixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blockFixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  blockFixedDisplacementZeroShearFvPatchVectorField
        //  onto a new patch
        blockFixedDisplacementZeroShearFvPatchVectorField
        (
            const blockFixedDisplacementZeroShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blockFixedDisplacementZeroShearFvPatchVectorField
        (
            const blockFixedDisplacementZeroShearFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blockFixedDisplacementZeroShearFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blockFixedDisplacementZeroShearFvPatchVectorField
        (
            const blockFixedDisplacementZeroShearFvPatchVectorField&,
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
                new blockFixedDisplacementZeroShearFvPatchVectorField(*this, iF)
            );
        }



    // Member functions


        // Access

            // virtual scalar normalDisp() const
            // {
            //     return timeSeries_(this->db().time().timeOutputValue());
            // }

        // Evaluation functions

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

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
