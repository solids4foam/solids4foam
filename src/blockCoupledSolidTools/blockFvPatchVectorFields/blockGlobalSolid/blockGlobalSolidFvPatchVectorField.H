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
    Foam::blockGlobalSolidFvPatchVectorField

Description
    Global patch for solidPolyMesh.

SourceFiles
    blockGlobalSolidFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef blockGlobalSolidFvPatchVectorField_H
#define blockGlobalSolidFvPatchVectorField_H

#include "fvPatchFields.H"
#include "coupledFvPatchFields.H"
#include "blockFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blockGlobalSolidFvPatch Declaration
\*---------------------------------------------------------------------------*/

class blockGlobalSolidFvPatchVectorField
:
    public coupledFvPatchField<vector>,
    public blockFvPatchVectorField
{
public:

    //- Runtime type information
    TypeName("blockGlobalSolid");


    // Constructors

        //- Construct from patch and internal field
        blockGlobalSolidFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blockGlobalSolidFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blockGlobalSolidFvPatchVectorField
        //  onto a new patch
        blockGlobalSolidFvPatchVectorField
        (
            const blockGlobalSolidFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        blockGlobalSolidFvPatchVectorField
        (
            const blockGlobalSolidFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blockGlobalSolidFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        blockGlobalSolidFvPatchVectorField
        (
            const blockGlobalSolidFvPatchVectorField&,
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
                new blockGlobalSolidFvPatchVectorField(*this, iF)
            );
        }



    // Member functions


        // Evaluation functions

            //- Return patchNeighbourField
            virtual tmp< Field< vector > > patchNeighbourField () const;

            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField &f, const direction cmpt
            ) const;

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField &psiInternal,
                scalarField &result,
                const lduMatrix &,
                const scalarField &coeffs,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const Field< vector > &,
                Field< vector > &,
                const BlockLduMatrix< vector > &,
                const CoeffField< vector > &,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

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
