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
    Foam::solidSymmetryFvPatchVectorField

Description
    symmetry plane with non-orthogonal correction
    same as symmetryDisplacement

SourceFiles
    solidSymmetryFvPatchVectorField.C

Author
    Zeljko Tukovic FSB Zagreb
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef solidSymmetryFvPatchVectorField_H
#define solidSymmetryFvPatchVectorField_H

#include "symmetryFvPatchFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solidSymmetryFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidSymmetryFvPatchVectorField
:
    public symmetryFvPatchVectorField
{
    // Private Data

        //- Second order accuracy correction
        Switch secondOrder_;

public:

    //- Runtime type information
    TypeName("solidSymmetry");


    // Constructors

        //- Construct from patch and internal field
        solidSymmetryFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidSymmetryFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidSymmetryFvPatchVectorField
        //  onto a new patch
        solidSymmetryFvPatchVectorField
        (
            const solidSymmetryFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidSymmetryFvPatchVectorField
        (
            const solidSymmetryFvPatchVectorField&
#ifdef OPENFOAMFOUNDATION
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidSymmetryFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidSymmetryFvPatchVectorField
        (
            const solidSymmetryFvPatchVectorField&,
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
                new solidSymmetryFvPatchVectorField(*this, iF)
            );
        }

        tmp<Field<vector> > snGrad() const;

        void evaluate(const Pstream::commsTypes);

        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
