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
    Foam::symmetryDisplacementFvPatchVectorField

Description
    Foam::symmetryDisplacementFvPatchVectorField

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    symmetryDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryDisplacementFvPatchVectorField_H
#define symmetryDisplacementFvPatchVectorField_H

#include "symmetryFvPatchFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class symmetryDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class symmetryDisplacementFvPatchVectorField
:
    public symmetryFvPatchVectorField
{
    // Private Data

        //- Second order accuracy correction
        Switch secondOrder_;

public:

    //- Runtime type information
    TypeName("symmetryDisplacement");


    // Constructors

        //- Construct from patch and internal field
        symmetryDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetryDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given symmetryDisplacementFvPatchVectorField
        //  onto a new patch
        symmetryDisplacementFvPatchVectorField
        (
            const symmetryDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        symmetryDisplacementFvPatchVectorField
        (
            const symmetryDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new symmetryDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        symmetryDisplacementFvPatchVectorField
        (
            const symmetryDisplacementFvPatchVectorField&,
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
                new symmetryDisplacementFvPatchVectorField(*this, iF)
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
