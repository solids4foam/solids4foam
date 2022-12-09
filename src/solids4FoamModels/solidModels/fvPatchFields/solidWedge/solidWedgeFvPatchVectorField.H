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
    Foam::solidWedgeFvPatchVectorField

Description
    Standard wedge boundary with non-orthogonal correction enabled.

SourceFiles
    solidWedgeFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef solidWedgeFvPatchVectorField_H
#define solidWedgeFvPatchVectorField_H

#include "wedgeFvPatchFields.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wedgeFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidWedgeFvPatchVectorField
:
    public wedgeFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("solidWedge");


    // Constructors

        //- Construct from patch and internal field
        solidWedgeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidWedgeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidWedgeFvPatchVectorField onto a new
        //- patch
        solidWedgeFvPatchVectorField
        (
            const solidWedgeFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidWedgeFvPatchVectorField
        (
            const solidWedgeFvPatchVectorField&
#ifdef OPENFOAMFOUNDATION
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidWedgeFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidWedgeFvPatchVectorField
        (
            const solidWedgeFvPatchVectorField&,
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
                new solidWedgeFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<vector> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
#ifdef OPENFOAMESIORFOUNDATION
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
            );

            //- Return face-gradient transform diagonal
            //virtual tmp<Field<vector> > snGradTransformDiag() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
