/*---------------------------------------------------------------------------*\
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
    solidDirectionMixedFvPatchVectorField

Description
    directionMixed with non-orthogonal correction for the diffusion term

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

SourceFiles
    solidDirectionMixedFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidDirectionMixedFvPatchVectorField_H
#define solidDirectionMixedFvPatchVectorField_H

#include "directionMixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidDirectionMixedFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Limiter coefficient for snGrad correction
        scalar limitCoeff_;


public:

    //- Runtime type information
    TypeName("solidDirectionMixed");


    // Constructors

        //- Construct from patch and internal field
        solidDirectionMixedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidDirectionMixedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidDirectionMixedFvPatchVectorField
        //  onto a new patch
        solidDirectionMixedFvPatchVectorField
        (
            const solidDirectionMixedFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

#ifdef OPENFOAMFOUNDATION
        solidDirectionMixedFvPatchVectorField
        (
            const solidDirectionMixedFvPatchVectorField&
        ) = delete;
#else
        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidDirectionMixedFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidDirectionMixedFvPatchVectorField
        (
            const solidDirectionMixedFvPatchVectorField&,
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
                new solidDirectionMixedFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~solidDirectionMixedFvPatchVectorField()
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

            //- evaluate the patch
            virtual void evaluate(const Pstream::commsTypes);

            //- snGrad with non-orthogonal correction
            tmp<Foam::Field<vector> > snGrad() const;


        //- Write

            virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
