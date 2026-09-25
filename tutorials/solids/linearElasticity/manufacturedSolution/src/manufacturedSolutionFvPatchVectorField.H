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
    Foam::manufacturedSolutionFvPatchVectorField

Description
    Applies a displacement field from a manufactured solution as a fixed value
    condition.

    The manufactured solution comes from the manufacturedSolution class.

SourceFiles
    manufacturedSolutionFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef manufacturedSolutionFvPatchVectorField_H
#define manufacturedSolutionFvPatchVectorField_H

#include "fixedDisplacementFvPatchVectorField.H"
#include "manufacturedSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class manufacturedSolutionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class manufacturedSolutionFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
    // Private Data

        //- Manufactured solution object pointer
        autoPtr<manufacturedSolution> mmsPtr_;

        //- A copy of the constructor dict
        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("manufacturedSolution");


    // Constructors

        //- Construct from patch and internal field
        manufacturedSolutionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        manufacturedSolutionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given manufacturedSolutionFvPatchVectorField
        //  onto a new patch
        manufacturedSolutionFvPatchVectorField
        (
            const manufacturedSolutionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        manufacturedSolutionFvPatchVectorField
        (
            const manufacturedSolutionFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new manufacturedSolutionFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        manufacturedSolutionFvPatchVectorField
        (
            const manufacturedSolutionFvPatchVectorField&,
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
                new manufacturedSolutionFvPatchVectorField(*this, iF)
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
