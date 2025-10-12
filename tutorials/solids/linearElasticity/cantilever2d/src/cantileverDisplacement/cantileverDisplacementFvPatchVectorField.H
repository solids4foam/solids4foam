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
    cantileverDisplacementFvPatchVectorField

Description
    The boundary condition applies the analytical face displacement to the
    fixed end of a simple bending cantilever test case, as described in:

    C.E. Augarde, A.J. Deeks, The use of Timoshenko’s exact solution for a
    cantilever beam in adaptive analysis. Finite Elements in Analysis and
    Design, 44, 2008, 595–601, 10.1016/j.finel.2008.01.010.

    It is assumed that the origin is at the centre of the fixed patch.

    The beam is L long and D deep (thick).


SourceFiles
    cantileverDisplacementFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cantileverDisplacementFvPatchVectorField_H
#define cantileverDisplacementFvPatchVectorField_H

#include "fixedDisplacementFvPatchVectorField.H"
#include "cantileverAnalytical.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class cantileverDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class cantileverDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
    // Private data

        //- Helper object for calculating the analytical displacements and
        //  traction
        cantileverAnalytical analyticalSol_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("cantileverDisplacement");


    // Constructors

        //- Construct from patch and internal field
        cantileverDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cantileverDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cantileverDisplacementFvPatchVectorField onto
        //  a new patch
        cantileverDisplacementFvPatchVectorField
        (
            const cantileverDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cantileverDisplacementFvPatchVectorField
        (
            const cantileverDisplacementFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new cantileverDisplacementFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        cantileverDisplacementFvPatchVectorField
        (
            const cantileverDisplacementFvPatchVectorField&,
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
                new cantileverDisplacementFvPatchVectorField(*this, iF)
            );
        }

    //- Destructor
    virtual ~cantileverDisplacementFvPatchVectorField();


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


        // Member functions

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
