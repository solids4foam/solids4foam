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
    Foam::extrapolatedFvPatchField

Description
    Foam::extrapolatedFvPatchField

SourceFiles
    extrapolatedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedFvPatchField_H
#define extrapolatedFvPatchField_H

#include "fixedGradientFvPatchField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class extrapolatedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class extrapolatedFvPatchField
:
    public fixedGradientFvPatchField<Type>
{
    // Private data

        //- Interpolation points
        PtrList<DynamicList<vector> > iPoints_;

        //- Return zero gradient
        Switch zeroGradient_;

public:

    //- Runtime type information
    TypeName("extrapolated");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given extrapolatedFvPatchField
        //  onto a new patch
        extrapolatedFvPatchField
        (
            const extrapolatedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedFvPatchField
        (
            const extrapolatedFvPatchField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new extrapolatedFvPatchField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        extrapolatedFvPatchField
        (
            const extrapolatedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new extrapolatedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

        // Return defining fields


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

            //- Evaluate the patch field
            virtual void evaluate
            (
#ifdef OPENFOAM_NOT_EXTEND
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

#ifdef NoRepository
#   include "extrapolatedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
