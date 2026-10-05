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
    Foam::neoHookeanMMSFvPatchVectorField

Description
    Fixed displacement condition that prescribes the neo-Hookean manufactured
    displacement at the patch face centres at the current time. The solution
    parameters are read from constant/neoHookeanMMSProperties.

Usage
    \verbatim
    patch0
    {
        type            neoHookeanMMS;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    neoHookeanMMSFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanMMSFvPatchVectorField_H
#define neoHookeanMMSFvPatchVectorField_H

#include "fixedDisplacementFvPatchVectorField.H"
#include "neoHookeanMMS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class neoHookeanMMSFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanMMSFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
    // Private Data

        //- Manufactured solution object pointer
        mutable autoPtr<neoHookeanMMS> mmsPtr_;


public:

    //- Runtime type information
    TypeName("neoHookeanMMS");


    // Constructors

        //- Construct from patch and internal field
        neoHookeanMMSFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neoHookeanMMSFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given neoHookeanMMSFvPatchVectorField
        //  onto a new patch
        neoHookeanMMSFvPatchVectorField
        (
            const neoHookeanMMSFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        neoHookeanMMSFvPatchVectorField
        (
            const neoHookeanMMSFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new neoHookeanMMSFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        neoHookeanMMSFvPatchVectorField
        (
            const neoHookeanMMSFvPatchVectorField&,
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
                new neoHookeanMMSFvPatchVectorField(*this, iF)
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
