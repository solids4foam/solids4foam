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
    Foam::neoHookeanManufacturedSolutionFvPatchVectorField

Description
    Fixed displacement condition that prescribes the neo-Hookean manufactured
    displacement at the patch faces at the current time. The solution
    parameters are read from constant/neoHookeanManufacturedSolutionProperties.

Usage
    \verbatim
    patch0
    {
        type            neoHookeanManufacturedSolution;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    neoHookeanManufacturedSolutionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanManufacturedSolutionFvPatchVectorField_H
#define neoHookeanManufacturedSolutionFvPatchVectorField_H

#include "fixedDisplacementFvPatchVectorField.H"
#include "neoHookeanManufacturedSolution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class neoHookeanManufacturedSolutionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanManufacturedSolutionFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
    // Private Data

        //- Manufactured solution object pointer
        mutable autoPtr<neoHookeanManufacturedSolution> mmsPtr_;


public:

    //- Runtime type information
    TypeName("neoHookeanManufacturedSolution");


    // Constructors

        //- Construct from patch and internal field
        neoHookeanManufacturedSolutionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neoHookeanManufacturedSolutionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  neoHookeanManufacturedSolutionFvPatchVectorField onto a new patch
        neoHookeanManufacturedSolutionFvPatchVectorField
        (
            const neoHookeanManufacturedSolutionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        neoHookeanManufacturedSolutionFvPatchVectorField
        (
            const neoHookeanManufacturedSolutionFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new neoHookeanManufacturedSolutionFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        neoHookeanManufacturedSolutionFvPatchVectorField
        (
            const neoHookeanManufacturedSolutionFvPatchVectorField&,
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
                new neoHookeanManufacturedSolutionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the exact displacement at the face quadrature points
            virtual autoPtr<CompactListList<vector>> evaluateQuadrature() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
