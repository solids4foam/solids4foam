/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam. If not, see <http://www.gnu.org/licenses/>.

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::movingDisplacementTractionFvPatchVectorField

Description
    Specifies contact traction for a traction boundary patch according to
    the Riemann solver. External boundary traction can also be specified.

Usage
    \table
        Property     | Description             | Required    | Default value
        density      | material density        | yes         |
        displacement | Maximum displacement    | yes         |
        endTime      | Displacement end time   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
	{
	    type           movingDisplacementTraction;
        density        1000;
        displacement   (0 0.01 0);
        endTime        1.0;
	}
	\endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    movingTractionFvPatchVectorField.C
    movingDisplacementTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingDisplacementTractionFvPatchVectorField_H
#define movingDisplacementTractionFvPatchVectorField_H

#include "movingTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class movingDisplacementTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingDisplacementTractionFvPatchVectorField
:
    public movingTractionFvPatchVectorField
{

private:

    // Density
    dimensionedScalar rho_;

    // Maximum displacement
    vector uMax_;

    // End time for maximum displacement
    scalar tEnd_;


public:

    // Runtime type information
    TypeName("movingDisplacementTraction");


    // Constructors

        //- Construct from patch and internal field
        movingDisplacementTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingDisplacementTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingDisplacementTractionFvPatchVectorField onto a new patch
        movingDisplacementTractionFvPatchVectorField
        (
            const movingDisplacementTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingDisplacementTractionFvPatchVectorField
        (
            const movingDisplacementTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingDisplacementTractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingDisplacementTractionFvPatchVectorField
        (
            const movingDisplacementTractionFvPatchVectorField&,
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
                new movingDisplacementTractionFvPatchVectorField(*this, iF)
            );
        }


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
                const fvPatchVectorField&,
                const labelList&
            );

        // Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //