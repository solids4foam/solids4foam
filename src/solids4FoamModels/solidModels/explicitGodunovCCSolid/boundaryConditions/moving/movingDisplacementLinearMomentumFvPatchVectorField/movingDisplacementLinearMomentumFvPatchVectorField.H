/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::movingDisplacementLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a moving boundary patch according
    to a bell shaped velocity profile.

Usage
    \table
        Property     | Description             | Required    | Default value
        density      | Material density        | yes         |
        displacement | Maximum displacement    | yes         |
        endTime      | Displacement end time   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    top
    {
        density         1000;
        displacement    (0 0.01 0);
        endTime         1.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    movingDisplacementLinearMomentumPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingDisplacementLinearMomentumFvPatchVectorField_H
#define movingDisplacementLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class movingDisplacementLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingDisplacementLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
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
    TypeName("movingDisplacementLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        movingDisplacementLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingDisplacementLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingDisplacementLinearMomentumFvPatchVectorField onto a new patch
        movingDisplacementLinearMomentumFvPatchVectorField
        (
            const movingDisplacementLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingDisplacementLinearMomentumFvPatchVectorField
        (
            const movingDisplacementLinearMomentumFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingDisplacementLinearMomentumFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingDisplacementLinearMomentumFvPatchVectorField
        (
            const movingDisplacementLinearMomentumFvPatchVectorField&,
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
                new movingDisplacementLinearMomentumFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        // Mapping functions

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