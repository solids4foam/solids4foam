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
    elasticWallPressureFvPatchScalarField

Description
    Robin boundary condition for pressure of a flow in a fluid-solid
    interaction problem which accounts for the solid inertia.

    Reference:
    Tuković et al. Added Mass Partitioned Fluid–Structure Interaction Solver
    Based on a Robin Boundary Condition for Pressure. OpenFOAM. 2019;1–22.

SourceFiles
    elasticWallPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef elasticWallPressureFvPatchScalarField_H
#define elasticWallPressureFvPatchScalarField_H

#include "robinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class elasticWallPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class elasticWallPressureFvPatchScalarField
:
    public robinFvPatchField<scalar>
{
    // Private data
    //
        //- Pressure applied at the solid side of the interface
        //  in the previous fsi iteration
        scalarField prevPressure_;

        //- Acceleration of the solid side of the interface
        //  in the previous fsi iteration
        vectorField prevAcceleration_;

public:

    //- Runtime type information
    TypeName("elasticWallPressure");


    // Constructors

        //- Construct from patch and internal field
        elasticWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        elasticWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given elasticWallPressureFvPatchScalarField
        //  onto a new patch
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new elasticWallPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        elasticWallPressureFvPatchScalarField
        (
            const elasticWallPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new elasticWallPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            //- Get previous pressure patch field
            const scalarField& prevPressure() const
            {
                return prevPressure_;
            }

            //- Get previous pressure patch field
            scalarField& prevPressure()
            {
                return prevPressure_;
            }

            //- Get previous pressure patch field
            const vectorField& prevAcceleration() const
            {
                return prevAcceleration_;
            }

            //- Get previous pressure patch field
            vectorField& prevAcceleration()
            {
                return prevAcceleration_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Calculate patch flux
            virtual void patchFlux
            (
                GeometricField<scalar, fvsPatchField, surfaceMesh>& flux,
                const fvMatrix<scalar>& matrix
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
