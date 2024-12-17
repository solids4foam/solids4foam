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
    thermalRobinFvPatchScalarField

Description
    Robin boundary condition for thermal fluid-solid interaction problem.
    For example, one can apply dirichlet bc by setting
    "dirichlet_ = true" and "neumann_ = false".

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    thermalRobinFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalRobinFvPatchScalarField_H
#define thermalRobinFvPatchScalarField_H

#include "robinFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class elasticWallPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class thermalRobinFvPatchScalarField
:
    public robinFvPatchField<scalar>
{
    // Private data

        //- Use Neumann bc
        Switch neumann_;

        //- Use dirichlet bc
        Switch dirichlet_;

        //- Neighbour patch temperature
        scalarField neiTemperature_;

        //- Neighbour patch heat flux
        scalarField neiHeatFlux_;

        //- Equivalent interface heat transfer coefficient (delta/lambda)_eq
        scalarField eqInterHeatTransferCoeff_;

        //- Name of the thermal conducitivity field
        word lambdaName_;

public:

    //- Runtime type information
    TypeName("thermalRobin");


    // Constructors

        //- Construct from patch and internal field
        thermalRobinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalRobinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given thermalRobinFvPatchScalarField
        //  onto a new patch
        thermalRobinFvPatchScalarField
        (
            const thermalRobinFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalRobinFvPatchScalarField
        (
            const thermalRobinFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalRobinFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        thermalRobinFvPatchScalarField
        (
            const thermalRobinFvPatchScalarField&,
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
                new thermalRobinFvPatchScalarField(*this, iF)
            );
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

            //- Get neighbour patch temperature field
            const scalarField& temperature() const
            {
                return neiTemperature_;
            }

            //- Get neighbour patch temperature field
            scalarField& temperature()
            {
                return neiTemperature_;
            }

            //- Get neighbour patch heat flux
            const scalarField& heatFlux() const
            {
                return neiHeatFlux_;
            }

            //- Get neighbour patch heat flux
            scalarField& heatFlux()
            {
                return neiHeatFlux_;
            }

            //- Get eq interface heat transfer coefficient
            const scalarField& eqInterHeatTransferCoeff() const
            {
                return eqInterHeatTransferCoeff_;
            }

            //- Get eq interface heat transfer coefficient
            scalarField& eqInterHeatTransferCoeff()
            {
                return eqInterHeatTransferCoeff_;
            }

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
