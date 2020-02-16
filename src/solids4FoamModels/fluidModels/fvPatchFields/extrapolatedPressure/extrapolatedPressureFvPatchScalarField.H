/*---------------------------------------------------------------------------* \
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
    Foam::extrapolatedPressureFvPatchScalarField

Description
    Set the pressure gradient boundary condition appropriately for buoyant flow.

    If the variable name is "pd" assume it is p - rho*g.h and set the gradient
    appropriately.  Otherwise assume the variable is the static pressure.

SourceFiles
    extrapolatedPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef extrapolatedPressureFvPatchScalarFields_H
#define extrapolatedPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class extrapolatedPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extrapolatedPressureFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
    // Private data

public:

    //- Runtime type information
    TypeName("extrapolatedPressure");


    // Constructors

        //- Construct from patch and internal field
        extrapolatedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extrapolatedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  extrapolatedPressureFvPatchScalarField onto a new patch
        extrapolatedPressureFvPatchScalarField
        (
            const extrapolatedPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        extrapolatedPressureFvPatchScalarField
        (
            const extrapolatedPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new extrapolatedPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extrapolatedPressureFvPatchScalarField
        (
            const extrapolatedPressureFvPatchScalarField&,
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
                new extrapolatedPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate the patch field
        virtual void evaluate
        (
#ifdef OPENFOAMESIORFOUNDATION
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
