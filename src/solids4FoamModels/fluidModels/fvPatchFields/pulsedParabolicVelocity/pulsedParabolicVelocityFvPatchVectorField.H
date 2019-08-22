/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    parabolicVelocityFvPatchVectorField

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

    The parabolic profile is scaled in time according to the following function:

        U = USteady*t^2/sqrt((t1 - t^2)^2 + (t2*t)^2)

    where USteady is the long term steady state centreline velocity, t is the
    time, and t1 and t2 are scaling parameters.
    This results in a profile that ramps up to a peak value then reduces to a
    steady state value.

SourceFiles
    pulsedParabolicVelocityFvPatchVectorField.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Modified by Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pulsedParabolicVelocityFvPatchVectorField_H
#define pulsedParabolicVelocityFvPatchVectorField_H

#include "foamTime.H"
#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class pulsedParabolicVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Steady state centreline velocity magnitude
        scalar steadyValue_;

        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;

        //- Time scaling parameter 1
        scalar t1_;

        //- Time scaling parameter 2
        scalar t2_;

        //- Bound box min (must be specified duo to bug in parallel run)
        vector boundBoxMin_;

        //- Bound box max (must be specified duo to bug in parallel run)
        vector boundBoxMax_;

public:

    //- Runtime type information
    TypeName("pulsedParabolicVelocity");


    // Constructors

        //- Construct from patch and internal field
        pulsedParabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pulsedParabolicVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given 
        //  pulsedParabolicVelocityFvPatchVectorField
        //  onto a new patch
        pulsedParabolicVelocityFvPatchVectorField
        (
            const pulsedParabolicVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pulsedParabolicVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pulsedParabolicVelocityFvPatchVectorField
        (
            const pulsedParabolicVelocityFvPatchVectorField&,
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
                new pulsedParabolicVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Return max value
        scalar& steadyValue()
        {
            return steadyValue_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
