/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::timeVaryingVelocityPointPatchVectorField

Description
    This boundary condition is simialr to the timeVaryingUniformFixedValue
    boundary condition, with the one difference being that in the current
    boundary condition the value is only not zero on the first occasion that
    initEvaluate is called per time-step. In an solver where there are multiple
    calls to the mesh motion solver within a time-step, this ensures that a
    non-zero condition is only enforced once per time-step.

SourceFiles
    timeVaryingVelocityPointPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingVelocityPointPatchVectorField_H
#define timeVaryingVelocityPointPatchVectorField_H

#ifdef OPENFOAMESIORFOUNDATION
    #include "fixedValuePointPatchField.H"
#else
    #include "FixedValuePointPatchField.H"
    #include "pointPatchField.H"
    #include "pointMesh.H"
    #include "pointPatch.H"
    #include "DummyMatrix.H"
#endif
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef OPENFOAMESIORFOUNDATION
    typedef fixedValuePointPatchField<vector>
        fixedValuePointPatchVectorField;
#else
    typedef FixedValuePointPatchField
    <
        pointPatchField, pointMesh, pointPatch, DummyMatrix, vector
    > fixedValuePointPatchVectorField;
#endif

/*---------------------------------------------------------------------------*\
        Class timeVaryingVelocityPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingVelocityPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{
    // Private data

        //- Time vs. velocity interpolation series
        interpolationTable<vector> timeSeries_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("timeVaryingVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        timeVaryingVelocityPointPatchVectorField
        (
            const timeVaryingVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
#ifdef OPENFOAMESIORFOUNDATION
            const pointPatchFieldMapper&
#else
            const PointPatchFieldMapper&
#endif
        );

#ifdef OPENFOAMFOUNDATION
        //- Construct as copy
        timeVaryingVelocityPointPatchVectorField
        (
            const timeVaryingVelocityPointPatchVectorField&
        ) = delete;
#else
        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new timeVaryingVelocityPointPatchVectorField
                (
                    *this
                )
            );
        }
#endif

        //- Construct as copy setting internal field reference
        timeVaryingVelocityPointPatchVectorField
        (
            const timeVaryingVelocityPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new timeVaryingVelocityPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Destructor

        ~timeVaryingVelocityPointPatchVectorField()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
#ifdef OPENFOAMESIORFOUNDATION
                const pointPatchFieldMapper&
#else
                const PointPatchFieldMapper&
#endif
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


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
