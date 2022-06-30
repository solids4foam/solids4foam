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
    fixedDisplacementZeroShearPointPatchVectorField

Description
    The point displacement version of fixedDisplacementZeroShearFvPatchField.

    Enforces a displacement condition in the patch normal direction and a zero
    traction condition in the patch tangential directions.

SourceFiles
    fixedDisplacementZeroShearPointPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementZeroShearPointPatchVectorField_H
#define fixedDisplacementZeroShearPointPatchVectorField_H

#ifdef OPENFOAMESIORFOUNDATION
    #include "slipPointPatchField.H"
    #include "pointPatchFieldMapper.H"
#else
    #include "componentMixedPointPatchFields.H"
    #include "PointPatchFieldMapper.H"
#endif
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef OPENFOAMESIORFOUNDATION
    typedef slipPointPatchField<vector> componentMixedPointPatchVectorField;
#endif

/*---------------------------------------------------------------------------* \
            Class fixedDisplacementZeroShearPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementZeroShearPointPatchVectorField
:
    public componentMixedPointPatchVectorField
{

    // Private data

        //- Displacment time series
        interpolationTable<vector> dispSeries_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("pointFixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementZeroShearPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementZeroShearPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        fixedDisplacementZeroShearPointPatchVectorField
        (
            const fixedDisplacementZeroShearPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
#ifdef OPENFOAMFOUNDATION
            const generalPointPatchFieldMapper&
#elif defined (OPENFOAMESI)
            const pointPatchFieldMapper&
#else
            const PointPatchFieldMapper&
#endif
        );

        //- Construct as copy
        fixedDisplacementZeroShearPointPatchVectorField
        (
            const fixedDisplacementZeroShearPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new fixedDisplacementZeroShearPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedDisplacementZeroShearPointPatchVectorField
        (
            const fixedDisplacementZeroShearPointPatchVectorField&,
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
                new fixedDisplacementZeroShearPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
#ifdef OPENFOAMFOUNDATION
                const generalPointPatchFieldMapper&
#elif defined (OPENFOAMESI)
                const pointPatchFieldMapper&
#else
                const PointPatchFieldMapper&
#endif
            );

            //- Reverse map the given pointPatchField onto
            //  this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        //- Insert boundary value into the internal field
        //  Fixed values need to be inserted first. In the second go,
        //  there's nothing to do
        virtual void initEvaluate
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
