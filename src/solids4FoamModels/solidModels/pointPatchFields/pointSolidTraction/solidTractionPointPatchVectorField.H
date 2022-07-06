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
    solidTractionPointPatchVectorField

Description
    The point displacement version of solidTractionFvPatchField.

    Used to enforce a traction condition for point based solid mechanics
    solvers.

SourceFiles
    solidTractionPointPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef solidTractionPointPatchVectorField_H
#define solidTractionPointPatchVectorField_H

#include "calculatedPointPatchFields.H"
#include "interpolationTable.H"
#ifdef OPENFOAMESIORFOUNDATION
    #include "pointPatchFieldMapper.H"
#else
    #include "PointPatchFieldMapper.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
#ifdef OPENFOAMFOUNDATION
    typedef generalPointPatchFieldMapper PointPatchFieldMapper;
#elif defined (OPENFOAMESI)
    typedef pointPatchFieldMapper PointPatchFieldMapper;
#endif

/*---------------------------------------------------------------------------*\
            Class solidTractionPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionPointPatchVectorField
:
    public calculatedPointPatchVectorField
{

    // Private data

        //- Traction field
        vectorField traction_;

        //- Pressure field
        scalarField pressure_;

        //- Traction time series
        interpolationTable<vector> tractionSeries_;

        //- Pressure time series
        interpolationTable<scalar> pressureSeries_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("pointSolidTraction");


    // Constructors

        //- Construct from patch and internal field
        solidTractionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        solidTractionPointPatchVectorField
        (
            const solidTractionPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        solidTractionPointPatchVectorField
        (
            const solidTractionPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new solidTractionPointPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidTractionPointPatchVectorField
        (
            const solidTractionPointPatchVectorField&,
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
                new solidTractionPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Const reference to traction
            const vectorField& traction() const
            {
                return traction_;
            }

            //- Non-const reference to traction
            vectorField& traction()
            {
                return traction_;
            }

            //- Const reference to pressure
            const scalarField& pressure() const
            {
                return pressure_;
            }

            //- Non-const reference to pressure
            scalarField& pressure()
            {
                return pressure_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
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
