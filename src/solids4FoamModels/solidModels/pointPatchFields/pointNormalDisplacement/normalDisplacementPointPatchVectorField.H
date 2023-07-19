/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
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
    normalDisplacementPointPatchVectorField

Description
    The point displacement version of normalDisplacementFvPatchField.

    Enforces a displacement condition in the patch normal direction and a zero
    traction condition in the patch tangential directions.

SourceFiles
    normalDisplacementPointPatchVectorField.C

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

\*---------------------------------------------------------------------------*/

#ifndef normalDisplacementPointPatchVectorField_H
#define normalDisplacementPointPatchVectorField_H

#ifdef OPENFOAMESIORFOUNDATION
    #include "fixedValuePointPatchField.H"
    #include "pointPatchFieldMapper.H"
#else
    #include "fixedValuePointPatchFields.H"
    #include "PointPatchFieldMapper.H"
#endif
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef OPENFOAMESIORFOUNDATION
    typedef fixedValuePointPatchField<vector> fixedValuePointPatchVectorField;

    #ifdef OPENFOAMESIORFOUNDATION
        typedef pointPatchFieldMapper PointPatchFieldMapper;
    #endif
#endif

/*---------------------------------------------------------------------------* \
            Class normalDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class normalDisplacementPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{

    // Private data

        //- Normal displacment time series
        interpolationTable<scalar> dispSeries_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("pointNormalDisplacement");


    // Constructors

        //- Construct from patch and internal field
        normalDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        normalDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        normalDisplacementPointPatchVectorField
        (
            const normalDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

#ifndef OPENFOAMFOUNDATION
        //- Construct as copy
        normalDisplacementPointPatchVectorField
        (
            const normalDisplacementPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new normalDisplacementPointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        normalDisplacementPointPatchVectorField
        (
            const normalDisplacementPointPatchVectorField&,
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
                new normalDisplacementPointPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
