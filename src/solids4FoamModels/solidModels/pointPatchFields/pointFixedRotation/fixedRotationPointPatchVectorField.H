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
    Foam::fixedRotationPointPatchVectorField

Description
    Fixed rotation boundary condition with boundary non-orthogonal
    correction for stress solvers.

    Optional: time-varying angle series.

SourceFiles
    fixedRotationPointPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Federico Mazzanti, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fixedRotationPointPatchVectorField_H
#define fixedRotationPointPatchVectorField_H

#ifdef OPENFOAM_NOT_EXTEND
    #include "fixedValuePointPatchField.H"
    #include "pointPatchFieldMapper.H"
#else
    #include "fixedValuePointPatchFields.H"
    #include "PointPatchFieldMapper.H"
#endif
#include "interpolationTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#ifdef OPENFOAM_NOT_EXTEND
    typedef fixedValuePointPatchField<vector> fixedValuePointPatchVectorField;

    #ifdef OPENFOAM_NOT_EXTEND
        typedef pointPatchFieldMapper PointPatchFieldMapper;
    #endif
#endif

/*---------------------------------------------------------------------------*\
                      Class fixedRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedRotationPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{

    // Private Data

        //- Enable non-orthogonal corrections
        const Switch nonOrthogonalCorrections_;

        //- Rotation in degrees
        scalar rotationAngle_;

        //- Rotation axis
        const vector rotationAxis_;

        //- Origin of rotation
        vector rotationOrigin_;

        //- Initial face centres
        //const vectorField origFaceCentres_;

        //- Initial patch points
        const vectorField origPatchPoints_;

        //- Optional angle time series
        interpolationTable<scalar> angleSeries_;

        //- Optional translation series
        interpolationTable<vector> dispSeries_;

        //- Optional origin series
        interpolationTable<vector> originSeries_;


public:

    //- Runtime type information
    TypeName("pointFixedRotation");


    // Constructors

        //- Construct from patch and internal field
        fixedRotationPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedRotationPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedRotationPointPatchVectorField
        //  onto a new patch
        fixedRotationPointPatchVectorField
        (
            const fixedRotationPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

        //- Construct as copy
        fixedRotationPointPatchVectorField
        (
            const fixedRotationPointPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new fixedRotationPointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        fixedRotationPointPatchVectorField
        (
            const fixedRotationPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new fixedRotationPointPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Evaluation functions

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
/*            virtual tmp<Field<vector> > snGrad() const;*/

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
/*            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;*/

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
