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
	linearSpatialDisplacementCantileverPointPatchVectorField

Description
    Displacement varying linearly in space as:

        u(x) = (v0*x/L)*t

    where v0 is the initial velocity, L is the length of the cantilever, x are the
    patch face centre positions, and t is the time of the simulation

SourceFiles
	linearSpatialDisplacementCantileverPointPatchVectorField.C

Author
    Federico Mazzanti, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef linearSpatialDisplacementCantileverPointPatchVectorField_H
#define linearSpatialDisplacementCantileverPointPatchVectorField_H

#ifdef OPENFOAM_NOT_EXTEND
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

#ifdef OPENFOAM_NOT_EXTEND
    typedef fixedValuePointPatchField<vector> fixedValuePointPatchVectorField;

    #ifdef OPENFOAM_NOT_EXTEND
        typedef pointPatchFieldMapper PointPatchFieldMapper;
    #endif
#endif

/*---------------------------------------------------------------------------* \
  Class linearSpatialDisplacementCantileverPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class linearSpatialDisplacementCantileverPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{
    // Private Data

        //- Initial velocity in m/s
        const scalar v0_;

        //- Length of cantilever in m
        const scalar L_;

        //- time of the simulation
        const scalar t_;


public:

    //- Runtime type information
    TypeName("pointLinearSpatialDisplacementCantilever");


    // Constructors

        //- Construct from patch and internal field
        linearSpatialDisplacementCantileverPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        linearSpatialDisplacementCantileverPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        linearSpatialDisplacementCantileverPointPatchVectorField
        (
            const linearSpatialDisplacementCantileverPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

#ifndef OPENFOAMFOUNDATION
        //- Construct as copy
        linearSpatialDisplacementCantileverPointPatchVectorField
        (
            const linearSpatialDisplacementCantileverPointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new linearSpatialDisplacementCantileverPointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        linearSpatialDisplacementCantileverPointPatchVectorField
        (
            const linearSpatialDisplacementCantileverPointPatchVectorField&,
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
                new linearSpatialDisplacementCantileverPointPatchVectorField(*this, iF)
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
