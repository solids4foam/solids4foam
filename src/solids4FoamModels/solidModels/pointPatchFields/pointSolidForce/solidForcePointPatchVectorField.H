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
    solidForcePointPatchVectorField

Description
    The point displacement version of solidForceFvPatchField.

    Force boundary condition for solid solvers. This is similar to solidTraction
    except the absolute force per point is specified.

    Note: this applies the force per point, not to the entire patch.

    Obviously, for a uniform force field, the total force applied to the patch
    increases as the mesh is refined, as the force per point stays constant.

    The primary motivation for this condition is for use FSI cases using the
    preCICE coupling library.

    The force field can be directly specified or specified via a "force"
    volVectorField, as used with preCICE.

SourceFiles
    solidForcePointPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef solidForcePointPatchVectorField_H
#define solidForcePointPatchVectorField_H

#include "solidTractionPointPatchVectorField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
#ifdef OPENFOAM_NOT_EXTEND
    typedef pointPatchFieldMapper PointPatchFieldMapper;
#endif

/*---------------------------------------------------------------------------*\
            Class solidForcePointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidForcePointPatchVectorField
:
    public solidTractionPointPatchVectorField
{

    // Private data

        //- Force
        vectorField force_;

        //- Force field given as a vol field
        autoPtr<volVectorField> forceFieldPtr_;

        //- Current time index
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("pointSolidForce");


    // Constructors

        //- Construct from patch and internal field
        solidForcePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidForcePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        solidForcePointPatchVectorField
        (
            const solidForcePointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const PointPatchFieldMapper&
        );

#ifndef OPENFOAM_ORG
        //- Construct as copy
        solidForcePointPatchVectorField
        (
            const solidForcePointPatchVectorField&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new solidForcePointPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidForcePointPatchVectorField
        (
            const solidForcePointPatchVectorField&,
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
                new solidForcePointPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Const reference to force
            const vectorField& force() const
            {
                return force_;
            }

            //- Non-const reference to force
            vectorField& force()
            {
                return force_;
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
#ifdef OPENFOAM_NOT_EXTEND
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
