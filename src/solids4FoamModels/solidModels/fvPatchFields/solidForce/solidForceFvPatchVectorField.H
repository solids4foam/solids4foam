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
    solidForceFvPatchVectorField

Description
    Force boundary condition for solid solvers. This is similar to solidTraction
    except the absolute force per face is specified.

    Note: this applies the force per face, not to the entire patch.

    Obviously, for a uniform force field, the total force applied to the patch
    as the mesh is refined as the force per face stays constant.

    The primary motivation for this condition is for use FSI cases using the
    preCICE coupling library.

    The force field can be directly specified or specified via a "force"
    volVectorField, as used with preCICE.

Author
    Philip Cardiff, UCD.

SourceFiles
    solidForceFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidForceFvPatchVectorField_H
#define solidForceFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidForceFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidForceFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{

    // Private Data

        //- Force
        vectorField force_;

        //- Force field given as a vol field
        autoPtr<volVectorField> forceFieldPtr_;

        //- Current time index
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("solidForce");


    // Constructors

        //- Construct from patch and internal field
        solidForceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidForceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidForceFvPatchVectorField onto a new patch
        solidForceFvPatchVectorField
        (
            const solidForceFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidForceFvPatchVectorField
        (
            const solidForceFvPatchVectorField&
#ifdef OPENFOAMFOUNDATION
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidForceFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidForceFvPatchVectorField
        (
            const solidForceFvPatchVectorField&,
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
                new solidForceFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
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
