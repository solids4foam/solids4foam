/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    followerForceBeamDisplacementNRFvPatchVectorField

Description
    Implementation of a special boundary condition where the end load follows the deformation of the beam (follower load).
    -- Derived from forceBeamDisplacementNR boundary condition

SourceFiles
    followerForceBeamDisplacementNRFvPatchVectorField.C

Author
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef followerForceBeamDisplacementNRFvPatchVectorField_H
#define followerForceBeamDisplacementNRFvPatchVectorField_H

#include "fvPatchFields.H"
#include "forceBeamDisplacementNRFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class followerForceBeamDisplacementNRFvPatch Declaration
\*---------------------------------------------------------------------------*/

class followerForceBeamDisplacementNRFvPatchVectorField
:
    public forceBeamDisplacementNRFvPatchVectorField
{
    // Private Data


        //- Current follower force
	vectorField followerForce_;

        //- Current follower force direction
        // vectorField followerForceDir_;

        //- Follower force time series
	interpolationTable<vector> followerForceSeries_;

public:

    //- Runtime type information
    TypeName("followerForceBeamDisplacementNR");

    // Constructors

        //- Construct from patch and internal field
        followerForceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        followerForceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  followerForceBeamDisplacementNRFvPatchVectorField onto a new patch
        followerForceBeamDisplacementNRFvPatchVectorField
        (
            const followerForceBeamDisplacementNRFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        followerForceBeamDisplacementNRFvPatchVectorField
        (
            const followerForceBeamDisplacementNRFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new followerForceBeamDisplacementNRFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        followerForceBeamDisplacementNRFvPatchVectorField
        (
            const followerForceBeamDisplacementNRFvPatchVectorField&,
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
                new followerForceBeamDisplacementNRFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

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

        // //- Evaluate the patch field, sets Updated to false
        // virtual void evaluate
        // (
        //     const Pstream::commsTypes commsType = Pstream::blocking
        // );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
