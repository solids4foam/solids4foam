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
    Foam::pressureInletOutletDdtCorrVelocityFvPatchVectorField

Description
    pressureInletOutletVelocity that reports fixesValue() false, so that
    ddtScheme::fvcDdtPhiCoeff keeps the ddtCorr flux correction on this
    patch. The parent inherits fixesValue() true from directionMixed, which
    zeroes the ddtCorr coefficient on the patch even though its normal flux
    comes from the pressure equation. Behaviour is otherwise identical.

    Used on the valve planes of a CFD domain coupled to a circulation network
    (circulationPressure on p), where the flux through the patch is the
    coupled quantity.

Usage
    \verbatim
    outlet
    {
        type            pressureInletOutletDdtCorrVelocity;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    pressureInletOutletDdtCorrVelocityFvPatchVectorField.C

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pressureInletOutletDdtCorrVelocityFvPatchVectorField_H
#define pressureInletOutletDdtCorrVelocityFvPatchVectorField_H

#include "pressureInletOutletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class pressureInletOutletDdtCorrVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureInletOutletDdtCorrVelocityFvPatchVectorField
:
    public pressureInletOutletVelocityFvPatchVectorField
{
public:

    //- Runtime type information
    TypeName("pressureInletOutletDdtCorrVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureInletOutletDdtCorrVelocityFvPatchVectorField
        (
            const fvPatch& p,
            const DimensionedField<vector, volMesh>& iF
        )
        :
            pressureInletOutletVelocityFvPatchVectorField(p, iF)
        {}

        //- Construct from patch, internal field and dictionary
        pressureInletOutletDdtCorrVelocityFvPatchVectorField
        (
            const fvPatch& p,
            const DimensionedField<vector, volMesh>& iF,
            const dictionary& dict
        )
        :
            pressureInletOutletVelocityFvPatchVectorField(p, iF, dict)
        {}

        //- Construct by mapping onto a new patch
        pressureInletOutletDdtCorrVelocityFvPatchVectorField
        (
            const pressureInletOutletDdtCorrVelocityFvPatchVectorField& ptf,
            const fvPatch& p,
            const DimensionedField<vector, volMesh>& iF,
            const fvPatchFieldMapper& mapper
        )
        :
            pressureInletOutletVelocityFvPatchVectorField(ptf, p, iF, mapper)
        {}

        //- Copy construct
        pressureInletOutletDdtCorrVelocityFvPatchVectorField
        (
            const pressureInletOutletDdtCorrVelocityFvPatchVectorField& ptf
        )
        :
            pressureInletOutletVelocityFvPatchVectorField(ptf)
        {}

        //- Copy construct setting internal field reference
        pressureInletOutletDdtCorrVelocityFvPatchVectorField
        (
            const pressureInletOutletDdtCorrVelocityFvPatchVectorField& ptf,
            const DimensionedField<vector, volMesh>& iF
        )
        :
            pressureInletOutletVelocityFvPatchVectorField(ptf, iF)
        {}

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureInletOutletDdtCorrVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureInletOutletDdtCorrVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- False: the normal flux is set by the pressure equation, so the
        //  ddtCorr flux correction is kept on this patch
        virtual bool fixesValue() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
