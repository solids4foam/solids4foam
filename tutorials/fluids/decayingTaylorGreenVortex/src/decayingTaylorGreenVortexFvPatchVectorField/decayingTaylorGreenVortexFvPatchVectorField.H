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
    Foam::decayingTaylorGreenVortexFvPatchVectorField

Description
    Prescribed fixed value velocity condition with boundary non-orthogonal
    correction for the decaying Taylor-Green vortex problem.

    The condition also gives the exact velocity at the face quadrature points
    of the high-order least-squares schemes (quadratureValuePatchField), so
    that a k-exact reconstruction of U gets exact boundary data.

SourceFiles
    decayingTaylorGreenVortexFvPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef decayingTaylorGreenVortexFvPatchVectorField_H
#define decayingTaylorGreenVortexFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"
#include "primitivePatchInterpolation.H"
#include "quadratureValuePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class decayingTaylorGreenVortexFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public quadratureValuePatchField<vector>
{
    // Private Data

        //- Enable non-orthogonal corrections
        const Switch nonOrthogonalCorrections_;


    // Private Member Functions

        //- Return the exact velocity at the given points at the current time
        tmp<vectorField> exactVelocity(const vectorField& points) const;

public:

    //- Runtime type information
    TypeName("decayingTaylorGreenVortexVelocity");


    // Constructors

        //- Construct from patch and internal field
        decayingTaylorGreenVortexFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        decayingTaylorGreenVortexFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given decayingTaylorGreenVortexFvPatchVectorField
        //  onto a new patch
        decayingTaylorGreenVortexFvPatchVectorField
        (
            const decayingTaylorGreenVortexFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        decayingTaylorGreenVortexFvPatchVectorField
        (
            const decayingTaylorGreenVortexFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new decayingTaylorGreenVortexFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        decayingTaylorGreenVortexFvPatchVectorField
        (
            const decayingTaylorGreenVortexFvPatchVectorField&,
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
                new decayingTaylorGreenVortexFvPatchVectorField(*this, iF)
            );
        }


    // Destructors

        //- Destructor
        virtual ~decayingTaylorGreenVortexFvPatchVectorField();


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return patch-normal gradient
            //- with non-orthogonal correction regardless of whether snGrad
            //- is corrected or not
            virtual tmp<Field<vector> > snGrad() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

            //- Return the exact velocity at the face quadrature points of quad
            virtual autoPtr<CompactListList<vector>> quadratureValues
            (
                const fvMeshQuadrature& quad
            ) const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
