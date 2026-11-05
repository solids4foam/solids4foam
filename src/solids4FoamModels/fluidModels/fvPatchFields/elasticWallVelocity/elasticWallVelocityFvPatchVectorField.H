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
    Foam::elasticWallVelocityFvPatchVectorField

Description
    Moving-wall velocity boundary condition that computes the tangential and
    normal velocity components from the normal volumetric flux.

    It must be used together with "elasticWallPressure" for correct
    reinforcement of velocity at the fluid-solid interface when Robin boundary
    conditions are used.

SourceFiles
    elasticWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef elasticWallVelocityFvPatchVectorField_H
#define elasticWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class elasticWallVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class elasticWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    //To keep track when to copy patch face centres
    label timeIndex_;

    //Old face centres of patch. Needed for (second order) time integration
    vectorField Fc_;
    vectorField oldFc_;
    vectorField oldOldFc_;

    //- Normal derivative of the normal velocity at the wall from continuity,
    //  -div_s(U_wall), instead of zero. Zero is exact only for a wall whose
    //  velocity has no surface divergence (a rigid or sliding wall); on a
    //  deforming wall it removes most of the viscous normal stress.
    //  Default yes; 'continuityNormalGradient no' restores the zero normal
    //  gradient used before (#670)
    Switch continuityNormalGradient_;


    // Private Member Functions

        //- Set Fc_ to the current patch face centres
        void setCurrentFaceCentres();

        //- Normal derivative of the normal velocity at the wall from
        //  continuity: minus the tangential divergence of the velocity, from
        //  the gradient in the wall cells
        tmp<scalarField> continuityNormalGradient
        (
            const fvPatchField<tensor>& gradU
        ) const;


public:

    //- Runtime type information
    TypeName("elasticWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        elasticWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        elasticWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given elasticWallVelocityFvPatchVectorField
        //  onto a new patch
        elasticWallVelocityFvPatchVectorField
        (
            const elasticWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        elasticWallVelocityFvPatchVectorField
        (
            const elasticWallVelocityFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new elasticWallVelocityFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        elasticWallVelocityFvPatchVectorField
        (
            const elasticWallVelocityFvPatchVectorField&,
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
                new elasticWallVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

    // Evaluation functions

        //- Return patch-normal gradient
        virtual tmp<Field<vector> > snGrad() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
