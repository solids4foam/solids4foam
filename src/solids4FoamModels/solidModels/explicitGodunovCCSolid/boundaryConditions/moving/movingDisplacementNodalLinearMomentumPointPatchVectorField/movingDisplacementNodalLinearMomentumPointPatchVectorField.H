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

Original Authors
    Jibran Haider - explicitSolidDynamics toolkit

    Original source code retrieved from:
    https://github.com/jibranhaider/explicitSolidDynamics/

Modifications by
    Khoder Alhamwi Alshaar, IITB - Modified to implement the explicit Godunov 
    scheme for solid mechanics as part of the solids4foam project.

Class
    Foam::movingDisplacementNodalLinearMomentumPointPatchVectorField

Description
    Specifies nodal linear momentum for a moving boundary patch according
    to a bell shaped velocity profile.

Usage
    \table
        Property     | Description             | Required    | Default value
        density      | material density        | yes         |
        displacement | Maximum displacement    | yes         |
        endTime      | Displacement end time   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    top
    {
        density         1000;
        displacement    (0 0.01 0);
        endTime         1.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    movingDisplacementNodalLinearMomentum.C

\*---------------------------------------------------------------------------*/

#ifndef movingDisplacementNodalLinearMomentumPointPatchVectorField_H
#define movingDisplacementNodalLinearMomentumPointPatchVectorField_H

#include "pointPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fvsPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class movingDisplacementNodalLinearMomentumPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingDisplacementNodalLinearMomentumPointPatchVectorField
:
    public fixedValuePointPatchVectorField
{

private:

    // Density
    dimensionedScalar rho_;

    // Maximum displacement
    vector uMax_;

    // End time for maximum displacement
    scalar tEnd_;

public:

    // Runtime type information
    TypeName("movingDisplacementNodalLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingDisplacementNodalLinearMomentumPointPatchVectorField onto a
        //  new patch
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const movingDisplacementNodalLinearMomentumPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const movingDisplacementNodalLinearMomentumPointPatchVectorField&
        );

        virtual autoPtr<pointPatchVectorField> clone() const
        {
            return autoPtr<pointPatchVectorField>
            (
                new movingDisplacementNodalLinearMomentumPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        movingDisplacementNodalLinearMomentumPointPatchVectorField
        (
            const movingDisplacementNodalLinearMomentumPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        virtual autoPtr<pointPatchVectorField> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchVectorField>
            (
                new movingDisplacementNodalLinearMomentumPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const pointPatchVectorField&,
                const labelList&
            );

        // Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //