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
    analyticalSphericalCavityTractionFvPatchVectorField

Description
    Traction boundary condition for a spherical cavity in an infinite elastic
    medium subjected to a distant stress test case to enforce analytical
    tractions at the boundaries.

Author
   Philip Cardiff, UCD.

SourceFiles
    analyticalSphericalCavityTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef analyticalSphericalCavityTractionFvPatchVectorField_H
#define analyticalSphericalCavityTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class analyticalSphericalCavityTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class analyticalSphericalCavityTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
private:

    // Private Data

        //- Far field traction
        //  It is assumed that it acts in the positive Z direction
        const scalar T0_;

        //- Hole radius
        const scalar cavityR_;

        //- Poisson's ration
        const scalar nu_;


public:

    //- Runtime type information
    TypeName("analyticalSphericalCavityTraction");


    // Constructors

        //- Construct from patch and internal field
        analyticalSphericalCavityTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticalSphericalCavityTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  analyticalSphericalCavityTractionFvPatchVectorField onto a new patch
        analyticalSphericalCavityTractionFvPatchVectorField
        (
            const analyticalSphericalCavityTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        analyticalSphericalCavityTractionFvPatchVectorField
        (
            const analyticalSphericalCavityTractionFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new analyticalSphericalCavityTractionFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        analyticalSphericalCavityTractionFvPatchVectorField
        (
            const analyticalSphericalCavityTractionFvPatchVectorField&,
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
                new analyticalSphericalCavityTractionFvPatchVectorField(*this, iF)
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
