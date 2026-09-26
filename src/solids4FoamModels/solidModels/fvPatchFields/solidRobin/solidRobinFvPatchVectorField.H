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
    solidRobinFvPatchVectorField

Description
    Spring-dashpot (Robin) boundary condition for solid solvers, with
    separate normal and tangential stiffness and damping coefficients.
    The spring and dashpot are treated implicitly through solidDirectionMixed.

    Optionally, the coefficients are scaled by a weight volScalarField.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    solidRobinFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidRobinFvPatchVectorField_H
#define solidRobinFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class solidRobinFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidRobinFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Normal spring stiffness per unit reference area [Pa/m]
        scalarField kNormal_;

        //- Tangential spring stiffness per unit reference area [Pa/m]
        scalarField kTangential_;

        //- Normal dashpot coefficient per unit reference area [Pa s/m]
        scalarField cNormal_;

        //- Tangential dashpot coefficient per unit reference area [Pa s/m]
        scalarField cTangential_;

        //- Additional traction, as in solidTraction [Pa]
        vectorField traction_;

        //- Additional pressure, as in solidTraction [Pa]
        scalarField pressure_;

public:

    //- Runtime type information
    TypeName("solidRobin");


    // Constructors

        //- Construct from patch and internal field
        solidRobinFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidRobinFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidRobinFvPatchVectorField
        //  onto a new patch
        solidRobinFvPatchVectorField
        (
            const solidRobinFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

#ifdef OPENFOAM_ORG
        //- Construct as copy
        solidRobinFvPatchVectorField
        (
            const solidRobinFvPatchVectorField&
        ) = delete;
#else
        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidRobinFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidRobinFvPatchVectorField
        (
            const solidRobinFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidRobinFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~solidRobinFvPatchVectorField()
        {}


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
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

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
