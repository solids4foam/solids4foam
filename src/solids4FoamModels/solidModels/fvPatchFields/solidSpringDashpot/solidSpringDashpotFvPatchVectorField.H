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
    solidSpringDashpotFvPatchVectorField

Description
    Spring-dashpot (Robin) boundary condition for solid solvers, with
    separate normal and tangential stiffness and damping coefficients.
    The spring and dashpot are treated implicitly through solidDirectionMixed.

    Optionally, the coefficients are scaled by a weight volScalarField.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    solidSpringDashpotFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidSpringDashpotFvPatchVectorField_H
#define solidSpringDashpotFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class solidSpringDashpotFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidSpringDashpotFvPatchVectorField
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

        //- Additional traction per unit current area [Pa]
        vectorField traction_;

        //- Additional pressure [Pa]
        scalarField pressure_;

        //- Old-time displacement extrapolated to the faces, for the dashpot
        //  velocity in springDashpotTraction; computed once per time step
        mutable vectorField DoldExtrap_;

        //- Time index at which DoldExtrap_ was computed
        mutable label DoldExtrapTimeIndex_;

public:

    //- Runtime type information
    TypeName("solidSpringDashpot");


    // Constructors

        //- Construct from patch and internal field
        solidSpringDashpotFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidSpringDashpotFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidSpringDashpotFvPatchVectorField
        //  onto a new patch
        solidSpringDashpotFvPatchVectorField
        (
            const solidSpringDashpotFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

#if defined(OPENFOAM_ORG) || (OPENFOAM >= 2606)
        //- Construct as copy
        solidSpringDashpotFvPatchVectorField
        (
            const solidSpringDashpotFvPatchVectorField&
        ) = delete;
#else
        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidSpringDashpotFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidSpringDashpotFvPatchVectorField
        (
            const solidSpringDashpotFvPatchVectorField&,
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
                new solidSpringDashpotFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~solidSpringDashpotFvPatchVectorField()
        {}


    // Member functions

        // Access

            //- Spring-dashpot traction per unit current area from the current
            //  boundary value, for solid models that impose boundary tractions
            //  directly; nPressure is the normal for the pressure and
            //  areaRatio the reference-to-current face area ratio
            tmp<vectorField> springDashpotTraction
            (
                const vectorField& nPressure,
                const scalarField& areaRatio
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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Diagonal of the valueFraction, used for the implicit matrix
            //  coefficient (the base class returns its square root)
            virtual tmp<Field<vector> > snGradTransformDiag() const;

            //- A Robin condition does not fix the value. The base class
            //  (directionMixed) returns true, which would make the high-order
            //  least-squares reconstruction use the patch value as a
            //  Dirichlet datum
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Evaluate the patch. With the high-order residual, the patch
            //  value is set to the high-order reconstruction at the face
            //  centres, which is the displacement the imposed spring-dashpot
            //  traction uses; the Robin closure of the base class is built
            //  from the low-order cell stress, which the high-order residual
            //  does not update
            virtual void evaluate
            (
#ifdef OPENFOAM_NOT_EXTEND
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
#else
                const Pstream::commsTypes commsType = Pstream::blocking
#endif
            );

            //- Displacement at the patch face centres reconstructed from the
            //  adjacent cells by the high-order least-squares scheme
            tmp<vectorField> highOrderFaceCentreValue
            (
                const volVectorField& vf,
                const solidModel& solMod
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
