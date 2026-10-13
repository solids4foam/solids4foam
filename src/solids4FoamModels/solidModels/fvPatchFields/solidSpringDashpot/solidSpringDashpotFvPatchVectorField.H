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
    Spring-dashpot (Robin) traction condition for solid solvers, as used to
    model the pericardium and the tissue above the base in cardiac mechanics
    benchmarks (e.g. Arostica et al. 2025, eq. 1).

    The nominal traction (per unit undeformed area) is

        P.N = -(alpha u + beta du/dt)                    (normalOnly no)
        P.N = -(alpha u.N + beta du/dt.N) N              (normalOnly yes)

    where N is the reference unit normal, alpha the spring stiffness [Pa/m]
    and beta the dashpot coefficient [Pa s/m]. With normalOnly the tangential
    traction is zero.

    It is a solidTraction condition whose traction is recomputed from the
    displacement every time the coefficients are updated, so that the
    residual a PETSc SNES solid model evaluates is a consistent function of
    the unknowns. The boundary displacement is extrapolated from the adjacent
    cell, using grad(D), in the same way as the SNES path evaluates traction
    boundaries; the velocity is a backward difference over the time step.

    useUndeformedArea is always on, because the traction is nominal.

Usage
    \verbatim
    epi
    {
        type            solidSpringDashpot;
        stiffness       1e8;     // alpha [Pa/m]
        damping         5e3;     // beta  [Pa s/m]
        normalOnly      yes;
        value           uniform (0 0 0);
    }
    \endverbatim

    A pressure may be added as for solidTraction (pressure, pressureSeries or
    pressureField). traction, tractionSeries and tractionField are not
    allowed: the traction is what this condition computes.

SourceFiles
    solidSpringDashpotFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidSpringDashpotFvPatchVectorField_H
#define solidSpringDashpotFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class solidSpringDashpotFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidSpringDashpotFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{

    // Private Data

        //- Spring stiffness alpha [Pa/m]
        scalar stiffness_;

        //- Dashpot coefficient beta [Pa s/m]
        scalar damping_;

        //- Act in the normal direction only
        Switch normalOnly_;


    // Private Member Functions

        //- The dictionary handed to solidTraction: a zero traction, which is
        //  replaced on every update, and useUndeformedArea on
        static dictionary tractionDict(const dictionary& dict);

        //- Boundary displacement extrapolated from the adjacent cells
        tmp<vectorField> extrapolatedDisplacement() const;


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

        //- Construct by mapping given
        //  solidSpringDashpotFvPatchVectorField onto a new patch
        solidSpringDashpotFvPatchVectorField
        (
            const solidSpringDashpotFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidSpringDashpotFvPatchVectorField
        (
            const solidSpringDashpotFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

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
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new solidSpringDashpotFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
