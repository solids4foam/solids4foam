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
    Foam::backflowStabilisedOutletVelocityFvPatchVectorField

Description
    Outlet velocity condition with the backflow-stabilised traction of
    Esmaily Moghadam et al. (2011), as used, with beta = 0.2, by the
    Hessenthaler et al. (2017) fluid-structure interaction benchmark:

        sigma & n = rho*beta*min(U & n, 0)*U

    to be used with a fixed pressure on the outlet. On outflow faces,
    U & n > 0, the traction is zero and the condition is zero gradient (a
    do-nothing outlet). On backflow faces, the traction cancels a fraction
    beta of the influx of kinetic energy, which causes backflow divergence;
    beta = 1 cancels it completely and beta = 0 gives zero gradient.

    Keeping the dominant viscous term of the traction, the condition is the
    Robin condition

        dU/dn = c*U,  c = beta*min(U & n, 0)/nu  (c <= 0)

    which is applied implicitly as a mixed condition with refValue = 0,
    refGrad = 0 and

        valueFraction = -c/(deltaCoeffs - c)

    so that U_b = U_P*deltaCoeffs/(deltaCoeffs - c). The coefficient c is
    evaluated from the velocity of the patch-adjacent cells, U_P, of the
    velocity field named by the optional "U" entry (default U), not from
    the field the condition belongs to, so that a field with the same
    boundary types, e.g. the correction field of the newtonIcoFluid
    preconditioner, gets the coefficients of the velocity.

    min(U & n, 0) = 0.5*((U & n) - |U & n|) has a kink at U & n = 0, so the
    tangential part of the traction has a derivative jump where a face
    changes between outflow and backflow. Unlike the switch of inletOutlet
    and totalPressure, the traction itself is continuous there, and Newton
    solvers converge with the exact condition in the tests so far; to keep
    the residual differentiable, |U & n| is nevertheless replaced by
    sqrt((U & n)^2 + eps^2), where the smoothing velocity eps is the
    fraction "smoothing" (default 0.01) of the mean speed of the
    patch-adjacent cells at the previous time step, so that it is fixed
    within a time step; smoothing 0 recovers the exact condition.

    The value fraction is recomputed at every evaluation, not only in
    updateCoeffs: updateCoeffs is skipped while the patch is flagged as
    updated, e.g. after an fvMatrix of the velocity has been assembled, and
    a stale value fraction spoils the finite-difference (matrix-free)
    Jacobians of Newton solvers.

    The kinematic viscosity, nu, is read from the optional "nu" entry or,
    if it is not given, from the "nu" entry of transportProperties.

    References:
        M. Esmaily Moghadam, Y. Bazilevs, T.-Y. Hsia, I. E. Vignon-Clementel,
        A. L. Marsden (2011). A comparison of outlet boundary treatments for
        prevention of backflow divergence with relevance to blood flow
        simulations. Comput. Mech. 48:277-291.

        A. Hessenthaler, O. Roehrle, D. Nordsletten (2017). Validation of a
        non-conforming monolithic fluid-structure interaction method using
        phase-contrast MRI. Int. J. Numer. Meth. Biomed. Engng. 33:e2845.

    Example of the boundary condition specification:
    \verbatim
    outlet
    {
        type            backflowStabilisedOutletVelocity;

        // Fraction of the backflow energy influx that is cancelled
        // (default 0.2); 0 gives zero gradient
        beta            0.2;

        // Smoothing velocity of min(U & n, 0), as a fraction of the mean
        // speed of the patch-adjacent cells (default 0.01); 0 gives the
        // exact, non-smooth condition
        smoothing       0.01;

        // Optional: name of the velocity field (default U)
        // U               U;

        // Optional: kinematic viscosity; read from transportProperties if
        // not given
        // nu              1e-06;

        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    backflowStabilisedOutletVelocityFvPatchVectorField.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef backflowStabilisedOutletVelocityFvPatchVectorField_H
#define backflowStabilisedOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class backflowStabilisedOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class backflowStabilisedOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Name of the velocity field, from which the coefficients are
        //  evaluated
        word UName_;

        //- Fraction of the backflow kinetic energy influx that is cancelled
        scalar beta_;

        //- Smoothing velocity of min(U & n, 0), as a fraction of the mean
        //  speed of the patch-adjacent cells at the previous time step
        scalar smoothing_;

        //- Kinematic viscosity, if given in the dictionary; otherwise it is
        //  read from transportProperties
        scalar nu_;

        //- Is nu given in the dictionary
        bool nuFromDict_;


    // Private Member Functions

        //- Return the kinematic viscosity
        scalar nu() const;

        //- Set refValue, refGrad and valueFraction from the current velocity
        void updateValueFraction();


public:

    //- Runtime type information
    TypeName("backflowStabilisedOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        backflowStabilisedOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        backflowStabilisedOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  backflowStabilisedOutletVelocityFvPatchVectorField onto a new
        //  patch
        backflowStabilisedOutletVelocityFvPatchVectorField
        (
            const backflowStabilisedOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        backflowStabilisedOutletVelocityFvPatchVectorField
        (
            const backflowStabilisedOutletVelocityFvPatchVectorField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new backflowStabilisedOutletVelocityFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        backflowStabilisedOutletVelocityFvPatchVectorField
        (
            const backflowStabilisedOutletVelocityFvPatchVectorField&,
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
                new backflowStabilisedOutletVelocityFvPatchVectorField
                (
                    *this, iF
                )
            );
        }


    // Member functions

        // Attributes

            //- Return false: the condition is a Robin condition whose value
            //  fraction is always below one, so it does not fix the value;
            //  e.g. the newtonIcoFluid Jacobian then linearises it as a
            //  Neumann (outlet) condition rather than as a fixed value
            virtual bool fixesValue() const
            {
                return false;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field, updating the value fraction from
            //  the current velocity
            virtual void evaluate
            (
#ifdef OPENFOAM_NOT_EXTEND
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
#else
                const Pstream::commsTypes commsType = Pstream::blocking
#endif
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
