/*---------------------------------------------------------------------------* \
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
    movingWallPressureFvPatchScalarField

Description
    Pressure condition for a moving no-slip wall, e.g. the fluid side of a
    Dirichlet-Neumann FSI interface. The normal pressure gradient is set from
    the wall-normal momentum equation, in which the convective term vanishes
    because the fluid moves with the wall:

        n & grad(p) = n & (viscous term - a)

    where a is the wall acceleration and, with viscousTerm yes,

        viscous term = laplacian(nuEff, U) + div(nuEff*dev2(T(grad(U))))

    is evaluated in the cells next to the wall at the current iterate, with
    nuEff from the fluid's turbulence (momentum transport) model, as in the
    fluid momentum equation; on OpenFOAM.com it includes the momentum source
    of the fvOptions (see wallNormalMomentum). The gradient is multiplied by
    the fluid density when the pressure has the dimensions of pressure; rho
    is then the rho field, if registered, or else the rho entry of
    transportProperties.

    The acceleration is that of a newMovingWallVelocity patch; for any other
    velocity condition it is the time derivative of the boundary velocity
    with the ddt scheme of U. This is the acceleration of the fluid at the
    wall only if the wall moves with the mesh, as at an FSI interface: the
    convective acceleration of a wall that moves tangentially relative to
    the mesh, e.g. a rotating wall on a static mesh, is not included.

    Enclosed domains (no boundary fixes the pressure) need no special
    treatment: the solver adds the flux of the prescribed gradient,
    rAU*snGrad(p)*|Sf|, after adjustPhi (see pimpleFluid), so that the net
    boundary flux that adjustPhi checks is that of the wall velocity.

Usage
    \verbatim
    interface
    {
        type            movingWallPressure;
        viscousTerm     yes;    // optional, default yes
        value           uniform 0;
    }
    \endverbatim

    With a fixed-pressure boundary (e.g. an outlet), viscousTerm no gives the
    condition of solids4foam v2.4 and earlier, which prescribes the
    acceleration term only; pimpleFluid then also keeps its earlier
    treatment, in which the predicted flux uses the gradient of the previous
    update of the condition.

SourceFiles
    movingWallPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef movingWallPressureFvPatchScalarFields_H
#define movingWallPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "volFieldsFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class movingWallPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class movingWallPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Include the viscous term of the wall-normal momentum equation
        Switch viscousTerm_;


    // Private Member Functions

        //- Return the wall acceleration
        tmp<vectorField> wallAcceleration() const;

        //- Return the fluid density on the patch
        tmp<scalarField> patchDensity() const;

        //- Return the normal gradient from the wall-normal momentum
        //  equation; viscousPtr is the viscous term field, or null if it is
        //  not needed
        tmp<scalarField> momentumGradient
        (
            const volVectorField* viscousPtr
        ) const;


public:

    //- Runtime type information
    TypeName("movingWallPressure");


    // Constructors

        //- Construct from patch and internal field
        movingWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingWallPressureFvPatchScalarField onto a new patch
        movingWallPressureFvPatchScalarField
        (
            const movingWallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingWallPressureFvPatchScalarField
        (
            const movingWallPressureFvPatchScalarField&
#ifdef OPENFOAM_ORG
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new movingWallPressureFvPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        movingWallPressureFvPatchScalarField
        (
            const movingWallPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new movingWallPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return true if the viscous term is included
        bool viscousTerm() const
        {
            return viscousTerm_;
        }

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
