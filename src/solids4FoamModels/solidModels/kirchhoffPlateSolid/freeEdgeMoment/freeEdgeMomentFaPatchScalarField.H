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
    Foam::freeEdgeMomentFaPatchScalarField

Description
    Free plate edge boundary condition for the moment sum of the Kirchhoff
    plate model.

    On a free edge the effective transverse shear force vanishes,
    V_n = Q_n + dM_nt/dt = 0, where Q_n = dM/dn is the transverse shear force
    and M_nt = -D (1 - nu) n.(grad grad w).t is the twisting moment. Along
    an edge of curvature kappa the mixed derivative is
    d(dw/dn)/dt - kappa dw/dt, so the normal gradient of the moment sum is

        dM/dn = D (1 - nu) d/dt [ d(dw/dn)/dt - kappa dw/dt ]

    which reduces to the second derivative of the normal slope along a
    straight edge. The normal slope is taken from the transverse
    displacement condition on the same patch and differentiated along the
    edge. At a corner with a clamped or supported edge the slopes vanish; at
    a corner with another free edge the twisting moment vanishes and the
    slope is symmetric about the corner.

    The transverse displacement is expected to use the freeEdgeDisplacement
    condition on the same patch.

Author
    Ivan Batistic, FSB.

SourceFiles
    freeEdgeMomentFaPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef freeEdgeMomentFaPatchScalarField_H
#define freeEdgeMomentFaPatchScalarField_H

#include "fixedGradientFaPatchField.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class freeEdgeMomentFaPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class freeEdgeMomentFaPatchScalarField
:
    public fixedGradientFaPatchField<scalar>
{
    // Private data

        //- Relaxation factor
        const scalar relaxFac_;


public:

    //- Runtime type information
    TypeName("freeEdgeMoment");


    // Constructors

        //- Construct from patch and internal field
        freeEdgeMomentFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeEdgeMomentFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        freeEdgeMomentFaPatchScalarField
        (
            const freeEdgeMomentFaPatchScalarField&,
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const faPatchFieldMapper&
        );

#if (OPENFOAM >= 2606)
        //- No copy construct without an internal field
        freeEdgeMomentFaPatchScalarField
        (
            const freeEdgeMomentFaPatchScalarField&
        ) = delete;
#else
        //- Construct as copy
        freeEdgeMomentFaPatchScalarField
        (
            const freeEdgeMomentFaPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<scalar> > clone() const
        {
            return tmp<faPatchField<scalar> >
            (
                new freeEdgeMomentFaPatchScalarField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        freeEdgeMomentFaPatchScalarField
        (
            const freeEdgeMomentFaPatchScalarField&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<scalar> > clone
        (
            const DimensionedField<scalar, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<scalar> >
            (
                new freeEdgeMomentFaPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
