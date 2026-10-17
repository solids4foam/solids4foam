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
    Foam::immersedBodyMotions::valveMorph

Description
    Deformation of the immersed surface of a heart valve by the linear
    interpolation between its reference (closed) surface and an open surface
    with the same vertices:

    \f[
        X = X_0 + \alpha(t) (X_{open} - X_0), \qquad
        \alpha = E_o (1 - E_c)
    \f]

    where the edges \f$ E_o \f$ and \f$ E_c \f$ rise smoothly from 0 to 1
    across the windows \c openWindow and \c closeWindow of the phase of the
    cardiac cycle (see valveTimeLaw), so that the valve opens during the first
    and closes during the second. The velocity of a vertex is
    \f$ \dot{\alpha} (X_{open} - X_0) \f$.

    The surface \c openSurface, relative to constant/triSurface, must have
    the vertices of the valve surface in the same order, for example when both
    are written by the same program from the same mesh.

    In the original, the open surface was given by \c valveOpenFile in the
    valve dictionary, and the closing window by \c close2Window.

    Usage:
    \verbatim
    motion
    {
        type            valveMorph;
        openSurface     "valveOpen.stl";
        period          0.8;            // [s]
        openWindow      (0.05 0.15);    // fractions of the period
        closeWindow     (0.45 0.55);
    }
    \endverbatim

Author
    Based on the closed-to-open valve morph of Sairam Pamulaparthi Venkata,
    UCD.
    Philip Cardiff, UCD.

SourceFiles
    valveMorph.C

\*---------------------------------------------------------------------------*/

#ifndef valveMorph_H
#define valveMorph_H

#include "immersedBodyMotion.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                         Class valveMorph Declaration
\*---------------------------------------------------------------------------*/

class valveMorph
:
    public immersedBodyMotion
{
    // Private Data

        //- Surface of the open valve
        const fileName openSurface_;

        //- Vertices of the open valve [m]
        const pointField openPoints_;

        //- Period [s]
        const scalar period_;

        //- Opening and closing windows, as fractions of the period
        const Pair<scalar> openWindow_;
        const Pair<scalar> closeWindow_;


    // Private Member Functions

        //- Check that points0 match the vertices of the open surface
        void checkPoints(const pointField& points0) const;

        //- Interpolation factor at time t
        scalar alpha(const scalar t) const;

        //- Time derivative of the interpolation factor at time t [1/s]
        scalar alphaDerivative(const scalar t) const;


public:

    //- Runtime type information
    TypeName("valveMorph");


    // Constructors

        //- Construct from dictionary
        explicit valveMorph(const dictionary& dict);


    //- Destructor
    virtual ~valveMorph() = default;


    // Member Functions

        //- Positions at time t of the surface points
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const;

        //- The motion is not rigid
        virtual bool rigid() const
        {
            return false;
        }

        //- Not defined for a deforming motion: use pointVelocities
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const;

        //- Velocity at time t of the surface points whose reference
        //  positions are points0
        virtual tmp<vectorField> pointVelocities
        (
            const pointField& points0,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBodyMotions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
