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
    Foam::immersedBodyMotions::customProfileBend

Description
    Deformation of the immersed body by a two-dimensional bending of a beam
    along the y direction, with the quintic smoothstep profile
    \f$ g(\xi) = 6 \xi^5 - 15 \xi^4 + 10 \xi^3 \f$ (g(0) = g'(0) = 0,
    g(1) = 1, g'(1) = 0) of the coordinate
    \f$ \xi = (y_0 - y_{min})/H \in [0, 1] \f$, \f$ H = y_{max} - y_{min} \f$:

    \f[
        x = x_0 + A_x g(\xi) s(t), \qquad
        y = y_0 - (x_0 - x_c) \frac{A_y}{H} g'(\xi) s(t)
    \f]

    with \f$ s(t) = \sin(\omega t + \phi) \f$, \f$ \omega = 2 \pi/T \f$,
    where \f$ (x_0, y_0) \f$ is the reference position of a surface point and
    \f$ x_c \f$ the x coordinate of the beam axis: the displacement in x
    bends the axis and the displacement in y rotates the sections.

    Usage:
    \verbatim
    motion
    {
        type        customProfileBend;
        amplitudeX  0.4;       // [m]
        amplitudeY  0.8;       // [m]
        period      1;         // [s]
        phase       0;         // [rad], optional, default 0
        yMin        0;         // [m]
        yMax        2.0295;    // [m]
        xCenter     12;        // [m]
    }
    \endverbatim

Author
    Based on the custom profile bend of the immersed boundary library contributed to
    cardiacFoam by Sairam Pamulaparthi Venkata, UCD.
    Philip Cardiff, UCD.

SourceFiles
    customProfileBend.C

\*---------------------------------------------------------------------------*/

#ifndef customProfileBend_H
#define customProfileBend_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                   Class customProfileBend Declaration
\*---------------------------------------------------------------------------*/

class customProfileBend
:
    public immersedBodyMotion
{
    // Private Data

        //- Amplitudes [m]
        const scalar amplitudeX_;
        const scalar amplitudeY_;

        //- Period [s]
        const scalar period_;

        //- Phase [rad]
        const scalar phase_;

        //- Range of y of the bending coordinate [m]
        const scalar yMin_;
        const scalar yMax_;

        //- x coordinate of the beam axis [m]
        const scalar xCenter_;


    // Private Member Functions

        //- The points displaced by the profile times the factor s, if
        //  addPoints, else the displacements alone
        tmp<pointField> displacedPoints
        (
            const pointField& points0,
            const scalar s,
            const bool addPoints
        ) const;


public:

    //- Runtime type information
    TypeName("customProfileBend");


    // Constructors

        //- Construct from dictionary
        explicit customProfileBend(const dictionary& dict);


    //- Destructor
    virtual ~customProfileBend() = default;


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
