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
    Foam::immersedBodyMotions::quadraticBend

Description
    Deformation of the immersed body by a bending displacement that grows
    quadratically along the y direction:

    \f[
        \vec{x}(t) = \vec{x}_0 + A \xi^2 \sin(\omega t + \phi) \vec{d}, \qquad
        \xi = \frac{y_0 - y_{min}}{y_{max} - y_{min}} \in [0, 1]
    \f]

    where \f$ \vec{x}_0 \f$ is the reference position of a surface point,
    \f$ \omega = 2 \pi/T \f$ and \f$ \vec{d} \f$ the (normalised) direction.

    Usage:
    \verbatim
    motion
    {
        type        quadraticBend;
        amplitude   0.4;       // [m]
        period      1;         // [s]
        phase       0;         // [rad], optional, default 0
        yMin        0;         // [m]
        yMax        2;         // [m]
        direction   (1 0 0);   // optional, default (1 0 0)
    }
    \endverbatim

Author
    Based on the quadratic bend of the immersed boundary library contributed to
    cardiacFoam by Sairam Pamulaparthi Venkata, UCD.
    Philip Cardiff, UCD.

SourceFiles
    quadraticBend.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticBend_H
#define quadraticBend_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                   Class quadraticBend Declaration
\*---------------------------------------------------------------------------*/

class quadraticBend
:
    public immersedBodyMotion
{
    // Private Data

        //- Amplitude [m]
        const scalar amplitude_;

        //- Period [s]
        const scalar period_;

        //- Phase [rad]
        const scalar phase_;

        //- Range of y of the bending coordinate [m]
        const scalar yMin_;
        const scalar yMax_;

        //- Direction (normalised)
        const vector direction_;


public:

    //- Runtime type information
    TypeName("quadraticBend");


    // Constructors

        //- Construct from dictionary
        explicit quadraticBend(const dictionary& dict);


    //- Destructor
    virtual ~quadraticBend() = default;


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
