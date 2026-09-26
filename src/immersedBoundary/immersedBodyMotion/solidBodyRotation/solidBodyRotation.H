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
    Foam::immersedBodyMotions::solidBodyRotation

Description
    Rotation of the immersed body at a constant angular velocity about an
    axis through a point (the origin):

    \f[
        \vec{x}(t) = \vec{x}_c + R(\omega (t - t_0)) (\vec{x}_0 - \vec{x}_c)
    \f]

    with velocity \f$ \vec{\Omega} \times (\vec{x} - \vec{x}_c) \f$, where
    \f$ \vec{\Omega} = \omega \vec{a} \f$, \f$ \vec{a} \f$ is the unit axis
    and R the rotation about it.

    Usage:
    \verbatim
    motion
    {
        type        solidBodyRotation;
        origin      (0 0 0);   // [m]
        axis        (0 0 1);
        omega       2;         // [rad/s]
        startTime   0;         // [s], optional, default 0
    }
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    solidBodyRotation.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyRotation_H
#define solidBodyRotation_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                   Class solidBodyRotation Declaration
\*---------------------------------------------------------------------------*/

class solidBodyRotation
:
    public immersedBodyMotion
{
    // Private Data

        //- Point on the axis [m]
        const point origin_;

        //- Unit axis
        const vector axis_;

        //- Angular velocity [rad/s]
        const scalar omega_;

        //- Time at which the body is at its initial position [s]
        const scalar startTime_;


public:

    //- Runtime type information
    TypeName("solidBodyRotation");


    // Constructors

        //- Construct from dictionary
        explicit solidBodyRotation(const dictionary& dict);


    //- Destructor
    virtual ~solidBodyRotation() = default;


    // Member Functions

        //- Positions at time t of the surface points
        virtual tmp<pointField> points
        (
            const pointField& points0,
            const scalar t
        ) const;

        //- Velocity of the body at time t at the points x
        virtual tmp<vectorField> velocity
        (
            const pointField& x,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBodyMotions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
