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
    Foam::immersedBodyMotions::uniformTranslation

Description
    Translation of the immersed body at a constant velocity:

    \f[
        \vec{x}(t) = \vec{x}_0 + \vec{U} (t - t_0)
    \f]

    where \f$ \vec{U} \f$ is the velocity and \f$ t_0 \f$ the time at which
    the body is at its initial position.

    Usage:
    \verbatim
    motion
    {
        type        uniformTranslation;
        velocity    (0.1 0 0); // [m/s]
        startTime   0;         // [s], optional, default 0
    }
    \endverbatim

Author
    Philip Cardiff, UCD.

SourceFiles
    uniformTranslation.C

\*---------------------------------------------------------------------------*/

#ifndef uniformTranslation_H
#define uniformTranslation_H

#include "immersedBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBodyMotions
{

/*---------------------------------------------------------------------------*\
                   Class uniformTranslation Declaration
\*---------------------------------------------------------------------------*/

class uniformTranslation
:
    public immersedBodyMotion
{
    // Private Data

        //- Velocity [m/s]
        const vector velocity_;

        //- Time at which the body is at its initial position [s]
        const scalar startTime_;


public:

    //- Runtime type information
    TypeName("uniformTranslation");


    // Constructors

        //- Construct from dictionary
        explicit uniformTranslation(const dictionary& dict);


    //- Destructor
    virtual ~uniformTranslation() = default;


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
